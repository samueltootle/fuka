/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "mpi.h"
#include "bco_utilities.hpp"
#include "Solvers/co_solver_utils.hpp"
#include "bh_3d_xcts_solver.hpp"
#include "bh_3d_xcts_regrid.hpp"

/**
 * \addtogroup BH_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template<typename config_t, typename space_t>
bh_3d_xcts_solver<config_t, space_t>::bh_3d_xcts_solver(config_t& config_in, 
  space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Vector& shift_in) :
      XCTS_Solver<config_t, space_t>(config_in, space_in, base_in), 
        conf(conf_in), lapse(lapse_in), shift(shift_in), 
          fmet(Metric_flat(space_in, base_in))
{
  // initialize only the vector fields we need
  coord_vectors = default_co_vector_ary(space);

  update_fields_co(cfields, coord_vectors, {}, 0.);
}

// standardized filename for each converged dataset at the end of each stage.
template<typename config_t, typename space_t>
std::string bh_3d_xcts_solver<config_t, space_t>::converged_filename(
  const std::string stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  std::stringstream ss;
  ss << "BH";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << bconfig(MCH) << "." 
     << bconfig(CHI)<< "."
     << std::setfill('0')  << std::setw(1) << bconfig(NSHELLS) << "."
     << std::setfill('0')  << std::setw(2) << res;
  return ss.str();
}

template<typename config_t, typename space_t>
int bh_3d_xcts_solver<config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled_no_throw(MSTAGE, stage_enabled);

  // Save in case of iterative solutions
  auto const & final_chi = bconfig.seq_setting(FINAL_CHI);
  
  if(stage_enabled[TOTAL]) {
    this->solver_stage = TOTAL;
    exit_status = fixed_lapse_stage();
  }

  if(stage_enabled[TOTAL_BC] && exit_status != RELOAD_FILE) {    
    this->solver_stage = TOTAL_BC;
    if(bconfig.control(ITERATIVE_CHI)) {      
      while(bconfig.control(ITERATIVE_CHI)) {        
        
        #ifdef DEBUG
        if(rank == 0) cout << bconfig << endl;
        #endif

        exit_status = von_Neumann_stage();
        stage_enabled[STAGES::TOTAL_BC] = true;
        if(exit_status == RELOAD_FILE)
          break;
        
        if((std::abs(final_chi) <= 0.8)
          || (std::abs(bconfig(CHI)) >= 0.8)) {
          bconfig.control(ITERATIVE_CHI) = false;
        } else if(std::abs(bconfig(CHI)) < 0.8)
          bconfig(CHI) = std::copysign(0.8, final_chi);
      }
      if(exit_status != RELOAD_FILE)
        bconfig(CHI) = final_chi;
    }
    if(exit_status != RELOAD_FILE)
      exit_status = von_Neumann_stage();
  }
  
  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
  if(last_stage_idx == BIN_BOOST && exit_status != RELOAD_FILE)
    exit_status = RUN_BOOST;
  return exit_status;

}

template<typename config_t, typename space_t>
void bh_3d_xcts_solver<config_t, space_t>::syst_init(System_of_eqs& syst) {
  
  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");
 
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  syst.add_cst("PI"  , M_PI) ;
  syst.add_cst("M"   , bconfig(MIRR)) ;
  syst.add_cst("CM"  , bconfig(MCH));
  syst.add_cst("chi" , bconfig(CHI)) ;
  syst.add_var("ome" , bconfig(OMEGA));

  // include the coordinate fields
  syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("sm"  , *coord_vectors[S_BCO1]);
  
  syst.add_cst("ex"  , *coord_vectors[EX])  ;
  syst.add_cst("ey"  , *coord_vectors[EY])  ;
  syst.add_cst("einf", *coord_vectors[S_INF]);
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P"   , conf);
  syst.add_var("N"   , lapse);
  syst.add_var("bet" , shift) ;
  
  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
  syst.add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");
  
  // definitions of integrals on the excision surface
  syst.add_def("intMsq= P^4 / 16. / PI") ;

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
}

template<typename config_t, typename space_t>
void bh_3d_xcts_solver<config_t, space_t>::print_diagnostics_norot(const System_of_eqs & syst, 
    const int ite, const double conv) const {

	int ndom = space.get_nbr_domains() ;
  double r = Kadath::bco_utils::get_radius(space.get_domain(1), OUTER_BC);
  
  Val_domain integMsq(syst.give_val_def("intMsq")()(2));
  double Mirrsq = space.get_domain(2)->integ(integMsq, INNER_BC);
  double Mirr = std::sqrt(Mirrsq);
  
  // compute the ADM mass as surface integral at infinity  
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl
            << FORMAT << "Mirr: " << Mirr << std::endl;
  std::cout << FORMAT << "R: " << r << std::endl;
  std::cout.flags(f);
} // end print diagnostics norot

// runtime diagnostics specific for rotating solutions
template<typename config_t, typename space_t>
void bh_3d_xcts_solver<config_t, space_t>::print_diagnostics(System_of_eqs const & syst, 
    const int ite, const double conv) const {

  // print all the diagnostics as in the non-rotating case first  
  print_diagnostics_norot(syst, ite, conv);

	int ndom = space.get_nbr_domains() ;
  
  Val_domain integS(syst.give_val_def("intS")()(2));
  double S = space.get_domain(2)->integ(integS, INNER_BC);

  Val_domain integMsq(syst.give_val_def("intMsq")()(2));
  double Mirrsq = space.get_domain(2)->integ(integMsq, INNER_BC);
  double Mch  = std::sqrt( Mirrsq + S * S / 4. / Mirrsq );

  // output the dimensionless spin and angular frequency parameter
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << FORMAT << "Omega: " << bconfig(OMEGA) << std::endl
            << FORMAT << "S: " << S << std::endl
            << FORMAT << "Mch: " << Mch << std::endl
            << FORMAT << "Chi: " << S / Mch / Mch << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << "\n\n";
} // end print diagnostics rot
/** @}*/
}}