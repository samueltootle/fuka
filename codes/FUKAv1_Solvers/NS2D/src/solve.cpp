/*
 * Copyright 2021
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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
// FUKA includes
#include "Solvers/solver_startup.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"

#include "kadath_adapted.hpp"
#include "kadath_adapted_polar.hpp"

#include "mpi.h"
#include <sstream>
#define O(X) std::cout << "var " #X ": " << X << std::endl;
using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Solvers;

// // forward declarations
template<typename space_t, typename syst_t, typename config_t>
void print_diagnostics_norot(space_t const & space, syst_t const & syst, const config_t & bconfig, int iter, double err);
template<class eos_t, class config_t>
void update_config(config_t& bconfig, Scalar& logh) ;
// template<typename space_t, typename syst_t, typename config_t>
// void print_diagnostics_rot(space_t const & space, syst_t const & syst, const config_t & bconfig, int const ite, double const conv);

template<typename config_t>
std::string converged_filename(const std::string&  stage, config_t bconfig);

template<class eos_t, typename config_t>
int driver(config_t& bconfig, std::string outputdir);

template<class eos_t, typename config_t>
int NS_solver_2d_norot (config_t& bconfig, bool fixed = false);

template<class eos_t, typename config_t>
int NS_solver_2d_uniform_rot (config_t& bconfig);

template<class eos_t, typename config_t>
int NS_solver_2d_differential_rot (config_t& bconfig);
// end forward declarations

int main(int argc, char **argv) {
  int rc = MPI_Init(&argc, &argv) ;
  if (rc!=MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl ;
    MPI_Abort(MPI_COMM_WORLD, rc) ;
  }
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

  using config_t = kadath_config_boost<BCO_NS_INFO>;
  using InitSolver = Initialize_Solver<config_t>;
  
  // Initialize static member variables
  InitSolver::input_configname = "initial_ns.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;

  // specify output directory  
  std::string outputdir = InitSolver::outputdir;

  if(!InitSolver::example_setup) {
    // We now have to assume bconfig is a minimal config
    // that contains sequences _init/_final
    bconfig.open_config();
  
    // load and setup the EOS
    const double h_cut = bconfig.eos<double>(HCUT);
    const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      driver<eos_t>(bconfig, outputdir);
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      driver<eos_t>(bconfig, outputdir);
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
  }
  #ifdef ENABLE_GPU_USE
    if(rank==0)
    {
      TESTING_CHECK(magma_finalize());
    }
  #endif
    MPI_Finalize();

    return EXIT_SUCCESS;
  } // end main()

template<class eos_t, typename config_t>
int driver(config_t& bconfig, std::string outputdir) {
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
  int exit_status = EXIT_SUCCESS;

  bconfig.set_outputdir(outputdir);

  // specify which stages are enabled, pulled from the configuration file
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  
  // last stage of the enabled stages
  auto [ last_stage, last_stage_idx ] 
    = get_last_enabled(MSTAGE, stage_enabled);

  if(rank == 0){
      std::cout << "Last Stage Enabled: " << last_stage << std::endl;
  }

  if(stage_enabled[STAGES::NOROT_BC]) {
    if(bconfig.control(CONTROLS::SEQUENCES))
      exit_status = NS_solver_2d_norot<eos_t>(bconfig, true);
    exit_status = NS_solver_2d_norot<eos_t>(bconfig);
  }
  if(stage_enabled[STAGES::TOTAL_BC])
    exit_status = NS_solver_2d_uniform_rot<eos_t>(bconfig);
  if(stage_enabled[STAGES::TESTING])
    exit_status = NS_solver_2d_differential_rot<eos_t>(bconfig);
  return exit_status;
}

template<class eos_t, typename config_t>
int NS_solver_2d_norot (config_t& bconfig, bool fixed) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #ifdef ENABLE_GPU_USE
    if(rank==0)
	{
		TESTING_CHECK(magma_init());
		magma_print_environment();
	}
#endif

  // convergence threshold
  // FIXME should this be part of the config?
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy, a variable in the system of equations 
  double loghc = std::log(bconfig(HC));
  std::string stage_name{"NOROT"};
  if(fixed) stage_name += "_FIXED";

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // load the fields defined on the space
	Scalar nulogA   (space, ff1) ;
	Scalar nu  (space, ff1) ;
  Scalar logh   (space, ff1) ;
	fclose(ff1) ;

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

  // std::cout << logh << std::endl;
  Scalar level(space);
  for (int d = 0; d < ndom - 1; d++)
    level.set_domain(d) = space.get_domain(d)->get_radius() * space.get_domain(d)->get_radius() 
                         -bconfig(RMID) * bconfig(RMID);
  level.set_domain(ndom - 1) = 1;
  level.std_base();
  
  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);

  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
//       syst.add_cst("Mb"  , bconfig(MB));
  
  syst.add_cst("Hc", loghc);
  syst.add_var("Madm", bconfig(MADM));
  

  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("nulogA", nulogA);

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(nulogA - nu)");
  syst.add_cst("lev", level); 

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = -dr(A) / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      // sources
      syst.add_def(d, "E = press * (1 + eps)");
      syst.add_def(d, "S = delta * 3 * press");
      syst.add_def(d, "Spp = press * delta");
      
      // constraint equations
      syst.add_def(d, "eqnu = delta * ( lap(nu) + scal(grad(nu), grad(nulogA)) ) - 4piG * A^2 * (E + S)") ;
      syst.add_def(d, "eqnulogA = delta * ( lap2(nulogA) + scal(grad(nu), grad(nu)) ) - 2 * 4piG * A^2 * Spp") ;
      // Extra...
      // syst.add_def(d, "eqNA = dr(drNA) + 3 * divr(drNA) - 4 * 4piG * NA * A^2 * press") ;


      // // definition for the baryonic mass integral
      // syst.add_def(d, "intMb = P^6 * rho");
      
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, "firstint = H + log(N)");

      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(nulogA))") ;
      syst.add_def(d, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))") ;
      break;
    }
  }
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqnulogA=0", "nulogA", "dn(nulogA)");
  space.add_eq_int_inf(syst, "integ(intMadm) - Madm = 0");
  
  
  if(fixed)
  syst.add_eq_bc(1, OUTER_BC, "lev=0");
  else
  syst.add_eq_bc(1, OUTER_BC, "H=0");
  
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");

  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nulogA=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // Kadath::Solver test_solver(
  // Verbosity = 1, 
  // Tolerance = conv_thres,
  // MaxNbIter = 15);

  // repeat until convergence is achieved
  while (!endloop) {        
    // do exactly one newton step
    endloop = syst.do_newton(conv_thres, conv);

    // output the data and diagnostics at this particular step
    std::stringstream ss;
    ss << "rot_3d_testing" << ite - 1 ;
    bconfig.set_filename(ss.str());

    if (rank == 0) {
      print_diagnostics_norot(space, syst, bconfig, ite, conv);
      // if(bconfig.control(CHECKPOINT))
        // bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
    // count the steps
    ite++;
  }
  // if(!fixed)  
  update_config<eos_t>(bconfig, logh);
  bconfig.set_filename(converged_filename(stage_name, bconfig));
  bconfig.control(CONTROLS::SEQUENCES) = false;
  Scalar bet(space);  
    {  
    Scalar B(exp(nulogA - nu));
    Scalar N(exp(nu));
    Scalar tmp(N * B - 1);
    bet = Scalar(tmp.mult_sin_theta().mult_r());
    }
  Scalar wrsint(space);
  wrsint.annule_hard();
  wrsint.std_base();
  if(rank == 0)
    bco_utils::save_to_file(space, bconfig, nulogA, nu, logh, bet, wrsint);
  MPI_Barrier(MPI_COMM_WORLD);
  return EXIT_SUCCESS;
}

template<class eos_t, typename config_t>
int NS_solver_2d_uniform_rot (config_t& bconfig) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #ifdef ENABLE_GPU_USE
  if(rank==0)
	{
		TESTING_CHECK(magma_init());
		magma_print_environment();
	}
  #endif

  // convergence threshold
  // FIXME should this be part of the config?
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy, a variable in the system of equations 
  double loghc = std::log(bconfig(HC));
  std::string stage_name{"ROT"};

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

  Scalar bet(space);
  Scalar wrsint(space);

  // load the fields defined on the space
	Scalar nulogA   (space, ff1) ;
	Scalar nu  (space, ff1) ;
  Scalar logh   (space, ff1) ;
  
  if(bconfig.set_field(BCO_FIELDS::SHIFT) && bconfig.set_field(BCO_FIELDS::NP)) {
    if(rank == 0)
    std::cout << "Starting from shift and omega fields from file\n";
    bet = Scalar(space, ff1);
    wrsint = Scalar(space, ff1);
  } else {
    if(rank == 0)
      std::cout << "Generating new shift and omega fields\n";
    Scalar w(space);
    w.annule_hard();
    w.std_base();
    wrsint = Scalar(w.mult_r().mult_sin_theta());
    
    Scalar B(exp(nulogA - nu));
    Scalar N(exp(nu));
    Scalar tmp(N * B - 1);
    bet = Scalar(tmp.mult_sin_theta().mult_r());
  }
  fclose(ff1) ;
  wrsint.affect_parameters();
  wrsint.set_parameters()->set_m_quant() = 1 ;
  wrsint.std_base();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);

  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  syst.add_cst("Omega", bconfig(BCO_PARAMS::OMEGA));

  // syst.add_cst("Mb"  , bconfig(MB));
  syst.add_var("Madm", bconfig(MADM));

  syst.add_cst("Hc", loghc);

  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("nulogA", nulogA);
  syst.add_var("bet", bet);
  syst.add_var("wrsint", wrsint);  

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(nulogA - nu)");
  syst.add_def("B = (divrsint(bet) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom-1; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      syst.add_def(d, "U = multrsint(B / N * (Omega - w))");
      syst.add_def(d, "Usq = U*U");
      syst.add_def(d, "Wsq = 1 / 1 - Usq");
      
      // sources
      syst.add_def(d, "edens = rho * (1 + eps)");
      syst.add_def(d, "E = Wsq * press * h - press * delta");
      syst.add_def(d, "Srrtt = press * delta");
      syst.add_def(d, "pressp = delta * multrsint(B * (E + Srrtt) * U)");
      syst.add_def(d, "Spp = delta * press * (1 + Usq) + E * Usq");
      syst.add_def(d, "S = 2 * Srrtt + Spp");
      
      syst.add_def(d, "eqnu  = delta * lap(nu) + delta * scal(grad(nu), grad(nu + log(B))) "
                            "- delta * multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w)) "
                            "- 4piG * A^2 * (E + S)");
      syst.add_def(d, "eqnulogA = delta * lap2(nulogA) + delta * scal(grad(nu), grad(nu))"
                      "- 3 * delta * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))"
                      "- 2 * 4piG * A^2 * Spp");
      syst.add_def(d, "eqbet = delta * lap2(bet) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
      syst.add_def(d, "eqw = delta * lap(wrsint) - delta * multrsint(scal(grad(w), grad(nu - 3 * log(B))))"
                          "+ 4 * 4piG * N * A^2 / B^2 * divrsint(pressp)");
      
      // definition for the baryonic mass integral - need a volume integral first
      // syst.add_def(d, "intMb = P^6 * rho");
      
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, "firstint = H + log(N) - 0.5 * log(Wsq)");
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");
      syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                      "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))"
                "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqbet = lap2(bet)");
      syst.add_def(d, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
      break;
    }
  }
    syst.add_eq_full(ndom-1, "H = 0");
    syst.add_def(ndom-1, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                            "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(multr(grad(w)), multr(grad(w))) ");
    syst.add_def(ndom-1, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))"
              "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(multr(grad(w)), multr(grad(w)))");
    syst.add_def(ndom-1, "eqbet = lap2(bet)");
    syst.add_def(ndom-1, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
  
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqnulogA=0", "nulogA", "dn(nulogA)");
  space.add_eq(syst, "eqbet=0", "bet", "dn(bet)");
  space.add_eq(syst, "eqw=0", "wrsint", "dn(wrsint)");
  
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");
  syst.add_eq_bc(1, OUTER_BC, "H=0");  

  space.add_eq_int_inf(syst, "integ(intMadm) - Madm = 0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nulogA=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "bet=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "wrsint=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // repeat until convergence is achieved
  while (!endloop) {        
    // do exactly one newton step
    endloop = syst.do_newton(conv_thres, conv);

    // output the data and diagnostics at this particular step
    std::stringstream ss;
    ss << "rot_2k_chkpt" << ite - 1 ;
    bconfig.set_filename(ss.str());

    if (rank == 0) {
      print_diagnostics_norot(space, syst, bconfig, ite, conv);
      // if(bconfig.control(CHECKPOINT))
        // bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
    // count the steps
    ite++;
  }
  update_config<eos_t>(bconfig, logh);
    // std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
    // stage_enabled.fill(false);
    // stage_enabled[STAGES::TOTAL_BC] = true;
    bconfig.set_filename(converged_filename(stage_name, bconfig));
    bconfig.set_field(BCO_FIELDS::SHIFT) = true;
    bconfig.set_field(BCO_FIELDS::NP) = true;
  if(rank == 0) {
    std::cout << "Success!\n";
    bco_utils::save_to_file(space, bconfig, nulogA, nu, logh, bet, wrsint);
  }  
  MPI_Barrier(MPI_COMM_WORLD);
  
  return EXIT_SUCCESS;
}

template<class eos_t, typename config_t>
int NS_solver_2d_differential_rot (config_t& bconfig) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #ifdef ENABLE_GPU_USE
  if(rank==0)
	{
		TESTING_CHECK(magma_init());
		magma_print_environment();
	}
  #endif

  // convergence threshold
  // FIXME should this be part of the config?
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy, a variable in the system of equations 
  double loghc = std::log(bconfig(HC));
  std::string stage_name{"DIFFROT"};

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

  Scalar bet(space);
  Scalar wrsint(space);
  Scalar Omega(space);
  Omega = bconfig(BCO_PARAMS::OMEGA);
  // Omega.annule_hard();
  // Omega.set_domain(0) = bconfig(BCO_PARAMS::OMEGA);
  // Omega.set_domain(1) = bconfig(BCO_PARAMS::OMEGA);
  Omega.std_base();

  // load the fields defined on the space
	Scalar nulogA   (space, ff1) ;
	Scalar nu  (space, ff1) ;
  Scalar logh   (space, ff1) ;
  
  if(bconfig.set_field(BCO_FIELDS::SHIFT) && bconfig.set_field(BCO_FIELDS::NP)) {
    if(rank == 0)
      std::cout << "Starting from shift and omega fields from file\n";
    bet = Scalar(space, ff1);
    wrsint = Scalar(space, ff1);
    if(bconfig.set_field(BCO_FIELDS::PHI)) {
      if(rank == 0)
        std::cout << "Starting from Omega field from file\n";
      Omega = Scalar(space, ff1);
    }
  } else {
    if(rank == 0)
      std::cout << "Generating new shift and omega fields\n";
    Scalar w(space);
    w.annule_hard();
    w.std_base();
    wrsint = Scalar(w.mult_r().mult_sin_theta());
    
    Scalar B(exp(nulogA - nu));
    Scalar N(exp(nu));
    Scalar tmp(N * B - 1);
    bet = Scalar(tmp.mult_sin_theta().mult_r());
  }
  fclose(ff1) ;
  wrsint.affect_parameters();
  wrsint.set_parameters()->set_m_quant() = 1 ;
  wrsint.std_base();

  auto npts = space.get_domain(1)->get_nbr_points();
  Index pos_eq (npts);
  pos_eq.set(0) = npts(0) - 1; /// Set to outer radius
  pos_eq.set(1) = npts(1) - 1; /// Set theta to be on the xy plane.

  Index pos_pole (npts);
  pos_pole.set(0) = npts(0) - 1; /// Set to outer radius

  auto adpt_dom = space.get_domain(1);
  double R0 = adpt_dom->get_radius()(pos_eq);
  double Rp = adpt_dom->get_radius()(pos_pole);

  // Differential rotation fixing parameters
  double diffA = 5.9;
  double Rratio = Rp / R0;
  int q = 1;

  std::string jint{};
  std::string jome{"diffA^2 * Omega * (omeratio^"+std::to_string(q)+" - 1)"};
  std::string F{"F = " + jome};
  // std::string eqOme{"eqOme = P^4 * Wsquare * f_ij * U^i * mg^j / N"};
  // std::string eqOme{"eqOme = (P^4 * Wsq * U * mg^j / N) / A^2 - " + jome};
  switch(q) {
    case 2:
      jint = "diffA^2 * Omega^2 * (omeratio^2 * log(Omega) - 0.5)"; // - omec^2 * (log(omec) - 0.5)";
      break;
    default:
      jint = "diffA^2 * Omega^2 * ((1 / (2-q)) * omeratio^"+std::to_string(q)+" - 0.5)"; // - omec^2 * q / (4 - 2 * q)";
      break;
  }
  std::string firstint{"firstint = (H + log(N) - 0.5 * log(Wsq)) + " + jint};
  if (rank == 0)
    std::cout << "###################################" << std::endl
              << "Differential Rotating models"      << std::endl
              << "Law: " << F << std::endl
              << firstint << std::endl
              // << eqOme << std::endl
              << "q: " << q << std::endl
              << "A: " << diffA << std::endl
              << "Rp/Re: " << Rratio << std::endl
              << "R0: " << R0 <<std::endl
              << "###################################" << std::endl;

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);

  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));

  // syst.add_cst("Mb"  , bconfig(MB));
  syst.add_var("Madm", bconfig(MADM));
  syst.add_cst("Hc", loghc);
  syst.add_cst("omec", bconfig(BCO_PARAMS::OMEGA));

  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("nulogA", nulogA);
  syst.add_var("bet", bet);
  syst.add_var("wrsint", wrsint);
  syst.add_var("Omega", Omega);

  syst.add_cst("q", q);
  syst.add_cst("diffA", diffA);
  syst.add_def("omeratio = omec / Omega");
  syst.add_def(F.c_str());

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(nulogA - nu)");
  syst.add_def("B = (divrsint(bet) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");
  syst.add_def("Fomega = B^2 * multrsint(multrsint(Omega - w)) "
                      "/ (N^2 - multrsint(B * (Omega - w))^2)");

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom-1; d++) {
    syst.add_eq_full(d, "Fomega - F = 0");
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      syst.add_def(d, "U = multrsint(B / N * (Omega - w))");
      syst.add_def(d, "Usq = U*U");
      syst.add_def(d, "Wsq = 1 / 1 - Usq");
      
      // sources
      syst.add_def(d, "edens = rho * (1 + eps)");
      syst.add_def(d, "E = Wsq * (edens + press) - press");
      syst.add_def(d, "Srrtt = press");
      syst.add_def(d, "pressp = multrsint(B * (E + Srrtt) * U)");
      syst.add_def(d, "Spp = press * (1 + Usq) + E * Usq");
      syst.add_def(d, "S = 2 * Srrtt + Spp");
      
      syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                            "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w)) "
                            "- 4piG * A^2 * (E + S)");
      syst.add_def(d, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))"
                      "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))"
                      "- 2 * 4piG * A^2 * Spp");
      syst.add_def(d, "eqbet = lap2(bet) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
      syst.add_def(d, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))"
                          "+ 4 * 4piG * N * A^2 / B^2 * divrsint(pressp)");
      
      
      // definition for the baryonic mass integral - need a volume integral first
      // syst.add_def(d, "intMb = P^6 * rho");
      
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, firstint.c_str());
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                      "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))"
                "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqbet = lap2(bet)");
      syst.add_def(d, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
      break;
    }
  }
    syst.add_eq_full(ndom-1, "Fomega - F = 0");
    syst.add_eq_full(ndom-1, "H = 0");

    syst.add_def(ndom-1, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                            "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(multr(grad(w)), multr(grad(w))) ");
    syst.add_def(ndom-1, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))"
              "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(multr(grad(w)), multr(grad(w)))");
    syst.add_def(ndom-1, "eqbet = lap2(bet)");
    syst.add_def(ndom-1, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
  
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqnulogA=0", "nulogA", "dn(nulogA)");
  space.add_eq(syst, "eqbet=0", "bet", "dn(bet)");
  space.add_eq(syst, "eqw=0", "wrsint", "dn(wrsint)");
  
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");
  syst.add_eq_bc(1, OUTER_BC, "H=0");  

  space.add_eq_int_inf(syst, "integ(intMadm) - Madm = 0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nulogA=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "bet=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "wrsint=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // repeat until convergence is achieved
  while (!endloop) {        
    // do exactly one newton step
    endloop = syst.do_newton(conv_thres, conv);

    // output the data and diagnostics at this particular step
    std::stringstream ss;
    ss << "rot_2k_chkpt" << ite - 1 ;
    bconfig.set_filename(ss.str());

    if (rank == 0) {
      print_diagnostics_norot(space, syst, bconfig, ite, conv);
      // if(bconfig.control(CHECKPOINT))
        // bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
    // count the steps
    ite++;
  }
      update_config<eos_t>(bconfig, logh);
    std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
    stage_enabled.fill(false);
    stage_enabled[STAGES::TESTING] = true;
    bconfig.set_filename(converged_filename(stage_name, bconfig));
    bconfig.set_field(BCO_FIELDS::SHIFT) = true;
    bconfig.set_field(BCO_FIELDS::NP) = true;
    bconfig.set_field(BCO_FIELDS::PHI) = true;
  if(rank == 0) {
    std::cout << "Success!\n";
    bco_utils::save_to_file(space, bconfig, nulogA, nu, logh, bet, wrsint, Omega);
  }  
  MPI_Barrier(MPI_COMM_WORLD);
  return EXIT_SUCCESS;
}

template<class eos_t, class config_t>
void update_config(config_t& bconfig, Scalar& logh) {
  auto const &  space = logh.get_space();

  auto rs = bco_utils::get_rmin_rmax(space, 1);
  auto loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);

  bconfig.set(BCO_PARAMS::RMID) = rs[0];
  bconfig.set(BCO_PARAMS::HC) = std::exp(loghc);
  bconfig.set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig(BCO_PARAMS::HC));
}

// diagnostics at runtime
#define FORMAT std::setw(10) << std::left << std::showpos 
template<typename space_t, typename syst_t, typename config_t>
void print_diagnostics_norot(space_t const & space, syst_t const & syst, 
    const config_t & bconfig, int iter, double err) {

  // total number of domains	
  int ndom = space.get_nbr_domains() ;

  // compute the baryonic mass at volume integral from the given integrant
  // double baryonic_mass =
  //     syst.give_val_def("intMb")()(0).integ_volume() +
  //     syst.give_val_def("intMb")()(1).integ_volume();

  // compute the ADM mass as surface integral at infinity  
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // get the maximum and minimum coordinate radius along the surface,
  // i.e. the adapted domain boundary
  auto rs = bco_utils::get_rmin_rmax(space, 1);

  // alternative, equivalent ADM mass integral
  // Val_domain integMadmalt(syst.give_val_def("intMadmalt")()(ndom - 1));
  // double Madmalt = space.get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << iter << std::endl
            << FORMAT << "Error: " << err << std::endl
            // << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            // << FORMAT << "Madm_ql: " << Madmalt 
            // << " [" << std::abs(Madm - Madmalt) / Madm << "]" << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n";
  std::cout.flags(f);
} // end print diagnostics norot


//standardized filename for each converged dataset at the end of each stage.
template<typename config_t>
std::string converged_filename(const std::string&  stage, config_t bconfig) {
  std::stringstream ss;
  ss << "NS2D";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << bconfig(MADM) << "." 
     << bconfig(CHI)<< "."   
     << std::setfill('0')  << std::setw(2) << bconfig(BCO_RES);
  return ss.str();
}
