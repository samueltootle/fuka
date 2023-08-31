#include "Solvers/solvers.hpp"
#include "mpi.h"
#include "bco_utilities.hpp"
#include "ns_3d_xcts_regrid.hpp"
#include "Solvers/sequences/sequence_utilities.hpp"
#include <cmath>

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
using namespace ::Kadath::Margherita;

template<class eos_t, typename config_t, typename space_t>
ns_3d_xcts_solver<eos_t, config_t, space_t>::ns_3d_xcts_solver(config_t& config_in, 
  space_t& space_in, Base_tensor& base_in,  
    Scalar& conf_in, Scalar& lapse_in, Scalar& logh_in, Vector& shift_in) :
      XCTS_Solver<config_t, space_t>(config_in, space_in, base_in), 
        conf(conf_in), lapse(lapse_in), logh(logh_in), shift(shift_in), 
          fmet(Metric_flat(space_in, base_in))
{
  // initialize only the vector fields we need
  coord_vectors = default_co_vector_ary(space);

  update_fields_co(cfields, coord_vectors,{}, 0.);
}

// standardized filename for each converged dataset at the end of each stage.
template<class eos_t, typename config_t, typename space_t>
std::string ns_3d_xcts_solver<eos_t, config_t, space_t>::converged_filename(
  const std::string stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  const std::string eosname{extract_eos_name()};
  std::stringstream ss;
  ss << "NS";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << eosname << ".";
  
  // Add mass fixing parameter to filename
  if(seq && seq->is_set()) {
    auto idx{std::get<0>(seq->get_indices())};
    auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, idx);

    switch(idx) {
      case BCO_PARAMS::HC:
      case BCO_PARAMS::NC:
      case BCO_PARAMS::MADM:
      case BCO_PARAMS::MB:      
        ss << seq_key << "." << bconfig(idx);
        break;
      default:
        ss << "madm." << bconfig(BCO_PARAMS::MADM) << ".";
    }
  } else {
    ss << "madm." << bconfig(BCO_PARAMS::MADM) << "."; 
  }
  
  // Add spin fixing parameter to filename
  if(stage != "NOROT_BC") {
    if(seq && seq->is_set()) {
      auto idx{std::get<0>(seq->get_indices())};
      auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, idx);
      switch(idx) {
        case BCO_PARAMS::OMEGA:
        case BCO_PARAMS::JADM:
        case BCO_PARAMS::CHI:
          ss << seq_key << "." << bconfig(idx);
          break;
        default:
          ss << "chi." << bconfig(BCO_PARAMS::CHI) << "."; 
          break;
        break;
      }
    }    
    else {
    ss << "chi." << bconfig(BCO_PARAMS::CHI) << "."; 
  }
  }
  else ss << "0.";
  ss << bconfig(BCO_PARAMS::NSHELLS) << "."
     <<std::setfill('0') << std::setw(2) << res;
  return ss.str();
}

template<class eos_t, typename config_t, typename space_t>
int ns_3d_xcts_solver<eos_t, config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled_no_throw(MSTAGE, stage_enabled);
  if(rank == 0) std::cout << "Last stage: " << last_stage << "\n";

  double const & final_chi = (!bconfig.control(CONTROLS::ITERATIVE_CHI)) ?
    bconfig(BCO_PARAMS::CHI) : bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI);
  double const initial_chi = bconfig(BCO_PARAMS::CHI);
  
  if(stage_enabled[STAGES::NOROT_BC]) {
    bconfig(BCO_PARAMS::CHI) = 0.;
    this->solver_stage = STAGES::NOROT_BC;
    exit_status = norot_stage(false);
    bconfig(BCO_PARAMS::CHI) = initial_chi;
  }

  if(stage_enabled[STAGES::TOTAL_BC] &&  exit_status != RELOAD_FILE) {
    this->solver_stage = STAGES::TOTAL_BC;
    if(bconfig.control(CONTROLS::ITERATIVE_CHI)) {
      exit_status = uniform_rot_stage();
      bconfig.control(CONTROLS::ITERATIVE_CHI) = false;
      stage_enabled[STAGES::TOTAL_BC] = true;
    }

    if(exit_status != RELOAD_FILE) {
      bconfig(BCO_PARAMS::CHI) = final_chi;
      exit_status = uniform_rot_stage();
    }
  }
  
  if(stage_enabled[STAGES::TESTING] &&  exit_status != RELOAD_FILE) {
    this->solver_stage = STAGES::TESTING;
    exit_status = differential_rot_stage();
  }

  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
  return exit_status;
}

template<class eos_t, typename config_t, typename space_t>
int ns_3d_xcts_solver<eos_t, config_t, space_t>::solve(Parameter_sequence<BCO_PARAMS> const * sequence_in) {
  if(sequence_in != nullptr) {
    this->seq.reset(new Parameter_sequence<BCO_PARAMS>(*sequence_in));
  }
  return this->solve();
}

template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::syst_init(System_of_eqs& syst) {
  using namespace ::Kadath::Margherita;
  
  auto& space = syst.get_space();
  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");
 
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  
  // include the coordinate fields
  syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("sm"  , *coord_vectors[S_BCO1]);
  syst.add_cst("einf", *coord_vectors[S_INF]);
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P"   , conf);
  syst.add_var("N"   , lapse);
  
  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
 
  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
  syst.add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
  
  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");
 
  // define the EOS operators
  Param p;
  syst.add_ope("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope("rho", &EOS<eos_t,DENSITY>::action, &p);
  syst.add_ope("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);
 
  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");
  syst.add_def("dHdlnrho = dHdlnrho(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

}

template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::print_diagnostics_norot(const System_of_eqs & syst, 
    const int ite, const double conv) const {

  // total number of domains	
  int ndom = space.get_nbr_domains() ;

  // compute the baryonic mass at volume integral from the given integrant
  double baryonic_mass =
      syst.give_val_def("intMb")()(0).integ_volume() +
      syst.give_val_def("intMb")()(1).integ_volume();

  // compute the ADM mass as surface integral at infinity  
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // get the maximum and minimum coordinate radius along the surface,
  // i.e. the adapted domain boundary
  auto rs = Kadath::bco_utils::get_rmin_rmax(space, 1);

  // alternative, equivalent ADM mass integral
  Val_domain integMadmalt(syst.give_val_def("intMadmalt")()(ndom - 1));
  double Madmalt = space.get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Madm_ql: " << Madmalt 
            << " [" << std::abs(Madm - Madmalt) / Madm << "]" << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n";
  std::cout.flags(f);
} // end print diagnostics norot

// runtime diagnostics specific for rotating solutions
template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::print_diagnostics(System_of_eqs const & syst, 
    const int ite, const double conv) const {

  // print all the diagnostics as in the non-rotating case first  
  print_diagnostics_norot(syst, ite, conv);

  // total number of domains
	int ndom = space.get_nbr_domains() ;

  // compute the ADM angular momentum as surface integral at infinity  
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

  // compute the ADM mass as surface integral at infinity    
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // output the dimensionless spin and angular frequency parameter
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << FORMAT << "Jadm: " << J << std::endl
            << FORMAT << "Chi: " << J / Madm / Madm << " [" << bconfig(CHI) << "]\n"
            << FORMAT << "Omega: " << bconfig(OMEGA) << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << "\n\n";
} // end print diagnostics rot

template<class eos_t, typename config_t, typename space_t>
void ns_3d_xcts_solver<eos_t, config_t, space_t>::update_config_quantities(const double loghc) {
  bconfig.set(HC) = std::exp(loghc);
  bconfig.set(NC) = EOS<eos_t,DENSITY>::get(bconfig(HC));
}
/** @}*/
}}
