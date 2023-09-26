#include "Solvers/solvers.hpp"
#include "mpi.h"
#include "bco_utilities.hpp"
#include <cmath>

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
using namespace ::Kadath::Margherita;

template<class eos_t, typename config_t, typename space_t>
ns_isotropic_uniform_rot_solver<eos_t, config_t, space_t>::ns_isotropic_uniform_rot_solver(config_t& config_in, 
  space_t& space_in, Scalar& nu_in, Scalar& lap_Aterm_in, Scalar& logh_in, Scalar& lap_Bterm_in, Scalar& lap_wterm_in) :
      Solver<config_t, space_t>(config_in, space_in), 
        nu(nu_in), lap_Aterm(lap_Aterm_in), logh(logh_in), lap_Bterm(lap_Bterm_in), lap_wterm(lap_wterm_in)
{ 
  lap_wterm.affect_parameters();
  lap_wterm.set_parameters()->set_m_quant() = 1 ;
  lap_wterm.std_base();

  keplerian = !std::isnan(bconfig.set(BCO_PARAMS::KEPLERIAN));
}

// standardized filename for each converged dataset at the end of each stage.
template<class eos_t, typename config_t, typename space_t>
std::string ns_isotropic_uniform_rot_solver<eos_t, config_t, space_t>::converged_filename(
  const std::string stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  const std::string eosname{extract_eos_name()};
  std::stringstream ss;
  ss << "NS_ISO";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << eosname << ".";
  
  // Add mass fixing parameter to filename
  auto default_idx = BCO_PARAMS::HC;
  if(seq) {
    update_filename_from_mass_fixing(bconfig, seq, ss);
  } else {
    auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
    ss << seq_key << "." << bconfig(default_idx) << "."; 
  }
  
  // Add spin fixing parameter to filename
  default_idx = BCO_PARAMS::CHI;
  if(seq) {
    update_filename_from_spin_fixing(bconfig, seq, ss);
  } else {
    auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
    ss << seq_key << "." << bconfig(default_idx) << "."; 
  }

  ss << bconfig(BCO_PARAMS::NSHELLS) << "."
     <<std::setfill('0') << std::setw(2) << res;
  return ss.str();
}

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_uniform_rot_solver<eos_t, config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  
  double const & final_chi = (!bconfig.control(CONTROLS::ITERATIVE_CHI)) ?
    bconfig(BCO_PARAMS::CHI) : bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI);
  double const initial_chi = bconfig(BCO_PARAMS::CHI);

  this->solver_stage = STAGES::UNIFORM_ROT;
  if(bconfig.control(CONTROLS::ITERATIVE_CHI) || keplerian) {
    exit_status = uniform_rot_stage();
    bconfig.control(CONTROLS::ITERATIVE_CHI) = false;
    stage_enabled[solver_stage] = true;
  }

  if(exit_status != RELOAD_FILE) {
    bconfig(BCO_PARAMS::CHI) = final_chi;
    exit_status = uniform_rot_stage();
  } 

  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
  return exit_status;
}

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_uniform_rot_solver<eos_t, config_t, space_t>::solve(ns_sequence const * sequence_in) {
  if(sequence_in != nullptr) {
    this->seq.reset(new ns_sequence(*sequence_in));
  }
  return this->solve();
}

template<class eos_t, typename config_t, typename space_t>
void ns_isotropic_uniform_rot_solver<eos_t, config_t, space_t>::syst_init(System_of_eqs& syst) {
  using namespace ::Kadath::Margherita;
   
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("lapAterm", lap_Aterm);
  syst.add_var("lapBterm", lap_Bterm);
  syst.add_var("wrsint", lap_wterm);  

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_def("B = (divrsint(lapBterm) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");
  syst.add_def("Brsint = multrsint(B)");
  syst.add_def("psi = log(Brsint)");
 
  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  
  // This expression does not give accurate results when compared to M_komar and the computed
  // ADM Mass from the 3D code.  The deviation from the correct answer gets large with increasing
  // differential rotation profiles, but the source of the error is unknown
  // eq. 4.21 arxiv.org/abs/1003.5015v2
  // FIXME - there must be a reason
  // syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
  // Instead, the following definition gives very accurate comparisons with M_komar
  syst.add_def(ndom - 1, "intMadm = - (dr(B)) / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");
  syst.add_def(ndom - 1, "intJ = -multrsint(multrsint(dr(w))) / 4 / 4piG");
  
  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");
 
  // define the EOS operators
  Param p;
  syst.add_ope("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope("rho", &EOS<eos_t,DENSITY>::action, &p);
 
  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

}

template<class eos_t, typename config_t, typename space_t>
void ns_isotropic_uniform_rot_solver<eos_t, config_t, space_t>::print_diagnostics(const System_of_eqs & syst, 
    const int ite, const double conv) const {

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

  // compute the ADM Angular Momentum as surface integral at infinity
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double Jadm = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

  // get the maximum and minimum coordinate radius along the surface,
  // i.e. the adapted domain boundary
  auto rs = Kadath::bco_utils::get_rmin_rmax(space, 1);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl
            << FORMAT << "Jadm: " << Jadm << endl
            << FORMAT << "CHI: " << Jadm / Madm / Madm << endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n\n";
  std::cout.flags(f);
} // end print diagnostics
/** @}*/
}}