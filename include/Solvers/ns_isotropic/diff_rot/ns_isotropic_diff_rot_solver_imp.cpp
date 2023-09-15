#include "Solvers/solvers.hpp"
#include "mpi.h"
#include "bco_utilities.hpp"
#include "name_tools.hpp"
#include <cmath>

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
using namespace ::Kadath::Margherita;

template<class eos_t, typename config_t, typename space_t>
ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::ns_isotropic_diff_rot_solver(config_t& config_in, 
  space_t& space_in, Scalar& nu_in, Scalar& lap_Aterm_in, Scalar& logh_in, Scalar& lap_Bterm_in, Scalar& lap_wterm_in, Scalar& Omega_in) :
      Solver<config_t, space_t>(config_in, space_in), 
        nu(nu_in), lap_Aterm(lap_Aterm_in), logh(logh_in), lap_Bterm(lap_Bterm_in), lap_wterm(lap_wterm_in), Omega(Omega_in)
{ 
  lap_wterm.affect_parameters();
  lap_wterm.set_parameters()->set_m_quant() = 1 ;
  lap_wterm.std_base();

  law = str_tolower(bconfig.template diffrot<std::string>(DIFFROT_PARAMS::DIFF_LAW));
}

// standardized filename for each converged dataset at the end of each stage.
template<class eos_t, typename config_t, typename space_t>
std::string ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::converged_filename(
  const std::string stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  const std::string eosname{extract_eos_name()};
  std::stringstream ss;
  ss << "NS_ISO"
     << "_" << stage << "." 
     << eosname << "."
     << law << ".";
  
  // Add mass fixing parameter to filename
  auto default_idx = BCO_PARAMS::HC;
  if(seq) {
    update_filename_from_mass_fixing(bconfig, seq, ss);
  } else {
    auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
    ss << seq_key << "." << bconfig(default_idx) << "."; 
  }
  ss << bconfig(BCO_PARAMS::OMEGA)<< ".";
  ss << bconfig(BCO_PARAMS::NSHELLS) << "."
     << std::setfill('0') << std::setw(2) << res;
  return ss.str();
}

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();

  this->solver_stage = STAGES::DIFF_ROT;
  std::string law = str_tolower(bconfig.template diffrot<std::string>(DIFFROT_PARAMS::DIFF_LAW));
  if(law == "keh")
    exit_status = keh_stage();  

  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);
  return exit_status;
}

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::solve(ns_sequence const * sequence_in) {
  if(sequence_in != nullptr) {
    this->seq.reset(new ns_sequence(*sequence_in));
  }
  return this->solve();
}

template<class eos_t, typename config_t, typename space_t>
void ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::syst_init(System_of_eqs& syst) {
  using namespace ::Kadath::Margherita;
   
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  
  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("lapAterm", lap_Aterm);
  syst.add_var("lapBterm", lap_Bterm);
  syst.add_var("wrsint", lap_wterm);
  syst.add_var("ome", Omega);

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_def("B = (divrsint(lapBterm) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");
  syst.add_def("Fomega = B^2 * multrsint(multrsint(ome - w)) "
                      "/ (N^2 - multrsint(B * (ome - w))^2)");
 
  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
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

  syst.add_def("U = multrsint(B / N * (ome - w))");
  syst.add_def("Usq = U*U");
  syst.add_def("Wsq = 1 / (1 - Usq)");
  syst.add_def("W = sqrt(Wsq)");
}

template<class eos_t, typename config_t, typename space_t>
void ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::print_diagnostics(const System_of_eqs & syst, 
    const int ite, const double conv) const {

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
  auto rs = Kadath::bco_utils::get_rmin_rmax(space, 1);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            // << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n\n";
  std::cout.flags(f);
} // end print diagnostics
/** @}*/
}}