#pragma once
#include <string>
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"

namespace Kadath::FUKA_Solvers {

/**
 * @brief Extract EOS name from the Config object
 *
 * @tparam bco_idx optional index type for binary Configs
 * @param bco optional index for binary Configs
 * @return std::string
 */
template <typename config_t, typename... bco_idx>
std::string extract_eos_name(config_t& bconfig, bco_idx... bco) {
  const std::string eos_file_abs =
      bconfig.template eos<std::string>(EOSFILE, bco...);
  const std::string eos_file = extract_filename(eos_file_abs);
  return eos_file.substr(0, eos_file.find("."));
}

/**
 * @brief Branch when maximum iterations are exceeded
 *
 * @param rank MPI rank
 * @param ite Iteration
 * @param conv Current convergence
 */
template <class solver_t>
void check_max_iter_exceeded(solver_t const& solver,
                             const int& ite,
                             const double& conv) {
  int rank = solver.get_rank();
  auto& bconfig = *solver.get_bconfig();
  bool exceeded = (ite > bconfig.seq_setting(MAX_ITER)) &&
                  conv >= bconfig.seq_setting(PREC);

  std::stringstream msg;
  if (exceeded && conv < 10. * bconfig.seq_setting(PREC)) {
    msg << "Max iterations exceeded at precison, " << conv
        << "\nFinishing since precision < 10. * PREC....\n"
        << "Running at higher resolution may help.\n";
  } else if (exceeded) {
    msg << "Max iterations exceeded at precison, " << conv << "\n";
  } else {
    return;
  }
  auto s = solver.converged_filename("termination_chkpt");
  bconfig.set_filename(s);
  if (rank == 0)
    solver.checkpoint(true);

  throw std::runtime_error(msg.str().c_str());
}

/**
 * @brief Setup the EOS operators in System_of_eqs
 *
 * @tparam eos_t C++ Polytrope/Table type
 * @param syst System of equations
 * @param p Parameter container (unused, but required by add_ope)
 */
template <class eos_t>
inline void set_eos_ope(System_of_eqs& syst, Param& p) {
  syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
  syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);
}

template <class solver_t>
inline void initialize_EOS(solver_t& solver) {
  auto& h_cut = solver.set_h_cut();
  auto& eos_file = solver.set_eos_file();
  auto& eos_type = solver.set_eos_type();
  auto& bconfig = solver.get_bconfig();

  // Initialize EOS
  h_cut =
      (*bconfig).template eos<double>(Kadath::FUKA_Config::EOS_PARAMS::HCUT);
  eos_file = (*bconfig).template eos<std::string>(
      Kadath::FUKA_Config::EOS_PARAMS::EOSFILE);
  eos_type = (*bconfig).template eos<std::string>(
      Kadath::FUKA_Config::EOS_PARAMS::EOSTYPE);

  if (eos_type == "Cold_Table") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts =
        ((*bconfig).template eos<int>(
             Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS) == 0)
            ? 2000
            : (*bconfig).template eos<int>(
                  Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS);

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut, interp_pts);
  } else if (eos_type == "Cold_PWPoly") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_PWPoly;
    EOS<eos_t, PRESSURE>::init(eos_file, h_cut);
  } else {
    throw std::invalid_argument("\nInvalid EOS Type\n)");
  }
}
}  // namespace Kadath::FUKA_Solvers