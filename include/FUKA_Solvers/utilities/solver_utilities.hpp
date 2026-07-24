#pragma once
#include <string>
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"
#include "EOS/FUKA_EOS_Utilities.hh"

using namespace Kadath::FUKA_Config;

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

template <class solver_t, typename... idx_t>
inline void initialize_EOS(solver_t& solver, idx_t... BCOidx) {
    auto& bconfig = *(solver.get_bconfig());

    ::Kadath::FUKA_EOS::EOS_initialize::init(bconfig, BCOidx...);
}

template <class bconfig_t, typename... idx_t>
inline void initialize_diffrot(bconfig_t& bconfig, idx_t... BCOidx) {
    using namespace Kadath::FUKA_Config;
    using namespace Kadath::FUKA_Solvers;
    bconfig.set_diffrot(DIFFROT_PARAMS::DIFF_LAW, BCOidx...) = "keh";
    bconfig.set_diffrot(DIFFROT_PARAMS::DIFF_ARATIO, BCOidx...) = 1e6;
    bconfig.set_diffrot(DIFFROT_PARAMS::DIFF_RRATIO, BCOidx...) = 1.0;
}
}  // namespace Kadath::FUKA_Solvers