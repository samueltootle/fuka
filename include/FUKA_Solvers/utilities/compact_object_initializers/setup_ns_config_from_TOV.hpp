#pragma once
#include <cmath>
#include <memory>
#include "Configurator/config_bco.hpp"
#include "EOS/standalone/tov.hh"
#include "bco_utilities.hpp"
using namespace Kadath::FUKA_Config;

namespace Kadath::FUKA_Solvers {

/**
 * @brief Set up the configuration for a neutron star based on a 1D TOV solution
 * @tparam eos_t EOS wrapper type
 * @tparam config_t configuration file type
 * @param bconfig Configuration object
 * @param mass_fixing_idx Config index used to fix stellar mass
 * @return Unique pointer to the TOV solution
 */
template <typename eos_t, typename config_t>
auto setup_ns_config_from_TOV(config_t& bconfig, size_t mass_fixing_idx) {
    using namespace Kadath::Margherita;
    auto tov = std::make_unique<MargheritaTOV<eos_t>>();

    bool adm_mass_fixing = (mass_fixing_idx == BCO_PARAMS::MADM);
    bool nc_mass_fixing = (mass_fixing_idx == BCO_PARAMS::NC);
    bool hc_mass_fixing = (mass_fixing_idx == BCO_PARAMS::HC);

    if (adm_mass_fixing) {
        bool use_Mmax = tov->solve_for_MADM(bconfig(mass_fixing_idx));

        if (use_Mmax)
            bconfig.set(BCO_PARAMS::MADM) = tov->mass;
    } else if (nc_mass_fixing) {
        tov->solve(bconfig(mass_fixing_idx));
    } else if (hc_mass_fixing) {
        bconfig.set(BCO_PARAMS::NC) =
            EOS<eos_t, DENSITY>::get(bconfig(BCO_PARAMS::HC));
        tov->solve(bconfig(BCO_PARAMS::NC));
    } else {
        std::string msg{
            "1D TOV solver not implemented for sequences using idx " +
            std::to_string(mass_fixing_idx)};
        throw std::runtime_error(msg.c_str());
    }

    bconfig.set(BCO_PARAMS::NC) = (nc_mass_fixing) ? bconfig(BCO_PARAMS::NC)
                                                   : tov->rhoc;
    bconfig.set(BCO_PARAMS::HC) =
        (hc_mass_fixing) ? bconfig(BCO_PARAMS::HC)
                         : EOS<eos_t, eos_var_t::PRESSURE>::h_cold__rho(
                               bconfig(BCO_PARAMS::NC));
    bconfig.set(BCO_PARAMS::MB) = tov->baryon_mass;
    bconfig.set(BCO_PARAMS::MADM) =
        (adm_mass_fixing) ? bconfig(BCO_PARAMS::MADM) : tov->mass;

    // update surface radius estimate
    bconfig.set(BCO_PARAMS::RMID) = tov->radius;
    bconfig.set(BCO_PARAMS::RIN) = 0.5 * bconfig(BCO_PARAMS::RMID);
    bconfig.set(BCO_PARAMS::ROUT) = bco_utils::gold_ratio *
                                    bconfig(BCO_PARAMS::RMID);

    return tov;
}
}  // namespace Kadath::FUKA_Solvers