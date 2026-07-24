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
#pragma once
#include "Solvers/solvers.hpp"
#include "kadath_bin_ns.hpp"

/**
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template <class eos_t,
          typename config_t,
          typename space_t = Kadath::Space_bin_ns>
class bns_xcts_solver : XCTS_Solver<config_t, space_t> {
   public:
    using typename XCTS_Solver<config_t, space_t>::base_config_t;
    using typename XCTS_Solver<config_t, space_t>::base_space_t;

   private:
    Scalar& conf;
    Scalar& lapse;
    Scalar& logh;
    Scalar& phi;
    Vector& shift;
    Metric_flat fmet;

    const double xc1;
    const double xc2;
    const double xo{0.};

    // Specify base class members used to avoid this->
    using XCTS_Solver<config_t, space_t>::space;
    using XCTS_Solver<config_t, space_t>::bconfig;
    using XCTS_Solver<config_t, space_t>::basis;
    using XCTS_Solver<config_t, space_t>::cfields;
    using XCTS_Solver<config_t, space_t>::coord_vectors;
    using XCTS_Solver<config_t, space_t>::ndom;
    using XCTS_Solver<config_t, space_t>::check_max_iter_exceeded;
    using XCTS_Solver<config_t, space_t>::solution_exists;
    using XCTS_Solver<config_t, space_t>::extract_eos_name;
    using XCTS_Solver<config_t, space_t>::checkpoint;
    using XCTS_Solver<config_t, space_t>::solver_stage;

   public:
    // solver is not trivially constructable since Kadath containers are not
    // trivially constructable
    bns_xcts_solver() = delete;

    bns_xcts_solver(config_t& config_in,
                    space_t& space_in,
                    Base_tensor& base_in,
                    Scalar& conf_in,
                    Scalar& lapse_in,
                    Vector& shift_in,
                    Scalar& logh_in,
                    Scalar& phi_in);

    // syst always requires the same initialization for the stages
    void syst_init(System_of_eqs& syst);
    void syst_init_Aterms(System_of_eqs& syst);

    // diagnostics at runtime
    void print_diagnostics(const System_of_eqs& syst,
                           const int ite = 0,
                           const double conv = 0) const override;

    std::string converged_filename(const std::string stage = "") const override;

    void save_to_file() const override {
        ::Kadath::bco_utils::save_to_file(space,
                                          bconfig,
                                          conf,
                                          lapse,
                                          shift,
                                          logh,
                                          phi);
    }

    // solve driver
    int solve();

    /**
   * hydrostatic_equilibrium stage
   *
   * The binary is solved using the force balance equations
   * in addition to solving the relativisitc Euler equation
   * to obtain a binary in hydrostatic equilibrium.
   */
    int hydrostatic_equilibrium_stage();
    /**
   * hydro_rescaling_stages
   *
   * The binary is solved using fixed orbital velocity
   * therefore the matter scalar fields are simply rescaled
   * based on the fixed baryonic mass of each star.
   *
   * @param[input] stage_text: Some changes are made based on TOTAL_BC or
   * ECC_RED stage.
   */
    int hydro_rescaling_stages(std::string stage_text);

    // Update bconfig(HC) and bconfig(NC)
    void update_config_quantities(System_of_eqs& syst);
};

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath

#include "bns_xcts_solver_imp.cpp"
#include "bns_xcts_stages.cpp"