/*
 * Copyright 2023
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
#include <memory>
#include "Solvers/sequences/ns_sequence.hpp"
#include "Solvers/solvers.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template <class eos_t,
          typename config_t,
          typename space_t = Space_polar_adapted>
class ns_isotropic_diff_rot_solver : public Solver<config_t, space_t> {
   public:
    using typename Solver<config_t, space_t>::base_config_t;
    using typename Solver<config_t, space_t>::base_space_t;

   private:
    Scalar& nu;
    Scalar& lap_Aterm;
    Scalar& logh;
    Scalar& lap_Bterm;
    Scalar& lap_wterm;
    Scalar& Omega;
    std::unique_ptr<ns_sequence const> seq;
    std::string law{};

    /// Specify base class members used to avoid this->
    using Solver<config_t, space_t>::space;
    using Solver<config_t, space_t>::bconfig;
    using Solver<config_t, space_t>::ndom;
    using Solver<config_t, space_t>::check_max_iter_exceeded;
    using Solver<config_t, space_t>::solution_exists;
    using Solver<config_t, space_t>::extract_eos_name;
    using Solver<config_t, space_t>::checkpoint;
    using Solver<config_t, space_t>::solver_stage;

   public:
    /// solver is not trivially constructable since Kadath containers are not
    /// trivially constructable
    ns_isotropic_diff_rot_solver() = delete;

    ns_isotropic_diff_rot_solver(config_t& config_in,
                                 space_t& space_in,
                                 Scalar& nu_in,
                                 Scalar& lap_Aterm_in,
                                 Scalar& logh_in,
                                 Scalar& lap_Bterm_in,
                                 Scalar& lap_wterm_in,
                                 Scalar& Omega_in);

    /// syst always requires the same initialization for the stages
    void syst_init(System_of_eqs& syst);

    /// diagnostics at runtime
    void print_diagnostics(const System_of_eqs& syst,
                           const int ite = 0,
                           const double conv = 0) const override;

    std::string converged_filename(const std::string stage = "") const override;

    void save_to_file() const override {
        Kadath::bco_utils::save_to_file(space,
                                        bconfig,
                                        lap_Aterm,
                                        nu,
                                        logh,
                                        lap_Bterm,
                                        lap_wterm,
                                        Omega);
    }

    /// solver driver
    int solve();
    int solve(ns_sequence const* sequence_in);

    /// solver stage - depending on rotation law
    int keh_stage();

    void update_config_quantities(System_of_eqs& syst) {
        // if(law != "keh") {
        //   auto rs = bco_utils::get_rmin_rmax(space, 1);
        //   bconfig.set(BCO_PARAMS::RMID) = rs[0];
        // }

        // compute the ADM mass as surface integral at infinity
        Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
        double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
        bconfig.set(BCO_PARAMS::MADM) = Madm;

        // compute the baryonic mass at volume integral from the given integrant
        double baryonic_mass = syst.give_val_def("intMb")()(0).integ_volume() +
                               syst.give_val_def("intMb")()(1).integ_volume();
        bconfig.set(BCO_PARAMS::MB) = baryonic_mass;

        auto loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);

        // compute the ADM Angular Momentum as surface integral at infinity
        Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
        double Jadm = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);
        double chi = Jadm / Madm / Madm;

        if (seq) {
            auto idx{seq->mass_idx()};
            switch (idx) {
                case BCO_PARAMS::HC:
                    bconfig.set(BCO_PARAMS::NC) =
                        EOS<eos_t, DENSITY>::get(bconfig(BCO_PARAMS::HC));
                    break;
                case BCO_PARAMS::NC:
                    bconfig.set(BCO_PARAMS::HC) = std::exp(loghc);
                    break;
                default:
                    bconfig.set(BCO_PARAMS::HC) = std::exp(loghc);
                    bconfig.set(BCO_PARAMS::NC) =
                        EOS<eos_t, DENSITY>::get(bconfig(BCO_PARAMS::HC));
                    break;
            }

            idx = seq->spin_idx();
            switch (idx) {
                // case BCO_PARAMS::CHI:
                //   break;
                default:
                    bconfig.set(BCO_PARAMS::CHI) = chi;
                    break;
            }
        } else {
            bconfig.set(BCO_PARAMS::NC) =
                EOS<eos_t, DENSITY>::get(bconfig(BCO_PARAMS::HC));
            bconfig.set(BCO_PARAMS::CHI) = chi;
        }
    }
};

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath

#include "ns_isotropic_diff_rot_solver_imp.cpp"
#include "ns_isotropic_diff_rot_stages.cpp"