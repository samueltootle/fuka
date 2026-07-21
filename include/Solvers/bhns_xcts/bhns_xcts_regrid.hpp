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
#pragma once
#include <math.h>
#include <sstream>
#include "Configurator/config_binary.hpp"
#include "FUKA_Solvers/utilities/scalar_calculations.hpp"
#include "Solvers/fuka_syst/fuka_syst.hpp"
#include "bco_utilities.hpp"
#include "kadath.hpp"

/**
 * \addtogroup BHNS_XCTS
 * \ingroup FUKA
 * @{*/
using namespace Kadath::FUKA_Config;

namespace Kadath {
namespace FUKA_Solvers {
using config_t = kadath_config_boost<BIN_INFO>;

inline int bhns_xcts_regrid(config_t& bconfig, std::string output_fname) {
    using namespace ::Kadath::bco_utils;
    bconfig.set(Q) = bconfig(MADM, BCO1) / bconfig(MCH, BCO2);
    bconfig.set(BIN_PARAMS::Q) = (bconfig.set(BIN_PARAMS::Q) > 1.0)
                                     ? 1.0 / bconfig.set(BIN_PARAMS::Q)
                                     : bconfig.set(BIN_PARAMS::Q);

    const double q = bconfig.set(BIN_PARAMS::Q);
    bconfig.set(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO1) =
        1.7 * std::pow(2.0, shell_factor / q);
    bconfig.set(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO2) =
        1.7 * std::pow(2.0, shell_factor / q);

    if (std::isnan(bconfig.set(OUTER_SHELLS)))
        bconfig.set(OUTER_SHELLS) = 0;

    std::string kadath_filename = bconfig.space_filename();

    FILE* fin = fopen(kadath_filename.c_str(), "r");
    Space_bhns old_space(fin);
    Scalar old_conf(old_space, fin);
    Scalar old_lapse(old_space, fin);
    Vector old_shift(old_space, fin);
    Scalar old_logh(old_space, fin);
    Scalar old_phi(old_space, fin);

    fclose(fin);

    const std::array<int, 2> old_adapted_doms{old_space.ADAPTEDNS,
                                              old_space.ADAPTEDBH};
    const std::array<int, 2> old_nuc_doms{old_space.NS, old_space.BH};

    const Domain_shell_outer_adapted* old_outer_adaptedNS =
        dynamic_cast<const Domain_shell_outer_adapted*>(
            old_space.get_domain(old_space.ADAPTEDNS));

    const Domain_shell_inner_adapted* old_inner_adaptedNS =
        dynamic_cast<const Domain_shell_inner_adapted*>(
            old_space.get_domain(old_space.ADAPTEDNS + 1));

    const Domain_shell_outer_homothetic* old_bh_outer =
        dynamic_cast<const Domain_shell_outer_homothetic*>(
            old_space.get_domain(old_space.ADAPTEDBH));

    int ndim = 3;
    int ndom = old_space.get_nbr_domains();

    int res = bconfig(BIN_RES);

    if ((res % 2) == 0) {
        std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
                  << std::endl;
        std::_Exit(EXIT_FAILURE);
    }

    int type_coloc = old_space.get_type_base();

    // Update config vars
    // This control was mainly for testing
    if (!bconfig.control(USE_CONFIG_VARS)) {
        update_config_NS_radii(old_space,
                               bconfig,
                               old_space.ADAPTEDNS,
                               NODES::BCO1);

        update_config_BH_radii(old_space,
                               bconfig,
                               old_space.ADAPTEDBH,
                               old_conf,
                               NODES::BCO2);

        double rmax = std::max(bconfig(RMID, BCO1), bconfig(RMID, BCO2));
        const double rout_sep_est = (bconfig(DIST) / 2. - rmax) / 3. + rmax;
        const double rout_max_est = gold_ratio * rmax;
        bconfig.set(ROUT, BCO1) = rout_sep_est;
        // (rout_sep_est > rout_max_est) ? rout_max_est : rout_sep_est;
        bconfig.set(ROUT, BCO2) = bconfig(ROUT, BCO1);
    }  // end updating config vars

    // create old radius scalar field
    Scalar old_space_radius(old_space);
    old_space_radius.annule_hard();
    for (int d = 0; d < ndom; ++d) {
        old_space_radius.set_domain(d) = old_space.get_domain(d)->get_radius();
    }

    // Get outer domain radii
    old_space_radius.set_domain(old_space.ADAPTEDNS) =
        old_outer_adaptedNS->get_outer_radius();
    old_space_radius.set_domain(old_space.ADAPTEDBH) =
        old_bh_outer->get_outer_radius();

    old_space_radius.std_base();
    // end creating old radius scalar fields

    // setup bounds for creating new space
    std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
    for (int e = 0; e < out_bounds.size(); ++e)
        out_bounds[e] = bconfig(REXT) * (1. + e * 0.25);

    std::vector<int> ns_interior_doms{
        FUKA_Syst_tools::vector_of_domains(old_space.NS, old_space.ADAPTEDNS)};
    std::vector<int> exclusion_doms{old_space.BH, old_space.BH + 1};
    // concat domain lists together
    std::for_each(ns_interior_doms.rbegin(),
                  ns_interior_doms.rend(),
                  [&exclusion_doms](auto e) {
                      auto it = exclusion_doms.begin();
                      exclusion_doms.insert(it, e);
                  });
    auto drPsi(compute_drPsi(old_space,
                             old_conf,
                             Metric_flat(old_space, old_shift.get_basis()),
                             exclusion_doms,
                             old_space.OUTER));
    std::vector<double> NS_bounds{
        set_arb_boundsv3(bconfig, drPsi, old_space.ADAPTEDNS + 1, NODES::BCO1)};

    std::vector<double> BH_bounds{
        set_arb_boundsv3(bconfig, drPsi, old_space.ADAPTEDBH + 1, NODES::BCO2)};
    // end setup bounds

    // print bounds to stdout - debugging only
    std::cout << "Local bounds:" << std::endl;
    print_bounds("NS-bounds", NS_bounds);
    print_bounds("BH-bounds", BH_bounds);
    print_bounds("outer-bounds", out_bounds);
    std::cout << std::endl;

    Space_bhns space(type_coloc,
                     bconfig(DIST),
                     NS_bounds,
                     BH_bounds,
                     out_bounds,
                     bconfig(BIN_RES),
                     bconfig(NINSHELLS, BCO1));
    Base_tensor basis(space, CARTESIAN_BASIS);

    std::cout << "Resolution of old space: ";
    print_constant_space_resolution(old_space);

    std::cout << "Resolution of new space: ";
    print_constant_space_resolution(space);

    std::cout << "\nOld bounds:" << std::endl;
    print_bounds_from_space(old_space);

    std::cout << "New bounds:" << std::endl;
    print_bounds_from_space(space);

    const Domain_shell_inner_adapted* new_ns_inner =
        dynamic_cast<const Domain_shell_inner_adapted*>(
            space.get_domain(space.ADAPTEDNS + 1));
    const Domain_shell_outer_adapted* new_ns_outer =
        dynamic_cast<const Domain_shell_outer_adapted*>(
            space.get_domain(space.ADAPTEDNS));

    // update BH fields to help with import
    update_adapted_field(old_conf,
                         old_space.ADAPTEDBH + 1,
                         old_space.ADAPTEDBH,
                         old_bh_outer,
                         OUTER_BC);
    update_adapted_field(old_lapse,
                         old_space.ADAPTEDBH + 1,
                         old_space.ADAPTEDBH,
                         old_bh_outer,
                         OUTER_BC);
    for (int i = 1; i < 4; ++i)
        update_adapted_field(old_shift.set(i),
                             old_space.ADAPTEDBH + 1,
                             old_space.ADAPTEDBH,
                             old_bh_outer,
                             OUTER_BC);

    update_adapted_field(old_phi,
                         old_space.ADAPTEDNS,
                         old_space.ADAPTEDNS + 1,
                         old_inner_adaptedNS,
                         INNER_BC);

    // Updated mapping for NS adapted fields
    interp_adapted_mapping(new_ns_inner, old_space.ADAPTEDNS, old_space_radius);
    interp_adapted_mapping(new_ns_outer, old_space.ADAPTEDNS, old_space_radius);

    // setup new fields
    Scalar conf(space);
    conf = 1.;

    Scalar lapse(space);
    lapse = 1.;

    Vector shift(space, CON, basis);
    for (int i = 1; i <= 3; i++)
        shift.set(i).annule_hard();

    Scalar logh(space);
    logh.annule_hard();

    Scalar phi(space);
    phi.annule_hard();
    // end setup

    // import old fields into new space
    conf.import(old_conf);
    lapse.import(old_lapse);
    logh.import(old_logh);
    phi.import(old_phi);

    shift.set(1).import(old_shift.set(1));
    shift.set(2).import(old_shift.set(2));
    shift.set(3).import(old_shift.set(3));
    // end import

    // make sure there is no matter or vel.pot. outside of the star
    logh.set_domain(space.ADAPTEDNS + 1).annule_hard();
    phi.set_domain(space.ADAPTEDNS + 1).annule_hard();
    for (int d = space.BH; d < space.get_nbr_domains(); ++d) {
        logh.set_domain(d).annule_hard();
        phi.set_domain(d).annule_hard();

        // make sure all fields are 0 inside the excision region
        if ((d >= space.BH) && (d < space.BH + 2)) {
            conf.set_domain(d).annule_hard();
            lapse.set_domain(d).annule_hard();
            for (int i = 1; i <= 3; i++)
                shift.set(i).set_domain(d).annule_hard();
        }
    }

    lapse.std_base();
    conf.std_base();
    logh.std_base();
    shift.std_base();
    phi.std_base();

    bconfig.set_filename(output_fname);

    save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    return EXIT_SUCCESS;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
