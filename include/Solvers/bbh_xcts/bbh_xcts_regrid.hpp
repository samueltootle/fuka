/*
 * Copyright 2022
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
#include "bco_utilities.hpp"
#include "kadath.hpp"

/**
 * \addtogroup BBH_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
namespace bco_u = ::Kadath::bco_utils;
using namespace Kadath::FUKA_Config;
using config_t = kadath_config_boost<BIN_INFO>;

inline void update_bin_config(config_t& bconfig,
                              const Space_bin_bh& space,
                              const Scalar& conf,
                              const Scalar& lapse) {
  if (std::isnan(bconfig.set(OUTER_SHELLS)))
    bconfig.set(OUTER_SHELLS) = 0;

  bconfig.set(Q) = bconfig(MCH, BCO1) / bconfig(MCH, BCO2);
  bconfig.set(BIN_PARAMS::Q) = (bconfig.set(BIN_PARAMS::Q) > 1.0)
                                   ? 1.0 / bconfig.set(BIN_PARAMS::Q)
                                   : bconfig.set(BIN_PARAMS::Q);

  if (!bconfig.control(USE_CONFIG_VARS)) {

    bconfig.set(REXT) = 2. * bconfig(DIST);

    const double r_max_tot =
        (bconfig(BCO_PARAMS::RMID, NODES::BCO1) > bconfig(RMID, NODES::BCO2))
            ? bconfig(BCO_PARAMS::RMID, NODES::BCO1)
            : bconfig(BCO_PARAMS::RMID, NODES::BCO2);

    const double rout_sep_est =
        (bconfig(BIN_PARAMS::DIST) / 2. - r_max_tot) / 3. + r_max_tot;
    const double rout_min_est = bco_u::gold_ratio * r_max_tot;

    bconfig.set(BCO_PARAMS::ROUT, NODES::BCO1) = rout_sep_est;
    //  (rout_sep_est < rout_min_est) ? rout_min_est : rout_sep_est;
    bconfig.set(BCO_PARAMS::ROUT, NODES::BCO2) =
        bconfig(BCO_PARAMS::ROUT, NODES::BCO1);

    const double q = bconfig.set(BIN_PARAMS::Q);
    bconfig.set(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO1) =
        1.7 * std::pow(2.0, bco_u::shell_factor / q);
    bconfig.set(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO2) =
        bconfig(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO1);
  }
}

inline int bbh_xcts_regrid(config_t& bconfig, std::string outputfile) {
  int exit_status = EXIT_SUCCESS;
  std::string kadath_filename = bconfig.space_filename();

  FILE* fin = fopen(kadath_filename.c_str(), "r");
  Space_bin_bh old_space(fin);
  Scalar old_conf(old_space, fin);
  Scalar old_lapse(old_space, fin);
  Vector old_shift(old_space, fin);
  fclose(fin);

  // Update config vars
  // This control was mainly for testing
  if (!bconfig.control(USE_CONFIG_VARS)) {
    bco_u::update_config_BH_radii(old_space, bconfig, old_space.BH1 + 1,
                                  old_conf, NODES::BCO1);
    bco_u::update_config_BH_radii(old_space, bconfig, old_space.BH2 + 1,
                                  old_conf, NODES::BCO2);
  }  // end updating config vars

  update_bin_config(bconfig, old_space, old_conf, old_lapse);

  int ndom = old_space.get_nbr_domains();

  int res = bconfig.set(BIN_RES);

  if ((res % 2) == 0) {
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
              << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  int type_coloc = old_space.get_type_base();

  // setup bounds
  std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) * (1. + e * 0.25);

  std::vector<int> exclusion_doms{old_space.BH1, old_space.BH1 + 1,
                                  old_space.BH2, old_space.BH2 + 1};
  auto drPsi(compute_drPsi(old_space, old_conf,
                           Metric_flat(old_space, old_shift.get_basis()),
                           exclusion_doms, old_space.OUTER));
  std::vector<double> BH1_bounds;
  std::vector<double> BH2_bounds;
  BH1_bounds =
      bco_u::set_arb_boundsv3(bconfig, drPsi, old_space.BH1 + 2, NODES::BCO1);
  BH2_bounds =
      bco_u::set_arb_boundsv3(bconfig, drPsi, old_space.BH2 + 2, NODES::BCO2);
  // end setup bounds

  Space_bin_bh space(type_coloc, bconfig(DIST), BH1_bounds, BH2_bounds,
                     out_bounds, bconfig(BIN_RES));
  Base_tensor basis(space, CARTESIAN_BASIS);

  std::cout << "Resolution of old space: ";
  bco_u::print_constant_space_resolution(old_space);

  std::cout << "Resolution of new space: ";
  bco_u::print_constant_space_resolution(space);

  std::cout << "\nold bounds:" << std::endl;
  bco_u::print_bounds_from_space(old_space);

  std::cout << "New bounds:" << std::endl;
  bco_u::print_bounds_from_space(space);

  // needed in some cases, since there is no data in 0,1 and the interpolation
  // can go crazy
  std::array<const int, 2> old_nuc_doms{old_space.BH1, old_space.BH2};
  std::array<const Domain_shell_outer_homothetic*, 2> old_outer_homothetic{
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          old_space.get_domain(old_space.BH1 + 1)),
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          old_space.get_domain(old_space.BH2 + 1))};
  for (auto& i : {0, 1}) {
    // update BH fields to help with interpolation later
    bco_u::update_adapted_field(old_conf, old_nuc_doms[i] + 2,
                                old_nuc_doms[i] + 1, old_outer_homothetic[i],
                                OUTER_BC);
    bco_u::update_adapted_field(old_lapse, old_nuc_doms[i] + 2,
                                old_nuc_doms[i] + 1, old_outer_homothetic[i],
                                OUTER_BC);
    for (int j = 1; j < 4; ++j)
      bco_u::update_adapted_field(old_shift.set(j), old_nuc_doms[i] + 2,
                                  old_nuc_doms[i] + 1, old_outer_homothetic[i],
                                  OUTER_BC);
  }

  // setup new fields
  Scalar conf(space);
  conf = 1.;

  Scalar lapse(space);
  lapse = 1.;

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();
  // end setup

  // import old fields into new space
  conf.import(old_conf);
  lapse.import(old_lapse);

  for (int c = 1; c <= 3; ++c)
    shift.set(c).import(old_shift(c));
  // end import fields

  // make sure excised domains are set to zero
  for (auto& dom : {space.BH1, space.BH2}) {
    lapse.set_domain(dom).annule_hard();
    lapse.set_domain(dom + 1).annule_hard();

    conf.set_domain(dom).annule_hard();
    conf.set_domain(dom + 1).annule_hard();

    for (int i = 1; i <= 3; ++i) {
      shift.set(i).set_domain(dom).annule_hard();
      shift.set(i).set_domain(dom + 1).annule_hard();
    }
  }

  lapse.std_base();
  conf.std_base();
  shift.std_base();

  bconfig.set_filename(outputfile);
  bco_u::save_to_file(space, bconfig, conf, lapse, shift);
  return exit_status;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath