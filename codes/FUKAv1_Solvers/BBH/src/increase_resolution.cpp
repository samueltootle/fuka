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
#include <math.h>
#include <cmath>
#include <sstream>
#include <vector>
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath_bin_bh.hpp"
using namespace Kadath;
using namespace Kadath::FUKA_Config;

template <typename config_t, typename space_t>
void setup_bin_config(config_t& bconfig,
                      const space_t& space,
                      const Scalar& conf,
                      const Scalar& lapse) {
  int res = bconfig(BIN_RES);

  // FIXME not sure if it's only about oddness...
  if ((res % 2) == 0) {
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
              << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  if (std::isnan(bconfig.set(OUTER_SHELLS)))
    bconfig.set(OUTER_SHELLS) = 0;

  bconfig.set(BCO_RES, BCO1) = res;
  bconfig.set(BCO_RES, BCO2) = res;

  std::cout << "Resolution of new space: " << res << " (r), " << res
            << " (theta), " << res - 1 << " (phi)" << std::endl;

  double Mtot = bconfig(MCH, BCO1) + bconfig(MCH, BCO2);

  bconfig.set(MCH, BCO1) = Mtot / (bconfig(Q) + 1);
  bconfig.set(MIRR, BCO1) =
      bco_utils::mirr_from_mch(bconfig(CHI, BCO1), bconfig(MCH, BCO1));
  bconfig.set(MCH, BCO2) = bconfig(Q) * bconfig(MCH, BCO1);
  bconfig.set(MIRR, BCO2) =
      bco_utils::mirr_from_mch(bconfig(CHI, BCO2), bconfig(MCH, BCO2));
  double kep_omega = std::sqrt((bconfig(MCH, BCO1) + bconfig(MCH, BCO2)) / 2. /
                               (bconfig(DIST) * bconfig(DIST) * bconfig(DIST)));

  std::cout << "\nTotal Mass of the System: " << Mtot << std::endl;
  std::cout << "Keplarian Omega (for reference): " << kep_omega << std::endl;

  if (!bconfig.control(USE_CONFIG_VARS)) {
    double rmax = 0;
    int i = BCO1;

    for (auto& d : {space.BH1 + 1, space.BH2 + 1}) {
      auto [rmin, trmax] = bco_utils::get_rmin_rmax(space, d);

      // estimate how small the inner radius should be based on relation
      // between conformal factor and numerical radius.
      // see https://arxiv.org/pdf/0805.4192, eq(64)
      double conf_inner = bco_utils::get_boundary_val(d + 1, conf, INNER_BC);
      double conf_i_sq = conf_inner * conf_inner;
      double est_r_div2 = bconfig(MCH, i) / conf_i_sq;
      bconfig.set(RIN, i) = est_r_div2;

      // set this here only for shell bounds to be calculated.
      bconfig.set(RMID, i) = 2. * est_r_div2;

      auto [fmin, fmax] = bco_utils::get_field_min_max(lapse, 2, INNER_BC);
      bconfig.set(FIXED_LAPSE, i) = fmin;

      i = BCO2;
      if (trmax > rmax)
        rmax = trmax;
    }

    bconfig.set(REXT) = 2. * bconfig(DIST);
  }
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout
        << "Missing arguments..." << std::endl
        << "Possible input options: ./increase_resol <filename.info> <base "
           "output filename> <resolution>\n"
        << "NOTES:\n"
        << "all remaining attributes are pulled from the input config file.\n"
        << "- shells per BH are pulled from bh<1/2>.nshells\n"
        << "- mass ratio is pull from binary.q\n"
        << "- total mass is determined from bh1.mch and bh2.mch\n"
        << "- enabled 'use_config_vars' for pulling numerical domain values "
           "from the config file\n";

    std::_Exit(EXIT_FAILURE);
  }

  // Name of config.info file
  std::string in_filename{argv[1]};
  std::string of_basename{argv[2]};
  kadath_config_boost<BIN_INFO> bconfig(in_filename);
  bconfig.set(BIN_RES) = stod(argv[3]);

  std::string in_spacefile = bconfig.space_filename();
  FILE* ff1 = fopen(in_spacefile.c_str(), "r");

  Space_bin_bh old_space(ff1);
  Scalar old_conf(old_space, ff1);
  Scalar old_lapse(old_space, ff1);
  Vector old_shift(old_space, ff1);
  fclose(ff1);
  setup_bin_config(bconfig, old_space, old_conf, old_lapse);

  std::cout << "Resolution of old space: "
            << old_space.get_domain(0)->get_nbr_points()(0) << " (r), "
            << old_space.get_domain(0)->get_nbr_points()(1) << " (theta), "
            << old_space.get_domain(0)->get_nbr_points()(2) << " (phi)"
            << std::endl;

  double Mtot = bconfig(MCH, BCO1) + bconfig(MCH, BCO2);

  std::cout << "\nTotal Mass of the System: " << Mtot << std::endl;

  int type_coloc = old_space.get_type_base();

  //setup bounds
  std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) * (1. + e * 0.25);

  std::vector<double> BH1_bounds(3 + bconfig(NSHELLS, BCO1));
  std::vector<double> BH2_bounds(3 + bconfig(NSHELLS, BCO2));
  bco_utils::set_BH_bounds(BH1_bounds, bconfig, BCO1);
  bco_utils::set_BH_bounds(BH2_bounds, bconfig, BCO2);

  // Set radius of the excision boundary to the current radius so that the solver
  // starts from the originial solution
  BH1_bounds[1] =
      bco_utils::get_radius(old_space.get_domain(old_space.BH1 + 1), OUTER_BC);
  BH2_bounds[1] =
      bco_utils::get_radius(old_space.get_domain(old_space.BH2 + 1), OUTER_BC);
  //end setup bounds

  Space_bin_bh space(type_coloc, bconfig(DIST), BH1_bounds, BH2_bounds,
                     out_bounds, bconfig(BIN_RES));
  Base_tensor basis(space, CARTESIAN_BASIS);

  std::cout << "New bounds:" << std::endl;
  for (int d = space.BH1; d < space.BH2; ++d)
    std::cout << bco_utils::get_radius(space.get_domain(d), OUTER_BC) << " ";
  std::cout << std::endl;
  for (int d = space.BH2; d < space.OUTER; ++d)
    std::cout << bco_utils::get_radius(space.get_domain(d), OUTER_BC) << " ";
  std::cout << std::endl;

  double xc1 = bco_utils::get_center(space, space.BH1);
  double xc2 = bco_utils::get_center(space, space.BH2);

  std::cout << "xc1: " << xc1 << std::endl;
  std::cout << "xc2: " << xc2 << std::endl;

  //setup new fields
  Scalar conf(space);
  conf = 1.;
  conf.std_base();

  Scalar lapse(space);
  lapse = 1.;
  lapse.std_base();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();
  shift.std_base();
  // end setup new fields

  // needed in some cases, since there is no data in 0,1 and the interpolation can go crazy
  std::array<const int, 2> old_nuc_doms{old_space.BH1, old_space.BH2};
  std::array<const Domain_shell_outer_homothetic*, 2> old_outer_homothetic{
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          old_space.get_domain(old_space.BH1 + 1)),
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          old_space.get_domain(old_space.BH2 + 1))};
  for (auto& i : {0, 1}) {
    // update BH fields to help with interpolation later
    bco_utils::update_adapted_field(old_conf, old_nuc_doms[i] + 2,
                                    old_nuc_doms[i] + 1,
                                    old_outer_homothetic[i], OUTER_BC);
    bco_utils::update_adapted_field(old_lapse, old_nuc_doms[i] + 2,
                                    old_nuc_doms[i] + 1,
                                    old_outer_homothetic[i], OUTER_BC);
    for (int j = 1; j < 4; ++j)
      bco_utils::update_adapted_field(old_shift.set(j), old_nuc_doms[i] + 2,
                                      old_nuc_doms[i] + 1,
                                      old_outer_homothetic[i], OUTER_BC);
  }

  // import fields
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
  // end set excised domains

  lapse.std_base();
  conf.std_base();
  shift.std_base();

  bconfig.set_filename(of_basename);
  bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
}
