/*
 * Copyright 2021
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
#include <sstream>
#include "coord_fields.hpp"
#include "kadath_bin_bh.hpp"

#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"

#include <iostream>
#include <string>
#include <utility>
using namespace Kadath;
using namespace Kadath::FUKA_Config;

auto set_bin_defaults() {
  /******* Equal Mass, irrot Initial Data Setup *********************/

  kadath_config_boost<BCO_BH_INFO> bhconfig;
  bhconfig.set_defaults();

  kadath_config_boost<BIN_INFO> bconfig;
  bconfig.initialize_binary({"bh", "bh"});

  // copy default BH configuration to both BHs in binary
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO1) = bhconfig.set(i);
    bconfig.set(i, BCO2) = bhconfig.set(i);
  }

  // setup binary parameters
  bconfig.set(BIN_RES) = 9;
  bconfig.set(DIST) = 10;
  bconfig.set(REXT) = 2 * bconfig(DIST);
  bconfig.set(QPIG) = 4 * M_PI;
  bconfig.set(Q) = 1;
  bconfig.set(ADOT) = 0;
  bconfig.set(OUTER_SHELLS) = 0;
  return bconfig;
}

int main(int argc, char** argv) {
  kadath_config_boost<BIN_INFO> bconfig;

  if (argc < 2) {
    std::cout << "Using default setup.\n";
    bconfig = set_bin_defaults();
  } else {
    std::string input_filename = std::string{argv[1]};
    std::cout << "Creating setup from: " << input_filename << "\n";
    bconfig.set_filename(input_filename);
    bconfig.open_config();
  }

  std::stringstream ifilename;
  ifilename << "./initbin.info";
  bconfig.set_filename(ifilename.str());

  std::vector<double> out_bounds(1);
  std::vector<double> BH1_bounds(3 + bconfig(NSHELLS, BCO2));
  std::vector<double> BH2_bounds(3 + bconfig(NSHELLS, BCO2));

  // for out_bounds.size > 1 - add equi-distant shells
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) + e * 0.25 * bconfig(REXT);

  bco_utils::set_BH_bounds(BH1_bounds, bconfig, BCO1);
  bco_utils::set_BH_bounds(BH2_bounds, bconfig, BCO2);

  // create space containing the domain decomposition
  int type_coloc = CHEB_TYPE;
  Space_bin_bh space(type_coloc, bconfig(DIST), BH1_bounds, BH2_bounds,
                     out_bounds, bconfig(BIN_RES));
  Base_tensor basis(space, CARTESIAN_BASIS);

  double xmin = bco_utils::get_center(space, space.BH1);
  double xmax = bco_utils::get_center(space, space.BH2);

  std::cout << "xmin/max: " << xmin << " " << xmax << std::endl;

  // remaining vars based on basic Newtonian estimates
  if (!bconfig.control(USE_CONFIG_VARS)) {
    bconfig.set(COMY) = 0;
    bconfig.set(COM) = 0;
    bco_utils::KadathPNOrbitalParams(bconfig, bconfig(MCH, BCO1),
                                     bconfig(MCH, BCO2));
    bconfig.set(OMEGA, BCO1) = bconfig(GOMEGA);  //irrotational
    bconfig.set(OMEGA, BCO2) = bconfig(GOMEGA);

    bconfig.set_field(CONF) = true;
    bconfig.set_field(LAPSE) = true;
    bconfig.set_field(SHIFT) = true;

    bconfig.set_stage(PRE) = true;
    bconfig.set_stage(FIXED_OMEGA) = true;
    bconfig.set_stage(COROT_EQUAL) = true;
    bconfig.set_stage(TOTAL) = true;
    bconfig.set_stage(TOTAL_BC) = true;
  }
  // END OF CONFIG FILE SETUP

  // setup fields to flat space values
  Scalar lapse(space);
  lapse = 1.;
  for (auto& dom : {space.BH1, space.BH2}) {
    lapse.set_domain(dom).annule_hard();
    lapse.set_domain(dom + 1).annule_hard();
  }
  Scalar conf(lapse);

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();

  conf.std_base();
  lapse.std_base();
  shift.std_base();
  // end setup of fields

  bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  return EXIT_SUCCESS;
}
