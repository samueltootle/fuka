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
#include "kadath_adapted_bh.hpp"

#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>
using namespace Kadath;
using namespace Kadath::FUKA_Config;

int main(int argc, char** argv) {
  kadath_config_boost<BCO_BH_INFO> bconfig;
  if (argc < 2) {
    std::cout << "Using default setup.\n";
    bconfig.set_defaults();
  } else {
    std::string input_filename = std::string{argv[1]};
    std::cout << "Creating setup from: " << input_filename << "\n";
    bconfig.set_filename(input_filename);
    bconfig.open_config();
  }
  bconfig.set(RIN) = bconfig(MCH) * bco_utils::invpsisq;
  bconfig.set(RMID) = 2 * bconfig(MCH) * bco_utils::invpsisq;
  bconfig.set(ROUT) = 4 * bconfig(RMID);

  bconfig.set(BVELX) = 0.;
  bconfig.set(BVELY) = 0.;

  // set new base filename for new setup
  std::string base_filename = "./initbh.info";
  bconfig.set_filename(base_filename);

  const int dim = 3;

  int type_coloc = CHEB_TYPE;
  Dim_array res(dim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = bconfig(BCO_RES);
  res.set(2) = bconfig(BCO_RES) - 1;

  Point center(dim);
  for (int i = 1; i <= dim; i++)
    center.set(i) = 0;

  int ndom = 4 + bconfig(NSHELLS);
  std::cout << "Number of Domains: " << ndom << std::endl;

  // setup bounds
  std::vector<double> bounds(ndom - 1);
  bco_utils::set_isolated_BH_bounds(bounds, bconfig);

  Space_adapted_bh space(type_coloc, center, res, bounds);
  Base_tensor basis(space, CARTESIAN_BASIS);

  // setup fields
  Scalar conf(space);
  conf = 1.;
  conf.set_domain(0).annule_hard();
  conf.set_domain(1).annule_hard();

  Scalar lapse(conf);

  // set a better estimate for PSI on the horizon
  conf.set_domain(2) = bco_utils::psi;
  conf.std_base();
  lapse.std_base();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();
  shift.std_base();
  // end setup fields

  std::cout << bconfig << std::endl;
  bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  return EXIT_SUCCESS;
}
