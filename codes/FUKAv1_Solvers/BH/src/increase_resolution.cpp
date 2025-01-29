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
#include <memory>
#include <sstream>
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "kadath_adapted_bh.hpp"
#include "mpi.h"

using namespace Kadath;
using namespace Kadath::FUKA_Config;

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 2) {
    std::cerr << "Usage: ./increase_resol /<path>/<ID base name>.info "
                 "./<new_path>/<new ID base name> <Resolution>"
              << endl;
    std::cerr << "Ex: ./increase_resol converged.9.info initbh 11" << endl;
    std::_Exit(EXIT_FAILURE);
  } else if (argc < 3) {
    std::cerr << "Output file base name missing...(e.g. initbh)" << endl;
    std::_Exit(EXIT_FAILURE);
  } else if (argc < 4) {
    std::cerr << "Missing new resolution...(e.g. 11)" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }
  std::string outputfile{argv[2]};

  std::string in_filename = argv[1];
  kadath_config_boost<BCO_BH_INFO> bconfig(in_filename);

  std::string in_spacefile = bconfig.space_filename();
  FILE* ff1 = fopen(in_spacefile.c_str(), "r");

  Space_adapted_bh old_space(ff1);
  Scalar old_conf(old_space, ff1);
  Scalar old_lapse(old_space, ff1);
  Vector old_shift(old_space, ff1);

  fclose(ff1);

  std::cout << "Resolution of old space: "
            << old_space.get_domain(0)->get_nbr_points()(0) << " (r), "
            << old_space.get_domain(0)->get_nbr_points()(1) << " (theta), "
            << old_space.get_domain(0)->get_nbr_points()(2) << " (phi)"
            << std::endl;

  int ndim = 3;
  int old_ndom = old_space.get_nbr_domains();

  bconfig.set(BCO_RES) = std::atoi(argv[3]);
  /**
   * Make sure the domains have adequate space between them based on
   * the previous solution
   */
  if (!bconfig.control(USE_CONFIG_VARS)) {
    bconfig.set(RMID) =
        bco_utils::get_radius(old_space.get_domain(1), OUTER_BC);

    //estimate how small the inner radius should be based on relation
    //between conformal factor and numerical radius.
    //see https://arxiv.org/pdf/0805.4192, eq(64)
    double conf_inner = bco_utils::get_boundary_val(2, old_conf, INNER_BC);
    double conf_i_sq = conf_inner * conf_inner;
    double est_r_div2 = bconfig(MCH) / conf_i_sq;
    bconfig.set(RIN) = est_r_div2;
    bconfig.set(ROUT) = 4. * est_r_div2 * 2.;
  }

  int ndom = 4 + bconfig(NSHELLS);
  std::cout << "Number of Domains: " << ndom << std::endl;
  std::vector<double> bounds(ndom - 1);
  bco_utils::set_isolated_BH_bounds(bounds, bconfig);

  // FIXME not sure if it's only about oddness...
  if (((int)bconfig(BCO_RES) % 2) == 0) {
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
              << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::cout << "Resolution of new space: " << bconfig(BCO_RES) << " (r), "
            << bconfig(BCO_RES) << " (theta), " << bconfig(BCO_RES) - 1
            << " (phi)" << std::endl;
  Dim_array res(ndim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = bconfig(BCO_RES);
  res.set(2) = bconfig(BCO_RES) - 1;

  int type_coloc = old_space.get_type_base();

  Point center(old_space.get_domain(0)->get_center());

  Space_adapted_bh space(type_coloc, center, res, bounds);

  Base_tensor basis(space, CARTESIAN_BASIS);

  std::cout << "old bounds:" << std::endl;
  for (int i = 0; i < old_ndom; ++i)
    std::cout << bco_utils::get_radius(old_space.get_domain(i), OUTER_BC)
              << " ";
  std::cout << std::endl;
  std::cout << "New bounds:" << std::endl;
  for (int i = 0; i < ndom; ++i)
    std::cout << bco_utils::get_radius(space.get_domain(i), OUTER_BC) << " ";
  std::cout << std::endl;

  //setup new fields
  Scalar conf(space);
  conf = 1.;
  conf.std_base();

  Scalar lapse(space);
  lapse = 1.;
  lapse.std_base();

  Vector shift(space, CON, basis);
  shift.annule_hard();
  shift.std_base();
  //end setup new fields

  // needed in some cases, since there is no data in domains [0,1] and the interpolation can go crazy
  const Domain_shell_outer_homothetic* old_outer_homothetic =
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          old_space.get_domain(1));

  //import fields
  bco_utils::update_adapted_field(old_conf, 2, 1, old_outer_homothetic,
                                  OUTER_BC);
  bco_utils::update_adapted_field(old_lapse, 2, 1, old_outer_homothetic,
                                  OUTER_BC);
  for (int i = 1; i <= 3; ++i)
    bco_utils::update_adapted_field(old_shift.set(i), 2, 1,
                                    old_outer_homothetic, OUTER_BC);

  conf.import(old_conf);
  lapse.import(old_lapse);

  shift.set(1).import(old_shift.set(1));
  shift.set(2).import(old_shift.set(2));
  shift.set(3).import(old_shift.set(3));
  //end import fields

  //reset fields to zero inside excised region
  for (auto& d : {0, 1}) {
    lapse.set_domain(d).annule_hard();
    conf.set_domain(d).annule_hard();
    for (int i = 1; i <= 3; ++i) {
      shift.set(i).set_domain(d).annule_hard();
    }
  }
  lapse.std_base();
  conf.std_base();
  shift.std_base();

  bconfig.set_filename(outputfile);
  if (rank == 0)
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
