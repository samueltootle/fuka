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
#include <sstream>
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"
#include "kadath_bin_ns.hpp"
#include "mpi.h"
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif
using namespace Kadath;
using namespace Kadath::FUKA_Config;
using bin_space_t = Space_bin_ns;

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
    std::cerr << "Ex: ./increase_resol converged.9.info initbin 11" << endl;
    std::_Exit(EXIT_FAILURE);
  } else if (argc < 3) {
    std::cerr << "Output file base name missing...(e.g. initbin)" << endl;
    std::_Exit(EXIT_FAILURE);
  } else if (argc < 4) {
    std::cerr << "Missing new resolution...(e.g. 11)" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::string input_filename = argv[1];
  std::string outputfile{argv[2]};

  kadath_config_boost<BIN_INFO> bconfig(input_filename);
  if (std::isnan(bconfig.set(OUTER_SHELLS)))
    bconfig.set(OUTER_SHELLS) = 0;

  std::string kadath_filename = bconfig.space_filename();
  if (!fs::exists(kadath_filename)) {
    if (rank == 0) {
      std::cerr << "File: " << kadath_filename << " not found.\n\n";
      std::_Exit(EXIT_FAILURE);
    }
  }

  FILE* fin = fopen(kadath_filename.c_str(), "r");
  bin_space_t old_space(fin);
  Scalar old_conf(old_space, fin);
  Scalar old_lapse(old_space, fin);
  Vector old_shift(old_space, fin);
  Scalar old_logh(old_space, fin);
  Scalar old_phi(old_space, fin);

  fclose(fin);

  const std::array<int, 2> old_adapted_doms{old_space.ADAPTED1,
                                            old_space.ADAPTED2};
  const std::array<int, 2> old_nuc_doms{old_space.NS1, old_space.NS2};

  std::array<const Domain_shell_outer_adapted*, 2> old_outer_adapted;
  std::array<const Domain_shell_inner_adapted*, 2> old_inner_adapted;

  for (int i = 0; i < 2; ++i) {
    int const d = old_adapted_doms[i];
    old_outer_adapted[i] = dynamic_cast<const Domain_shell_outer_adapted*>(
        old_space.get_domain(d));
    old_inner_adapted[i] = dynamic_cast<const Domain_shell_inner_adapted*>(
        old_space.get_domain(d + 1));
  }

  bconfig.set(BIN_RES) = std::atoi(argv[3]);
  int res = bconfig.set(BIN_RES);
  if ((res % 2) == 0) {
    if (rank == 0)
      std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
                << std::endl;
    abort();
  }
  int ndim = 3;
  int ndom = old_space.get_nbr_domains();
  if (rank == 0) {
    std::cout << "Resolution of old space: "
              << old_space.get_domain(0)->get_nbr_points()(0) << " (r), "
              << old_space.get_domain(0)->get_nbr_points()(1) << " (theta), "
              << old_space.get_domain(0)->get_nbr_points()(2) << " (phi)"
              << std::endl;
    std::cout << "Resolution of new space: " << res << " (r), " << res
              << " (theta), " << res - 1 << " (phi)" << std::endl;
  }

  int type_coloc = old_space.get_type_base();

  //start Update config vars
  std::array<double, 2> r_min;
  double r_max_tot = 0.;

  if (rank == 0)
    std::cout << "Rmin/max: " << std::endl;
  for (int i = 0; i < 2; ++i) {
    int const dom = old_adapted_doms[i];

    auto [rmin, rmax] = bco_utils::get_rmin_rmax(old_space, dom);
    if (rank == 0)
      std::cout << rmin << " " << rmax << std::endl;

    bconfig.set(RIN, i) = 0.5 * rmin;
    bconfig.set(RMID, i) = rmin;

    r_max_tot = std::max(rmax, r_max_tot);
  }
  bconfig.set(ROUT, BCO1) = (bconfig(DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  bconfig.set(ROUT, BCO2) = (bconfig(DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  //end updating config vars

  //create old radius scalar field
  Scalar old_space_radius(old_space);
  old_space_radius.annule_hard();

  for (int d = 0; d < ndom; ++d) {
    old_space_radius.set_domain(d) = old_space.get_domain(d)->get_radius();
  }

  //Get adapted domain radii
  for (int i = 0; i < 2; ++i) {
    int const dom = old_adapted_doms[i];
    old_space_radius.set_domain(dom) = old_outer_adapted[i]->get_outer_radius();
  }
  old_space_radius.std_base();
  //end create old radius scalar fields

  //Setup necessary bounds
  std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
  std::vector<double> NS1_bounds(3 + bconfig(NSHELLS, BCO1));
  std::vector<double> NS2_bounds(3 + bconfig(NSHELLS, BCO2));

  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) * (1. + e * 0.25);

  bco_utils::set_NS_bounds(NS1_bounds, bconfig, BCO1);
  bco_utils::set_NS_bounds(NS2_bounds, bconfig, BCO2);

  if (rank == 0) {
    std::cout << "Bounds:" << std::endl;
    bco_utils::print_bounds("NS1", NS1_bounds);
    bco_utils::print_bounds("NS2", NS2_bounds);
  }
  //end set bounds

  bin_space_t space(type_coloc, bconfig(DIST), NS1_bounds, NS2_bounds,
                    out_bounds, res);
  ndom = space.get_nbr_domains();

  Base_tensor basis(space, CARTESIAN_BASIS);

  const std::array<int, 2> new_adapted_doms{space.ADAPTED1, space.ADAPTED2};
  const std::array<int, 2> new_nuc_doms{space.NS1, space.NS2};
  const std::array<double, 2> xc{bco_utils::get_center(space, space.NS1),
                                 bco_utils::get_center(space, space.NS2)};

  std::array<const Domain_shell_outer_adapted*, 2> new_outer_adapted;
  std::array<const Domain_shell_inner_adapted*, 2> new_inner_adapted;

  for (int i = 0; i < 2; ++i) {
    auto& d = new_adapted_doms[i];
    new_outer_adapted[i] =
        dynamic_cast<const Domain_shell_outer_adapted*>(space.get_domain(d));
    new_inner_adapted[i] = dynamic_cast<const Domain_shell_inner_adapted*>(
        space.get_domain(d + 1));
  }

  for (int i = 0; i < 2; ++i) {
    int const dom = old_adapted_doms[i];

    //Updated mapping for NS
    bco_utils::interp_adapted_mapping(new_inner_adapted[i], dom,
                                      old_space_radius);
    bco_utils::interp_adapted_mapping(new_outer_adapted[i], dom,
                                      old_space_radius);

    //interpolate old_phi field outside of the star for import
    bco_utils::update_adapted_field(old_phi, dom, dom + 1, old_inner_adapted[i],
                                    INNER_BC);
  }

  if (rank == 0) {
    std::cout << "xc1: " << xc[0] << std::endl;
    std::cout << "xc2: " << xc[1] << std::endl;
  }

  //start Create new fields
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

  Scalar logh(space);
  logh.annule_hard();
  logh.std_base();

  Scalar phi(space);
  phi.annule_hard();
  phi.std_base();
  //end Create new fields

  //start Import fields
  conf.import(old_conf);
  lapse.import(old_lapse);
  logh.import(old_logh);
  phi.import(old_phi);

  shift.set(1).import(old_shift.set(1));
  shift.set(2).import(old_shift.set(2));
  shift.set(3).import(old_shift.set(3));
  //end Import fields

  //start - set matter and velocity potential to zero outside stars
  for (auto& dom : new_adapted_doms) {
    phi.set_domain(dom + 1).annule_hard();
    logh.set_domain(dom + 1).annule_hard();
  }

  for (int i = space.OUTER; i < ndom; ++i) {
    phi.set_domain(i).annule_hard();
    logh.set_domain(i).annule_hard();
  }
  //end - set phi,logh to zero outside

  lapse.std_base();
  conf.std_base();
  logh.std_base();
  shift.std_base();
  phi.std_base();

  //in case output filename excludes a path
  bconfig.set_outputdir("./");

  bconfig.set_filename(outputfile);

  if (rank == 0)
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
