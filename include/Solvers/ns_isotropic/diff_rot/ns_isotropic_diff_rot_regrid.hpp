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
#include <math.h>
#include <sstream>
#include "Configurator/config_bco.hpp"
#include "FUKA_Solvers/utilities/compact_object_initializers/setup_2dns_isotropic.hpp"
#include "bco_utilities.hpp"
#include "kadath.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template <typename config_t>
int ns_isotropic_diff_rot_regrid(config_t& bconfig, std::string outputfile) {
  int exit_status = EXIT_SUCCESS;
  using space_t = Space_polar_adapted;

  std::string kadath_filename = bconfig.space_filename();

  if (!fs::exists(kadath_filename)) {
    std::cerr << "File: " << kadath_filename << " not found.\n\n";
    std::_Exit(EXIT_FAILURE);
  }

  FILE* ff1 = fopen(kadath_filename.c_str(), "r");
  space_t old_space(ff1);
  Scalar old_lap_Aterm(old_space, ff1);
  Scalar old_nu(old_space, ff1);
  Scalar old_logh(old_space, ff1);
  Scalar old_lap_Bterm(old_space, ff1);
  Scalar old_lap_wterm(old_space, ff1);
  Scalar old_Omega(old_space, ff1);
  fclose(ff1);

  std::cout << "Resolution of old space: "
            << old_space.get_domain(0)->get_nbr_points()(0) << " (r), "
            << old_space.get_domain(0)->get_nbr_points()(1) << " (theta)\n";

  int ndim = 2;

  // get the adapted domain and cast it to its correct type to be able to call
  // its member functions
  const Domain_polar_shell_outer_adapted* old_outer_adapted =
      dynamic_cast<const Domain_polar_shell_outer_adapted*>(
          old_space.get_domain(old_space.ADAPTED_OUTER));

  // setup a scalar field representing the old radius
  Scalar old_space_radius(old_space);
  old_space_radius = 0.;

  // get the radius from each domain
  for (int i = 0; i < old_space.get_nbr_domains(); ++i)
    old_space_radius.set_domain(i) = old_space.get_domain(i)->get_radius();
  // get the adapted radius of the adapted domain
  old_space_radius.set_domain(old_space.ADAPTED_OUTER) =
      old_outer_adapted->get_outer_radius();

  // define a standard decomposition, compatible with the parity of this field
  old_space_radius.std_base();
  // end setup old radius field

  // get the minimal and maximal radius from the adapted domain
  auto [r_min, r_max] = Kadath::bco_utils::get_rmin_rmax(old_space, 1);

  std::cout << "Rmin/max: " << r_min << " " << r_max << std::endl;

  // set new resolutions in each spatial dimension
  Dim_array res(ndim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = res(0);

  // FIXME not sure if it's only about oddness...
  if (res(0) % 2 == 0) {
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
              << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::cout << "Resolution of new space: " << res(0) << " (r), " << res(1)
            << " (theta)" << std::endl;

  // get the type of the colocation points
  int type_coloc = old_space.get_type_base();

  // ignore domain radii scaling from config if needed
  if (!bconfig.control(USE_CONFIG_VARS)) {
    bconfig.set(RIN) = 0.5 * r_min;
    bconfig.set(ROUT) = bco_utils::gold_ratio * r_max;
    bconfig.set(RMID) = r_max;
  }
  // end update config

  // setup radius bounds of the domains

  int const shells = bconfig(BCO_PARAMS::NSHELLS);
  int const ndom = 4 + shells;
  Array<double> bounds(ndom - 1);
  bounds.set(0) = bconfig(RIN);
  bounds.set(1) = bconfig(RMID);
  bounds.set(2) = bconfig(ROUT);

  for (int shell = 1, b = 3; shell <= shells; ++shell, ++b) {
    bounds.set(b) = bounds(b - 1) * bco_utils::gold_ratio;
  }

  // get origin of nucleus domain
  Point center = old_space.get_domain(0)->get_center();

  // initialize space with new resolution and domain decomposition
  Space_polar_adapted space(type_coloc, center, res, bounds);

  // get adapted domains to update the radius
  const Domain_polar_shell_outer_adapted* new_outer_adapted =
      dynamic_cast<const Domain_polar_shell_outer_adapted*>(
          space.get_domain(space.ADAPTED_OUTER));
  const Domain_polar_shell_inner_adapted* new_inner_adapted =
      dynamic_cast<const Domain_polar_shell_inner_adapted*>(
          space.get_domain(space.ADAPTED_INNER));

  // update adapted domain mapping
  Kadath::bco_utils::interp_adapted_mapping(new_outer_adapted, 1,
                                            old_space_radius);
  Kadath::bco_utils::interp_adapted_mapping(new_inner_adapted, 1,
                                            old_space_radius);

  // setup new fields
  // initialize to one or zero first
  Scalar lap_Aterm(space);
  lap_Aterm.annule_hard();
  lap_Aterm.std_base();

  Scalar one(space);
  one = 1.;
  one.std_base();

  Scalar rsint = Scalar(one.mult_sin_theta().mult_r());

  Scalar lap_Bterm(rsint);

  Scalar lap_wterm(lap_Bterm);

  Scalar nu(space);
  nu.annule_hard();
  nu.std_base();

  Scalar logh(space);
  logh.annule_hard();
  logh.std_base();

  Scalar Omega(space);
  Omega.annule_hard();
  // end setup new fields

  // import data from fields in the old space
  lap_Aterm.import(old_lap_Aterm);
  nu.import(old_nu);
  logh.import(old_logh);
  lap_Bterm.import(old_lap_Bterm);
  lap_wterm.import(old_lap_wterm);
  Omega.import(old_Omega);
  // end import old fields

  // enforce spectral decomposition compatible with the parities
  lap_Aterm.std_base();
  nu.std_base();
  logh.std_base();

  for (int d = 0; d < ndom; ++d) {
    lap_Bterm.set_domain(d).set_base() = rsint(d).get_base();
    lap_wterm.set_domain(d).set_base() = rsint(d).get_base();
  }

  Omega.std_base();

  // Fix outer boundary value for Omega
  auto npts_compact = space.get_domain(ndom - 1)->get_nbr_points();
  Index pos(npts_compact);
  pos.set(0) = npts_compact(0) - 1;
  Index posrminus1(pos);
  posrminus1.set(0) = npts_compact(0) - 2;
  for (auto i = 0; i < npts_compact(1); ++i) {
    pos.set(1) = i;
    posrminus1.set(1) = i;
    Omega.set_domain(ndom - 1).set(pos) = Omega(ndom - 1)(posrminus1);
  }

  // output data
  bconfig.set_filename(outputfile);
  bco_utils::save_to_file(space, bconfig, lap_Aterm, nu, logh, lap_Bterm,
                          lap_wterm, Omega);

  return exit_status;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath