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
#include "kadath.hpp"
#include <math.h> 
#include <sstream>
#include "Configurator/config_bco.hpp"
#include "Solvers/co_solver_utils.hpp"
#include "bco_utilities.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
using space_t = Space_spheric_adapted;

template<typename config_t>
int ns_3d_xcts_regrid(config_t& bconfig, std::string outputfile) {
  int exit_status = 0;

  std::string kadath_filename = bconfig.space_filename();

  if(!fs::exists(kadath_filename)) {
    std::cerr << "File: " << kadath_filename << " not found.\n\n";
    std::_Exit(EXIT_FAILURE);
  }

	FILE* ff1 = fopen(kadath_filename.c_str(), "r") ;
	space_t old_space(ff1) ;
  Scalar old_conf(old_space, ff1) ;
	Scalar old_lapse(old_space, ff1) ;
	Vector old_shift(old_space, ff1) ;
	Scalar old_logh(old_space, ff1) ;
  fclose(ff1) ;

  std::cout << "Resolution of old space: "
    << old_space.get_domain(0)->get_nbr_points()(0) << " (r), "
    << old_space.get_domain(0)->get_nbr_points()(1) << " (theta), "
    << old_space.get_domain(0)->get_nbr_points()(2) << " (phi)" << std::endl;

  int ndim = 3;

  // get the adapted domain and cast it to its correct type to be able to call its member functions
  const Domain_shell_outer_adapted* old_outer_adapted =
        dynamic_cast<const Domain_shell_outer_adapted*>(old_space.get_domain(1));

  // setup a scalar field representing the old radius
  Scalar old_space_radius(old_space);
  old_space_radius = 0.;

  // get the radius from each domain
  for(int i = 0; i < old_space.get_nbr_domains(); ++i)
    old_space_radius.set_domain(i) = old_space.get_domain(i)->get_radius();
  // get the adapted radius of the adapted domain
  old_space_radius.set_domain(1) = old_outer_adapted->get_outer_radius();

  // define a standard decomposition, compatible with the parity of this field
  old_space_radius.std_base();
  //end setup old radius field

  // get the minimal and maximal radius from the adapted domain
  auto [r_min, r_max] = Kadath::bco_utils::get_rmin_rmax(old_space, 1);

  std::cout << "Rmin/max: " << r_min << " " << r_max << std::endl;

  // set new resolutions in each spatial dimension
  Dim_array res(ndim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = res(0);
  res.set(2) = res(0) - 1;

	// FIXME not sure if it's only about oddness...
  if(res(0) % 2 == 0 || res(2) % 2 != 0){
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::cout << "Resolution of new space: "
    << res(0) << " (r), "
    << res(1) << " (theta), "
    << res(2) << " (phi)" << std::endl;

  // get the type of the colocation points
  int type_coloc = old_space.get_type_base();

  // ignore domain radii scaling from config if needed
  if(!bconfig.control(USE_CONFIG_VARS)){
    bconfig.set(RIN)  = 0.5 * r_min;
    bconfig.set(ROUT) = 1.5 * r_max;
    bconfig.set(RMID) = r_max;
  }
  // end update config
  
  // setup radius bounds of the domains

  int ndom = 4 + bconfig(NSHELLS);
  std::vector<double> bounds(ndom-1);
  Kadath::bco_utils::set_NS_bounds(bounds, bconfig);
  
  Kadath::bco_utils::print_bounds("New bounds: ", bounds);

  // get origin of nucleus domain
  Point center = old_space.get_domain(0)->get_center();

  // initialize space with new resolution and domain decomposition
  Space_spheric_adapted space(type_coloc, center, res, bounds);
  Base_tensor basis(space, CARTESIAN_BASIS);

  // get adapted domains to update the radius
  const Domain_shell_outer_adapted* new_outer_adapted = dynamic_cast<const Domain_shell_outer_adapted*>(space.get_domain(1));
	const Domain_shell_inner_adapted* new_inner_adapted = dynamic_cast<const Domain_shell_inner_adapted*>(space.get_domain(2));

  // update adapted domain mapping
  Kadath::bco_utils::interp_adapted_mapping(new_outer_adapted, 1, old_space_radius);
  Kadath::bco_utils::interp_adapted_mapping(new_inner_adapted, 1, old_space_radius);

  // setup new fields
  // initialize to one or zero first
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

  // end setup new fields
  
  // import data from fields in the old space
  conf.import(old_conf);
  lapse.import(old_lapse);
  logh.import(old_logh);

  shift.set(1).import(old_shift.set(1));
  shift.set(2).import(old_shift.set(2));
  shift.set(3).import(old_shift.set(3));

  // end import old fields

  // enforce spectral decomposition compatible with the parities
  lapse.std_base();
  conf.std_base();
  logh.std_base();
  shift.std_base();
  
  // output data  
  bconfig.set_filename(outputfile);
  Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);

  return EXIT_SUCCESS;
}
/** @}*/
}}