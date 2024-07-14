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
#include "kadath.hpp"
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"
#include <math.h>
#include <sstream>

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
  
  if(std::isnan(bconfig.set(OUTER_SHELLS)))
    bconfig.set(OUTER_SHELLS) = 0;

  std::string kadath_filename = bconfig.space_filename();

	FILE* fin = fopen(kadath_filename.c_str(), "r") ;
	Space_bhns old_space(fin) ;
  Scalar old_conf  (old_space, fin) ;
  Scalar old_lapse (old_space, fin) ;
  Vector old_shift (old_space, fin) ;
  Scalar old_logh  (old_space, fin) ;
  Scalar old_phi   (old_space, fin) ;

	fclose(fin) ;

  const std::array<int, 2> old_adapted_doms{old_space.ADAPTEDNS, old_space.ADAPTEDBH};
  const std::array<int, 2> old_nuc_doms{old_space.NS, old_space.BH};

  const Domain_shell_outer_adapted* old_outer_adaptedNS =
        dynamic_cast<const Domain_shell_outer_adapted*>(old_space.get_domain(old_space.ADAPTEDNS));
  
  const Domain_shell_inner_adapted* old_inner_adaptedNS =
        dynamic_cast<const Domain_shell_inner_adapted*>(old_space.get_domain(old_space.ADAPTEDNS+1));

  const Domain_shell_outer_homothetic* old_bh_outer = 
        dynamic_cast<const Domain_shell_outer_homothetic*>(old_space.get_domain(old_space.ADAPTEDBH));

  std::cout << "Resolution of old space: "
  	<< old_space.get_domain(0)->get_nbr_points()(0) << " (r), "
    << old_space.get_domain(0)->get_nbr_points()(1) << " (theta), "
    << old_space.get_domain(0)->get_nbr_points()(2) << " (phi)" << std::endl;

  int ndim = 3;
	int ndom = old_space.get_nbr_domains() ;

  int res = bconfig(BIN_RES);

  if((res % 2) == 0){
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::cout << "Resolution of new space: "
  	<< res << " (r), "
    << res << " (theta), "
    << res - 1 << " (phi)" << std::endl;

  int type_coloc = old_space.get_type_base();

  // Update config vars
  // This control was mainly for testing
  if(!bconfig.control(USE_CONFIG_VARS)) {
    // Loop for ease of initialization
    for(int i = 0; i < 2; ++i){
      set_radius(old_nuc_doms[i]    , old_space, bconfig, RIN , i);
      set_radius(old_adapted_doms[i], old_space, bconfig, FIXED_R, i);
    }

    // update NS radii
    auto [ rmin, rmax ] = get_rmin_rmax(old_space, old_space.ADAPTEDNS);

    std::cout << "Rmin/max: " << std::endl
              << rmin << " " << rmax << std::endl;
    bconfig.set(RMID, BCO1) = rmin;
    bconfig.set(RIN , BCO1) = 0.5 * rmin;
    // end update NS radii

    // estimate how small the inner radius should be based on relation
    // between conformal factor and numerical radius.
    // see https://arxiv.org/pdf/0805.4192, eq(64)
    auto [ cmin, cmax ] = get_field_min_max(old_conf, old_space.ADAPTEDBH+1, INNER_BC);
    double conf_inner = cmin;
    double conf_i_sq  = conf_inner * conf_inner;
    double est_r_div2 = bconfig(MCH, BCO2) / conf_i_sq;
    bconfig.set(RIN, BCO2) =  est_r_div2;
    
    // this estimate is critical for calculating domain bounds
    // especially when attempting large changes in M_BH
    bconfig.set(RMID, BCO2) = 2 * est_r_div2; 
    
    rmax = (rmax > bconfig(RMID, BCO2)) ? rmax : bconfig(RMID, BCO2);
    bconfig.set(ROUT, BCO1) = (bconfig(DIST) / 2. - rmax) / 3. + rmax;
    bconfig.set(ROUT, BCO2) = bconfig(ROUT, BCO1);
  }// end updating config vars

  // create old radius scalar field
  Scalar old_space_radius(old_space);
  old_space_radius.annule_hard();
  for(int d = 0; d < ndom; ++d){
      old_space_radius.set_domain(d) = old_space.get_domain(d)->get_radius();
  }

  // Get outer domain radii
  old_space_radius.set_domain(old_space.ADAPTEDNS) = old_outer_adaptedNS->get_outer_radius();
  old_space_radius.set_domain(old_space.ADAPTEDBH) = old_bh_outer->get_outer_radius();

  old_space_radius.std_base();
  // end creating old radius scalar fields

  // setup bounds for creating new space
  std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
  std::vector<double> NS_bounds(3 + bconfig(NINSHELLS, BCO1) + bconfig(NSHELLS,BCO1));
  std::vector<double> BH_bounds(3 + bconfig(NSHELLS,BCO2));

  for(int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) * (1. + e * 0.25);

  set_NS_bounds(NS_bounds, bconfig, BCO1);
  set_BH_bounds(BH_bounds, bconfig, BCO2);
  
  // Set radius of the excision boundary to the current radius so that the solver
  // starts from the originial solution
  BH_bounds[1] = get_radius(old_space.get_domain(old_space.BH+1), OUTER_BC) ;
  // end setup bounds

  // print bounds to stdout - debugging only
  //std::cout << "Bounds:" << std::endl;
	//print_bounds("NS", NS_bounds);
	//print_bounds("BH", BH_bounds);
  //std::cout << std::endl;

  Space_bhns space (type_coloc, bconfig(DIST), NS_bounds, BH_bounds, out_bounds, bconfig(BIN_RES), bconfig(NINSHELLS, BCO1));
  Base_tensor basis(space, CARTESIAN_BASIS);

  const Domain_shell_inner_adapted* new_ns_inner = 
    dynamic_cast<const Domain_shell_inner_adapted*>(space.get_domain(space.ADAPTEDNS+1));
  const Domain_shell_outer_adapted* new_ns_outer = 
    dynamic_cast<const Domain_shell_outer_adapted*>(space.get_domain(space.ADAPTEDNS));

  // update BH fields to help with import
  update_adapted_field(old_conf , old_space.ADAPTEDBH+1, old_space.ADAPTEDBH, old_bh_outer, OUTER_BC);
  update_adapted_field(old_lapse, old_space.ADAPTEDBH+1, old_space.ADAPTEDBH, old_bh_outer, OUTER_BC);
  for(int i = 1; i < 4; ++i)
    update_adapted_field(old_shift.set(i), old_space.ADAPTEDBH+1, old_space.ADAPTEDBH, 
      old_bh_outer, OUTER_BC);

  update_adapted_field(old_phi, old_space.ADAPTEDNS, old_space.ADAPTEDNS+1, 
    old_inner_adaptedNS, INNER_BC);

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

  //import old fields into new space
  conf.import(old_conf);
  lapse.import(old_lapse);
  logh.import(old_logh);
  phi.import(old_phi);

  shift.set(1).import(old_shift.set(1));
  shift.set(2).import(old_shift.set(2));
  shift.set(3).import(old_shift.set(3));
  // end import

  // make sure there is no matter or vel.pot. outside of the star
  logh.set_domain(space.ADAPTEDNS+1).annule_hard();
  phi.set_domain(space.ADAPTEDNS+1).annule_hard();
  for(int d = space.BH; d < space.get_nbr_domains(); ++d){
    logh.set_domain(d).annule_hard();
    phi.set_domain(d).annule_hard();
    
    // make sure all fields are 0 inside the excision region
    if( (d >= space.BH) && (d < space.BH+2) ) {
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
}}