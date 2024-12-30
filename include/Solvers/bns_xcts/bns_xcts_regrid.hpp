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

#include "Configurator/config_binary.hpp"
#include "Solvers/fuka_syst/fuka_syst.hpp"
#include "bco_utilities.hpp"
#include "Solvers/bco_solver_utils.hpp"
#include "kadath_bin_ns.hpp"

/**
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

using bin_space_t = Space_bin_ns;

template <typename config_t>
int bns_xcts_regrid(config_t& bconfig, std::string output_fname) {
  int exit_status = EXIT_SUCCESS;
  using namespace ::Kadath::bco_utils;

  // file containing KADATH fields must have same name as config file with only
  // the extension being different
  std::string kadath_filename = bconfig.space_filename();
  if (!fs::exists(kadath_filename)) {
    std::cerr << "File: " << kadath_filename << " not found.\n\n";
    std::_Exit(EXIT_FAILURE);
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

  std::array<const Domain_shell_outer_adapted*, 2> old_outer_adapted;
  std::array<const Domain_shell_inner_adapted*, 2> old_inner_adapted;

  for (int i = 0; i < 2; ++i) {
    int const d = old_adapted_doms[i];
    old_outer_adapted[i] = dynamic_cast<const Domain_shell_outer_adapted*>(
        old_space.get_domain(d));
    old_inner_adapted[i] = dynamic_cast<const Domain_shell_inner_adapted*>(
        old_space.get_domain(d + 1));
  }

  int res = bconfig.set(BIN_PARAMS::BIN_RES);
  if ((res % 2) == 0) {
    std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)"
              << std::endl;
    std::_Exit(EXIT_FAILURE);
  }
  int ndim = 3;
  int ndom = old_space.get_nbr_domains();

  int type_coloc = old_space.get_type_base();

  // start Update config vars
  update_config_NS_radii(old_space, bconfig, old_space.ADAPTED1, NODES::BCO1);
  update_config_NS_radii(old_space, bconfig, old_space.ADAPTED2, NODES::BCO2);
  double r_max_tot = std::max(bconfig(BCO_PARAMS::RMID, BCO1), bconfig(BCO_PARAMS::RMID, BCO2));

  const double rout_sep_est =
      (bconfig(BIN_PARAMS::DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  const double rout_max_est = gold_ratio * r_max_tot;
  bconfig.set(BCO_PARAMS::ROUT, NODES::BCO1) = rout_sep_est;
      // (rout_sep_est > rout_max_est) ? rout_max_est : rout_sep_est;
  bconfig.set(BCO_PARAMS::ROUT, NODES::BCO2) =
      bconfig(BCO_PARAMS::ROUT, NODES::BCO1);
  // end updating config vars

  // create old radius scalar field
  Scalar old_space_radius(old_space);
  old_space_radius.annule_hard();

  for (int d = 0; d < ndom; ++d) {
    old_space_radius.set_domain(d) = old_space.get_domain(d)->get_radius();
  }

  // Get adapted domain radii
  for (int i = 0; i < 2; ++i) {
    int const dom = old_adapted_doms[i];
    old_space_radius.set_domain(dom) = old_outer_adapted[i]->get_outer_radius();
  }
  old_space_radius.std_base();
  // end create old radius scalar fields

  std::vector<int> ns_interior_doms{FUKA_Syst_tools::vector_of_domains(
      old_space.NS1, old_space.ADAPTED1)};
  std::vector<int> exclusion_doms{old_space.NS2, old_space.ADAPTED1};
  // concat domain lists together
  std::for_each(ns_interior_doms.rbegin(), ns_interior_doms.rend(),
                [&exclusion_doms](auto e) {
                  auto it = exclusion_doms.begin();
                  exclusion_doms.insert(it, e);
                });
  Scalar drPsi(compute_drPsi(old_space, old_conf,
                           Metric_flat(old_space, old_shift.get_basis()),
                           exclusion_doms, old_space.OUTER));

  std::vector<double> out_bounds(1 + bconfig(BIN_PARAMS::OUTER_SHELLS));
  // set reasonable radii to each stellar domain
  std::vector<double> NS1_bounds{
      set_arb_boundsv3(bconfig, drPsi, old_space.ADAPTED1 + 1, NODES::BCO1)};
  std::vector<double> NS2_bounds{
      set_arb_boundsv3(bconfig, drPsi, old_space.ADAPTED2 + 1, NODES::BCO2)};

  // space needs to be able fixed to add shells
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(BIN_PARAMS::REXT) * (1. + e * 0.25);

  std::cout << "Bounds:" << std::endl;
  print_bounds("NS1-bounds", NS1_bounds);
  print_bounds("NS2-bounds", NS2_bounds);
  print_bounds("Outer", out_bounds);

  bin_space_t space(type_coloc, bconfig(BIN_PARAMS::DIST), NS1_bounds,
                    NS2_bounds, out_bounds, res);

  std::cout << "Resolution of old space: ";
  print_constant_space_resolution(old_space);

  std::cout << "Resolution of new space: ";
  print_constant_space_resolution(space);

  std::cout << "\nOld bounds:" << std::endl;
  print_bounds_from_space(old_space);

  std::cout << "New bounds:" << std::endl;
  print_bounds_from_space(space);

  ndom = space.get_nbr_domains();

  Base_tensor basis(space, CARTESIAN_BASIS);

  const std::array<int, 2> new_adapted_doms{space.ADAPTED1, space.ADAPTED2};
  const std::array<int, 2> new_nuc_doms{space.NS1, space.NS2};
  const std::array<double, 2> xc{space.get_domain(space.NS1)->get_center()(1),
                                 space.get_domain(space.NS2)->get_center()(1)};

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

    // Updated mapping for NS
    interp_adapted_mapping(new_inner_adapted[i], dom, old_space_radius);
    interp_adapted_mapping(new_outer_adapted[i], dom, old_space_radius);

    // Interpolate old_phi field outside of the star for import
    update_adapted_field(old_phi, dom, dom + 1, old_inner_adapted[i], INNER_BC);
  }

  std::cout << "xc1: " << xc[0] << std::endl;
  std::cout << "xc2: " << xc[1] << std::endl;

  // start Create new fields
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
  // end Create new fields

  // start Import fields
  conf.import(old_conf);
  lapse.import(old_lapse);
  logh.import(old_logh);
  phi.import(old_phi);

  shift.set(1).import(old_shift.set(1));
  shift.set(2).import(old_shift.set(2));
  shift.set(3).import(old_shift.set(3));
  // end Import fields

  // start - set matter and velocity potential to zero outside stars
  for (auto& dom : new_adapted_doms) {
    phi.set_domain(dom + 1).annule_hard();
    logh.set_domain(dom + 1).annule_hard();
  }

  for (int i = space.OUTER; i < ndom; ++i) {
    phi.set_domain(i).annule_hard();
    logh.set_domain(i).annule_hard();
  }
  // end - set phi,logh to zero outside

  lapse.std_base();
  conf.std_base();
  logh.std_base();
  shift.std_base();
  phi.std_base();

  bconfig.set_filename(output_fname);
  save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
  return exit_status;
}
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath