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
 *
 * This implementation is designed to receive a NS to start from - i.e.
 * we are not creating a star from scratch.  This is useful for 
 * running sequences of NSs.
*/
#include <string>
#include <typeinfo>
#include <utility>
#include "Configurator/config_binary.hpp"
#include "EOS/EOS.hh"
#include "bco_utilities.hpp"
#include "kadath.hpp"
#include "kadath_adapted_bh.hpp"
const double M2km = 1.4769994423016508;
double distkm = 45.;

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

using config_t = kadath_config_boost<BIN_INFO>;

// forward declarations
template <typename eos_t>
void setup_binbhns(kadath_config_boost<BCO_NS_INFO>& NSconfig,
                   kadath_config_boost<BCO_BH_INFO>& BHconfig);

void setup_binbhns_config(config_t& bconfig);

template <typename eos_t, typename space_t>
void bin_config_import_NS(kadath_config_boost<BCO_NS_INFO>& NSconfig,
                          config_t& bconfig,
                          const space_t& space);

template <typename space_t>
void bin_config_import_BH(kadath_config_boost<BCO_BH_INFO>& BHconfig,
                          config_t& bconfig,
                          const space_t& space,
                          Scalar& conf);

// end forward declarations

int main(int argc, char** argv) {
  // expecting one 3D NS solution and one 3D BH input
  if (argc < 3) {
    std::cerr << "Usage: ./setup /<path>/<str: NS ID basename>.info "
                 "/<path>/<str: BH ID basename>.info"
                 "<float: separation distance in KM>"
              << std::endl;
    std::cerr << "e.g. ./setup converged.NS.9.info converged.BH.9.info 45"
              << std::endl;
    std::cerr << "the separation is an optional parameter." << std::endl;
    std::_Exit(EXIT_FAILURE);
  }
  if (argc > 3) {
    distkm = std::atof(argv[3]);
    std::cout << "Separation Distance set to: " << distkm << "km" << std::endl;
  } else
    std::cout << "Separation Distance set to default " << distkm << "km"
              << std::endl;

  // readin NS and BH configs
  std::string bhfilename{argv[1]};
  std::string nsfilename{argv[2]};
  kadath_config_boost<BCO_NS_INFO> NSconfig(nsfilename);
  kadath_config_boost<BCO_BH_INFO> BHconfig(bhfilename);

  // setup eos and update central density
  const double h_cut = NSconfig.eos<double>(HCUT);
  const std::string eos_file = NSconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = NSconfig.eos<std::string>(EOSTYPE);

  if (eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut);
    setup_binbhns<eos_t>(NSconfig, BHconfig);
  } else if (eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (NSconfig.eos<int>(INTERP_PTS) == 0)
                               ? 2000
                               : NSconfig.eos<int>(INTERP_PTS);

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut, interp_pts);
    setup_binbhns<eos_t>(NSconfig, BHconfig);
  }

  return EXIT_SUCCESS;
}

/**
 * setup_binbhns
 * Setup and initialization of BHNS config, numerical space and fields.
 *
 * @tparam nsconfig_t BNS config file type
 * @tparam bhconfig_t BH config file type
 * @tparam bconfig_t BHNS config file type
 * @param[input] BNSconfig BNS config file containing needed NS parameters
 * @param[input] BHconfig BH config file containing needed BH parameters
 * @param[input] bconfig BHNS config file being updated
 */
template <typename eos_t>
void setup_binbhns(kadath_config_boost<BCO_NS_INFO>& NSconfig,
                   kadath_config_boost<BCO_BH_INFO>& BHconfig) {
  // Initialize Config File
  std::string ifilename = "./initbin.info";
  kadath_config_boost<BIN_INFO> bconfig;
  bconfig.set_filename(ifilename);
  bconfig.initialize_binary({"ns", "bh"});
  // open previous ns solution
  std::string nsspaceinf = NSconfig.space_filename();
  FILE* ff1 = fopen(nsspaceinf.c_str(), "r");
  Space_spheric_adapted nsspacein(ff1);
  Scalar nsconf(nsspacein, ff1);
  Scalar nslapse(nsspacein, ff1);
  Vector nsshift(nsspacein, ff1);
  Scalar nslogh(nsspacein, ff1);
  Scalar nsphi(nsspacein, ff1);
  fclose(ff1);
  // end opening bns solution

  // Update HC before passing BNSconfig to setup_binbhns_config
  NSconfig.set(HC) = std::exp(bco_utils::get_boundary_val(0, nslogh));
  bin_config_import_NS<eos_t>(NSconfig, bconfig, nsspacein);

  // obtain adapted NS shells for radius information and copying adapted mapping later
  const Domain_shell_outer_adapted* old_outer_adapted1 =
      dynamic_cast<const Domain_shell_outer_adapted*>(nsspacein.get_domain(1));
  const Domain_shell_inner_adapted* old_inner_adapted1 =
      dynamic_cast<const Domain_shell_inner_adapted*>(nsspacein.get_domain(2));

  // setup radius field - needed for copying the adapted domain mappings.
  Scalar old_space_radius(nsspacein);
  old_space_radius.annule_hard();

  int ndominns = nsspacein.get_nbr_domains();

  for (int d = 0; d < ndominns; ++d)
    old_space_radius.set_domain(d) = nsspacein.get_domain(d)->get_radius();

  old_space_radius.set_domain(1) = old_outer_adapted1->get_outer_radius();
  old_space_radius.std_base();
  // end setup radius field

  // open old BH solution
  std::string bhspaceinf = BHconfig.space_filename();
  FILE* ff2 = fopen(bhspaceinf.c_str(), "r");
  Space_adapted_bh bhspacein(ff2);
  Scalar bhconf(bhspacein, ff2);
  Scalar bhlapse(bhspacein, ff2);
  Vector bhshift(bhspacein, ff2);
  fclose(ff2);
  // end open BH solution

  bin_config_import_BH(BHconfig, bconfig, bhspacein, bhconf);
  setup_binbhns_config(bconfig);

  std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
  std::vector<double> NS_bounds(3 + bconfig(NSHELLS, BCO1));
  std::vector<double> BH_bounds(3 + bconfig(NSHELLS, BCO2));

  bco_utils::set_BH_bounds(BH_bounds, bconfig, BCO2, true);

  //for out_bounds.size > 1 - add equi-distance shells
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) + e * 0.25 * bconfig(REXT);

  bco_utils::set_NS_bounds(NS_bounds, bconfig, BCO1);

  std::cout << bconfig << std::endl;

  bco_utils::print_bounds("NS-bounds", NS_bounds);
  bco_utils::print_bounds("BH-bounds", BH_bounds);
  bco_utils::print_bounds("outer-bounds", out_bounds);

  //Setup actual space
  int typer = CHEB_TYPE;
  Space_bhns space(typer, bconfig(DIST), NS_bounds, BH_bounds, out_bounds,
                   bconfig(BIN_RES));
  Base_tensor basis(space, CARTESIAN_BASIS);

  const Domain_shell_inner_adapted* new_ns_inner =
      dynamic_cast<const Domain_shell_inner_adapted*>(space.get_domain(2));
  const Domain_shell_outer_adapted* new_ns_outer =
      dynamic_cast<const Domain_shell_outer_adapted*>(space.get_domain(1));

  const Domain_shell_outer_homothetic* old_bh_outer =
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          bhspacein.get_domain(1));

  //Update BH fields based to help with interpolation later
  bco_utils::update_adapted_field(bhconf, 2, 1, old_bh_outer, OUTER_BC);
  bco_utils::update_adapted_field(bhlapse, 2, 1, old_bh_outer, OUTER_BC);

  //Updated mapping for NS
  bco_utils::interp_adapted_mapping(new_ns_inner, 1, old_space_radius);
  bco_utils::interp_adapted_mapping(new_ns_outer, 1, old_space_radius);

  double xc1 = bco_utils::get_center(space, space.NS);
  double xc2 = bco_utils::get_center(space, space.BH);

  std::cout << "xc1: " << xc1 << std::endl;
  std::cout << "xc2: " << xc2 << std::endl;

  if (NSconfig.set_field(PHI) == true)
    bco_utils::update_adapted_field(nsphi, 1, 2, old_inner_adapted1, INNER_BC);

  Scalar logh(space);
  logh.annule_hard();
  logh.std_base();

  Scalar conf(space);
  conf.annule_hard();

  Scalar lapse(space);
  lapse.annule_hard();

  //After much testing, it was found that discarding the shift provided a better initial guess
  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();

  Scalar phi(space);
  phi.annule_hard();

  const double ns_overw4 = bco_utils::set_decay(bconfig, BCO1);
  const double bh_overw4 = bco_utils::set_decay(bconfig, BCO2);

  std::cout << "WeightNS: " << bconfig(DECAY, BCO1) << ", "
            << "WeightBH: " << bconfig(DECAY, BCO2) << std::endl;
  int ndom = space.get_nbr_domains();
  for (int dom = 0; dom < ndom; dom++) {
    if (dom != space.BH && dom != space.BH + 1) {
      Index new_pos(space.get_domain(dom)->get_nbr_points());

      do {
        double x = space.get_domain(dom)->get_cart(1)(new_pos);
        double y = space.get_domain(dom)->get_cart(2)(new_pos);
        double z = space.get_domain(dom)->get_cart(3)(new_pos);

        // define a point shifted suitably to the stellar centers in the binary
        Point absol1(3);
        absol1.set(1) = (x - xc1);
        absol1.set(2) = y;
        absol1.set(3) = z;
        double r2 = y * y + z * z;
        double r2_1 = (x - xc1) * (x - xc1) + r2;
        double r4_1 = r2_1 * r2_1;
        double r4_overw4_1 = r4_1 * ns_overw4;
        double decay_1 = std::exp(-r4_overw4_1);

        Point absol2(3);
        absol2.set(1) = (x - xc2);
        absol2.set(2) = y;
        absol2.set(3) = z;
        double r2_2 = (x - xc2) * (x - xc2) + r2;
        double r4_2 = r2_2 * r2_2;
        double r4_overw4_2 = r4_2 * ns_overw4;
        double decay_2 = std::exp(-r4_overw4_2);

        if (dom != ndom - 1) {
          conf.set_domain(dom).set(new_pos) =
              1. + decay_1 * (nsconf.val_point(absol1) - 1.) +
              decay_2 * (bhconf.val_point(absol2) - 1.);

          lapse.set_domain(dom).set(new_pos) =
              1. + decay_1 * (nslapse.val_point(absol1) - 1.) +
              decay_2 * (bhlapse.val_point(absol2) - 1.);

          logh.set_domain(dom).set(new_pos) =
              decay_1 * nslogh.val_point(absol1);

          if (NSconfig.set_field(PHI) == true)
            phi.set_domain(dom).set(new_pos) =
                decay_1 * nsphi.val_point(absol1);

          for (int i = 1; i <= 3; i++)
            shift.set(i).set_domain(dom).set(new_pos) =
                decay_1 * nsshift(i).val_point(absol1) +
                decay_2 * bhshift(i).val_point(absol2);
        } else {
          conf.set_domain(dom).set(new_pos) = 1.;
          lapse.set_domain(dom).set(new_pos) = 1.;
        }
      } while (new_pos.inc());
    }
  }
  // Want to make sure there is no matter around the BH
  for (int d = space.BH; d < space.OUTER; ++d) {
    logh.set_domain(d).annule_hard();
    phi.set_domain(d).annule_hard();
  }
  for (int d = space.BH; d < space.BH + 2; ++d) {
    conf.set_domain(d).annule_hard();
    lapse.set_domain(d).annule_hard();
  }

  conf.std_base();
  lapse.std_base();
  logh.std_base();
  shift.std_base();
  phi.std_base();

  bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
}

/**
 * setup_binbhns_config
 * Setup and initialization of basic binary details, based on the imported objects.
 *
 * @param[input] bconfig BHNS config file being updated
 */
void setup_binbhns_config(config_t& bconfig) {

  // setup binary parameters
  bconfig.set(DIST) = distkm / M2km;
  bconfig.set(REXT) = 2 * bconfig(DIST);
  bconfig.set(CHECKPOINT) = 0;
  bconfig.set(QPIG) = bconfig(BCO_QPIG, BCO1);
  bconfig.set(Q) = bconfig(MADM, BCO1) / bconfig(MCH, BCO2);
  bconfig.set(COM) = bco_utils::com_estimate(bconfig(DIST), bconfig(MADM, BCO1),
                                             bconfig(MCH, BCO2));
  bconfig.set(COM) =
      (bconfig(MADM, BCO1) > bconfig(MCH, BCO2)) ? bconfig(COM) : -bconfig(COM);
  bconfig.set(COMY) = 0.;
  bconfig.set(OUTER_SHELLS) = 0;
  bconfig.set(BIN_RES) = bconfig(BCO_RES, BCO1);

  // obtain 3PN estimate for the global, orbital omega
  bco_utils::KadathPNOrbitalParams(bconfig, bconfig(MADM, BCO1),
                                   bconfig(MCH, BCO2));

  // delete ADOT, this can always be recalculated during
  // the eccentricity reduction stage
  bconfig.set(ADOT) = std::nan("1");

  // update to ensure consistent ROUT
  double r_max_tot = std::max(bconfig(RMID, BCO1), bconfig(RMID, BCO2));
  bconfig.set(ROUT, BCO1) = (bconfig(DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  bconfig.set(ROUT, BCO2) = (bconfig(DIST) / 2. - r_max_tot) / 3. + r_max_tot;

  // catalog fields that are implemented
  bconfig.set_field(CONF) = true;
  bconfig.set_field(LAPSE) = true;
  bconfig.set_field(SHIFT) = true;
  bconfig.set_field(LOGH) = true;
  bconfig.set_field(PHI) = true;

  // set stages
  bconfig.set_stage(TOTAL) = false;
  bconfig.set_stage(TOTAL_BC) = true;
  bconfig.control(FIXED_GOMEGA) = true;
}

/**
 * bin_config_import_NS
 * update BHNS config based on imported NS specific parameters
 *
 * @tparam nsconfig_t config file type containing NS parameters
 * @tparam space_t space file type of imported BNS
 * @param[input] BNSconfig configuration file containing NS parameters
 * @param[input] bconfig bhns configuration file being updated
 * @param[input] BNSspace BNS space to pull radius information from 
 */
template <typename eos_t, typename space_t>
void bin_config_import_NS(kadath_config_boost<BCO_NS_INFO>& NSconfig,
                          config_t& bconfig,
                          const space_t& space) {
  bco_utils::update_config_NS_radii(space, NSconfig, NSconfig(NSHELLS) + 1);
  NSconfig.set(NC) = EOS<eos_t, DENSITY>::get(NSconfig(HC));

  //Copy NS Parameters
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO1) = NSconfig.set(i);
  }

  //Copy  EOS Parameters
  for (int i = 0; i < NUM_EOS_PARAMS; ++i) {
    bconfig.set_eos(i, BCO1) = NSconfig.set_eos(i);
  }
}

/**
 * bin_config_import_BH
 * update BHNS config based on imported BH specific parameters
 *
 * @tparam bhconfig_t config file type containing BH parameters
 * @tparam space_t space file type of imported BH 
 * @param[input] BHconfig configuration file containing BH parameters
 * @param[input] bconfig bhns configuration file being updated
 * @param[input] space BH space to pull radius information from 
 */
template <typename space_t>
void bin_config_import_BH(kadath_config_boost<BCO_BH_INFO>& BHconfig,
                          config_t& bconfig,
                          const space_t& space,
                          Scalar& conf) {
  bco_utils::update_config_BH_radii(space, BHconfig, 1, conf);
  // Copy BH parameters
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO2) = BHconfig.set(i);
  }
  bconfig.set(NSHELLS, BCO2) = 1;
}
