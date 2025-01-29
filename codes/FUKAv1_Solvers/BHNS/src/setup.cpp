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
#include "mpi.h"
const double M2km = 1.4769994423016508;
double distkm = 45.;

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

using config_t = kadath_config_boost<BIN_INFO>;

// forward declarations
template <typename nsconfig_t, typename bhconfig_t>
void setup_binbhns(nsconfig_t& BNSconfig,
                   bhconfig_t& BHconfig,
                   config_t& bconfig);

void setup_binbhns_config(config_t& bconfig);

template <typename nsconfig_t, typename space_t>
void bin_config_import_NS(nsconfig_t& BNSconfig,
                          config_t& bconfig,
                          const space_t& space);

template <typename bhconfig_t, typename space_t>
void bin_config_import_BH(bhconfig_t& BHconfig,
                          config_t& bconfig,
                          const space_t& space);

// end forrward declarations

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (argc < 3) {
    std::cerr << "Missing input config files *.info\n"
              << "Usage: setup <bns.info> <bh.info>\n";
    std::_Exit(EXIT_FAILURE);
  }

  // readin BNS and BH configs
  std::string nsfilename{argv[1]};
  std::string bhfilename{argv[2]};
  kadath_config_boost<BIN_INFO> BNSconfig(nsfilename);
  kadath_config_boost<BCO_BH_INFO> BHconfig(bhfilename);

  // initialize bhns config file
  std::string ifilename = "./initbin.info";
  kadath_config_boost<BIN_INFO> bconfig;
  bconfig.set_filename(ifilename);
  bconfig.initialize_binary({"ns", "bh"});

  //always start at low res
  bconfig.set(BIN_RES) = 9;
  bconfig.set(DIST) = BNSconfig(DIST);
  setup_binbhns(BNSconfig, BHconfig, bconfig);

  MPI_Finalize();

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
template <typename nsconfig_t, typename bhconfig_t>
void setup_binbhns(nsconfig_t& BNSconfig,
                   bhconfig_t& BHconfig,
                   config_t& bconfig) {
  // open previous bns solution
  std::string nsspaceinf = BNSconfig.space_filename();
  FILE* ff1 = fopen(nsspaceinf.c_str(), "r");
  Space_bin_ns nsspacein(ff1);
  Scalar nsconf(nsspacein, ff1);
  Scalar nslapse(nsspacein, ff1);
  Vector nsshift(nsspacein, ff1);
  Scalar nslogh(nsspacein, ff1);
  Scalar nsphi(nsspacein, ff1);
  fclose(ff1);
  // end opening bns solution

  //delete secondary star
  for (int d = nsspacein.NS2; d < nsspacein.OUTER; ++d) {
    nsconf.set_domain(d) = 1.;
    nslapse.set_domain(d) = 1.;
  }

  //Update HC before passing BNSconfig to setup_binbhns_config
  BNSconfig.set(HC) = std::exp(bco_utils::get_boundary_val(0, nslogh));

  bin_config_import_NS(BNSconfig, bconfig, nsspacein);

  int typer = CHEB_TYPE;

  // obtain adapted NS shells for radius information and copying adapted mapping later
  const Domain_shell_outer_adapted* old_outer_adapted1 =
      dynamic_cast<const Domain_shell_outer_adapted*>(
          nsspacein.get_domain(nsspacein.ADAPTED1));
  const Domain_shell_inner_adapted* old_inner_adapted1 =
      dynamic_cast<const Domain_shell_inner_adapted*>(
          nsspacein.get_domain(nsspacein.ADAPTED1 + 1));

  // setup radius field - needed for copying the adapted domain mappings.
  Scalar old_space_radius(nsspacein);
  old_space_radius.annule_hard();

  int ndominns = nsspacein.get_nbr_domains();

  for (int d = 0; d < ndominns; ++d)
    old_space_radius.set_domain(d) = nsspacein.get_domain(d)->get_radius();

  old_space_radius.set_domain(nsspacein.ADAPTED1) =
      old_outer_adapted1->get_outer_radius();
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

  // setup BH config settings based on old solution
  bin_config_import_BH(BHconfig, bconfig, bhspacein);

  setup_binbhns_config(bconfig);
  std::cout << bconfig << "\n\n";
  bconfig.set(GOMEGA) = BNSconfig(GOMEGA);
  std::vector<double> out_bounds(1);
  std::vector<double> NS_bounds(3 + bconfig(NSHELLS, BCO1));
  std::vector<double> BH_bounds(3 + bconfig(NSHELLS, BCO2));

  bco_utils::set_BH_bounds(BH_bounds, bconfig, BCO2, true);

  //for out_bounds.size > 1 - add equi-distance shells
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) + e * 0.25 * bconfig(REXT);

  bco_utils::set_NS_bounds(NS_bounds, bconfig, BCO1);
  bco_utils::print_bounds("NS-bounds", NS_bounds);
  bco_utils::print_bounds("BH-bounds", BH_bounds);
  bco_utils::print_bounds("outer-bounds", out_bounds);
  //Setup actual space
  Space_bhns space(typer, bconfig(DIST), NS_bounds, BH_bounds, out_bounds,
                   bconfig(BIN_RES));
  Base_tensor basis(space, CARTESIAN_BASIS);

  const Domain_shell_inner_adapted* new_ns_inner =
      dynamic_cast<const Domain_shell_inner_adapted*>(
          space.get_domain(space.ADAPTEDNS + 1));
  const Domain_shell_outer_adapted* new_ns_outer =
      dynamic_cast<const Domain_shell_outer_adapted*>(
          space.get_domain(space.ADAPTEDNS));

  const Domain_shell_outer_homothetic* old_bh_outer =
      dynamic_cast<const Domain_shell_outer_homothetic*>(
          bhspacein.get_domain(1));

  //Update BH fields based to help with interpolation later
  bco_utils::update_adapted_field(bhconf, 2, 1, old_bh_outer, OUTER_BC);
  bco_utils::update_adapted_field(bhlapse, 2, 1, old_bh_outer, OUTER_BC);

  //Updated mapping for NS
  bco_utils::interp_adapted_mapping(new_ns_inner, nsspacein.ADAPTED1,
                                    old_space_radius);
  bco_utils::interp_adapted_mapping(new_ns_outer, nsspacein.ADAPTED1,
                                    old_space_radius);

  double xc1 = bco_utils::get_center(space, space.NS);
  double xc2 = bco_utils::get_center(space, space.BH);

  std::cout << "xc1: " << xc1 << std::endl;
  std::cout << "xc2: " << xc2 << std::endl;

  /** for BNS import, we don't shift the coordinates
    * Note BNS must (currently) be the same separation as BHNS space 
    * */
  xc1 = 0;

  //Only done for BNS since phi is not applicable to single star
  bco_utils::update_adapted_field(nsphi, nsspacein.ADAPTED1,
                                  nsspacein.ADAPTED1 + 1, old_inner_adapted1,
                                  INNER_BC);

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
  int ndom = space.get_nbr_domains();
  for (int dom = 0; dom < ndom; dom++) {
    if (dom != space.BH && dom != space.BH + 1) {
      Index new_pos(space.get_domain(dom)->get_nbr_points());

      do {
        double x = space.get_domain(dom)->get_cart(1)(new_pos);
        double y = space.get_domain(dom)->get_cart(2)(new_pos);
        double z = space.get_domain(dom)->get_cart(3)(new_pos);

        Point absol1(3);
        absol1.set(1) = x - xc1;
        absol1.set(2) = y;
        absol1.set(3) = z;

        Point absol2(3);
        absol2.set(1) = x - xc2;
        absol2.set(2) = y;
        absol2.set(3) = z;
        if (dom != ndom - 1) {
          conf.set_domain(dom).set(new_pos) =
              nsconf.val_point(absol1) * bhconf.val_point(absol2);
          lapse.set_domain(dom).set(new_pos) =
              nslapse.val_point(absol1) * bhlapse.val_point(absol2);
          phi.set_domain(dom).set(new_pos) = nsphi.val_point(absol1);
          logh.set_domain(dom).set(new_pos) = nslogh.val_point(absol1);
        } else {
          conf.set_domain(dom).set(new_pos) = 1.;
          lapse.set_domain(dom).set(new_pos) = 1.;
          logh.set_domain(dom).annule_hard();
          phi.set_domain(dom).annule_hard();
        }
      } while (new_pos.inc());
    }
  }
  /* Want to make sure there is no matter or vel.pot. around the BH initially - safety when importing from BNS*/
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
 * Setup and initialization of basic binary details, independent of the imported objects.
 *
 * @param[input] bconfig BHNS config file being updated
 */
void setup_binbhns_config(config_t& bconfig) {
  // setup binary parameters
  bconfig.set(REXT) = 2 * bconfig(DIST);
  bconfig.set(CHECKPOINT) = 0;
  bconfig.set(QPIG) = bconfig(BCO_QPIG, BCO1);
  bconfig.set(Q) = 1;
  bconfig.set(COM) = 0.;
  bconfig.set(COMY) = 0.;
  bconfig.set(OUTER_SHELLS) = 0;

  // catalog fields that are implemented
  bconfig.set_field(CONF) = true;
  bconfig.set_field(LAPSE) = true;
  bconfig.set_field(SHIFT) = true;
  bconfig.set_field(LOGH) = true;
  bconfig.set_field(PHI) = true;

  // set stages
  bconfig.set_stage(TOTAL) = true;
  bconfig.set_stage(TOTAL_BC) = true;
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
template <typename nsconfig_t, typename space_t>
void bin_config_import_NS(nsconfig_t& BNSconfig,
                          config_t& bconfig,
                          const space_t& BNSspace) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //Copy NS Parameters
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO1) = BNSconfig.set(i, BCO1);
  }

  //Copy  EOS Parameters
  for (int i = 0; i < NUM_EOS_PARAMS; ++i) {
    bconfig.set_eos(i, BCO1) = BNSconfig.set_eos(i, BCO1);
  }

  // setup eos to update central density
  const double h_cut = bconfig.eos<double>(HCUT, BCO1);

  const std::string eos_file = bconfig.eos<std::string>(EOSFILE, BCO1);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);

  if (eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut);
    bconfig.set(NC, BCO1) = EOS<eos_t, DENSITY>::get(bconfig(HC, BCO1));
  } else if (eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS, BCO1) == 0)
                               ? 2000
                               : bconfig.eos<int>(INTERP_PTS, BCO1);

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut, interp_pts);
    bconfig.set(NC, BCO1) = EOS<eos_t, DENSITY>::get(bconfig(HC, BCO1));
  } else {
    if (rank == 0)
      std::cerr << eos_type << " is not recognized.\n";
    std::_Exit(EXIT_FAILURE);
  }
  // end eos setup

  auto [r_min, r_max] =
      bco_utils::get_rmin_rmax(BNSspace, 1 + BNSconfig(NSHELLS, BCO1));
  bconfig.set(RIN, BCO1) = 0.5 * r_min;
  bconfig.set(RMID, BCO1) = r_max;
  bconfig.set(ROUT, BCO1) = (bconfig(DIST) / 2. - r_max) / 3. + r_max;
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
template <typename bhconfig_t, typename space_t>
void bin_config_import_BH(bhconfig_t& BHconfig,
                          config_t& bconfig,
                          const space_t& space) {
  // Copy BH parameters
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO2) = BHconfig.set(i);
  }
  bconfig.set(RIN, BCO2) = bconfig(MCH, BCO2) * bco_utils::invpsisq;

  // Make sure router is the same for both - set to the larger of the two.
  if (bconfig(ROUT, BCO2) > bconfig(ROUT, BCO1))
    bconfig.set(ROUT, BCO1) = bconfig(ROUT, BCO2);
  else
    bconfig.set(ROUT, BCO2) = bconfig(ROUT, BCO1);

  // update config FIXED_R based on AH Surface radius
  bco_utils::set_radius(1, space, bconfig, FIXED_R, BCO2);
}
