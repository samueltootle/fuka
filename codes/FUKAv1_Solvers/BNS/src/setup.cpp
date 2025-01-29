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
 *
 * This implementation is designed to receive a NS to start from - i.e.
 * we are not creating a star from scratch.  This is useful for 
 * running sequences of NSs.
*/
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>
#include "Configurator/config_binary.hpp"
#include "EOS/EOS.hh"
#include "bco_utilities.hpp"
#include "kadath.hpp"
#include "mpi.h"

// conversion from solar mass to km
const double M2km = 1.4769994423016508;
// distance in km, potential input parameter
// declared global to be used in a few places easily
double distkm = 45.;

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

// forward declarations
template <typename config_t, typename space_t>
kadath_config_boost<BIN_INFO> setup_bin_config(config_t& NSconfig,
                                               const space_t& spacein1);

template <typename config_t>
void setup_3d(config_t NSconfig);

// end forrward declarations

int main(int argc, char** argv) {
  // initialize MPI
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // expecting one 3D NS solution as input
  if (argc < 2) {
    std::cerr << "Usage: ./setup /<path>/<str: NS ID basename>.info "
                 "<float: separation distance in KM>"
              << std::endl;
    std::cerr << "e.g. ./setup converged.NS.9.info 45" << std::endl;
    std::cerr << "here converged.NS.9.info is a 3D single NS solution and"
              << std::endl;
    std::cerr << "the separation is an optional parameter." << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  if (argc > 2) {
    distkm = std::atof(argv[2]);
    if (rank == 0)
      std::cout << "Separation Distance set to: " << distkm << "km"
                << std::endl;
  } else if (rank == 0)
    std::cout << "Separation Distance set to default " << distkm << "km"
              << std::endl;

  // load single NS configuration
  std::string nsfilename{argv[1]};
  kadath_config_boost<BCO_NS_INFO> NSconfig(nsfilename);

  // print imported star parameters
  std::cout << NSconfig << std::endl;

  // call the binary setup routine
  if (NSconfig(DIM) == 3)
    setup_3d(NSconfig);
  else
    std::cerr << "ns.dim not found. Either add manually or rerun NS solver\n "
                 "Only 3D stars are supported (yet).";

  MPI_Finalize();
  return EXIT_SUCCESS;
}

/**
 * setup_bin_config
 * Setup and initialization of shared binary parameters independent of 2d or 3d star import
 *
 * @tparam config_t config file type
 * @tparam space_t space file type
 * @param[input] NSconfig neutron star configuration file
 * @param[input] spacein1 space from the isolated NS
 */
template <typename config_t, typename space_t>
kadath_config_boost<BIN_INFO> setup_bin_config(config_t& NSconfig,
                                               const space_t& spacein1) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Initialize Config File */
  std::string ifilename = "./initbin.info";
  kadath_config_boost<BIN_INFO> bconfig;
  bconfig.set_filename(ifilename);
  bconfig.initialize_binary({"ns", "ns"});

  /* copy BCO Parameters from NS config file*/
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO1) = NSconfig.set(i);
  }

  bconfig.set(NSHELLS, BCO1) = 0;

  /* copy EOS Parameters */
  for (int i = 0; i < NUM_EOS_PARAMS; ++i) {
    bconfig.set_eos(i, BCO1) = NSconfig.set_eos(i);
  }

  // setup eos and update central density
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

  // setup boundaries of the stellar domains
  auto [r_min, r_max] = bco_utils::get_rmin_rmax(spacein1, 1);
  bconfig.set(RIN, BCO1) = 0.5 * r_min;
  bconfig.set(RMID, BCO1) = r_max;
  bconfig.set(ROUT, BCO1) = 1.5 * r_max;

  // copy parameters so both stars are the same - EOS is always read from BCO1
  for (int i = 0; i < NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, BCO2) = bconfig.set(i, BCO1);
  }

  /* Binary Parameters */
  bconfig.set(BIN_RES) = 9;  ///> Always solve at low res first
  bconfig.set(OUTER_SHELLS) = 0;
  bconfig.set(DIST) = distkm / M2km;
  // radius of the exterior domain, matching to the compactified domain
  bconfig.set(REXT) = 2 * bconfig(DIST);
  bconfig.set(CHECKPOINT) = 0;
  bconfig.set(QPIG) = NSconfig(BCO_QPIG);
  // equal mass system with "center of mass" at the origin
  bconfig.set(Q) = 1;
  bconfig.set(COM) = 0.;
  bconfig.set(OUTER_SHELLS) = 0;

  // obtain 3PN estimate for the global, orbital omega
  bco_utils::KadathPNOrbitalParams(bconfig, bconfig(MADM, BCO1),
                                   bconfig(MADM, BCO2));

  // delete ADOT, this can always be recalculated during
  // the eccentricity reduction stage
  bconfig.set(ADOT) = std::nan("1");

  // catalog fields that are implemented
  bconfig.set_field(CONF) = true;
  bconfig.set_field(LAPSE) = true;
  bconfig.set_field(SHIFT) = true;
  bconfig.set_field(LOGH) = true;
  bconfig.set_field(PHI) = true;

  // set the solver stages
  // only total is set, as that is needed first
  // before attempting more interesting setups q!=1, chi !=0
  bconfig.set_stage(TOTAL) = true;

  // output configuration
  if (rank == 0)
    std::cout << "\n" << bconfig;

  return bconfig;
}

/* Binary setup from a single 3D star import */
template <typename config_t>
void setup_3d(config_t NSconfig) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // read single NS configuration and data
  std::string nsspacein = NSconfig.space_filename();

  FILE* ff1 = fopen(nsspacein.c_str(), "r");
  Space_spheric_adapted spacein1(ff1);
  Scalar confin1(spacein1, ff1);
  Scalar lapsein1(spacein1, ff1);
  Vector shiftin1(spacein1, ff1);
  Scalar loghin1(spacein1, ff1);
  fclose(ff1);
  int ndomin1 = spacein1.get_nbr_domains();

  // overwrite HC before passing NSconfig to setup_bin_config
  NSconfig.set(HC) =
      std::exp(bco_utils::get_boundary_val(0, loghin1, INNER_BC));

  auto bconfig = setup_bin_config(NSconfig, spacein1);
  //end updating config vars

  // use Chebychev
  int typer = CHEB_TYPE;
  // get the adapted domain of the single star
  const Domain_shell_outer_adapted* old_outer_adapted1 =
      dynamic_cast<const Domain_shell_outer_adapted*>(spacein1.get_domain(1));

  // create scalar field representing the radius in the single star space
  Scalar old_space_radius(spacein1);
  old_space_radius.annule_hard();

  old_space_radius.set_domain(0) = spacein1.get_domain(0)->get_radius();
  old_space_radius.set_domain(1) = old_outer_adapted1->get_outer_radius();
  for (int i = 2; i < ndomin1; ++i) {
    old_space_radius.set_domain(i) = spacein1.get_domain(i)->get_radius();
  }
  old_space_radius.std_base();
  // end create old radius scalar fields

  // output the equitorial and polar radii
  auto [rmin, rmax] = bco_utils::get_rmin_rmax(spacein1, 1);

  std::cout << "Rmin/max: " << std::endl << rmin << " " << rmax << std::endl;

  // setup domain boundaries
  std::vector<double> out_bounds(1 + bconfig(OUTER_SHELLS));
  std::vector<double> NS1_bounds(3 + bconfig(NSHELLS, BCO1));
  std::vector<double> NS2_bounds(3 + bconfig(NSHELLS, BCO2));

  // scale outer shells by constant steps of 1/4 for the time being
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(REXT) * (1. + e * 0.25);

  // set reasonable radii to each stellar domain
  bco_utils::set_NS_bounds(NS1_bounds, bconfig, BCO1);
  bco_utils::set_NS_bounds(NS2_bounds, bconfig, BCO2);
  // end setup domain boundaries

  // create a binary neutron star space
  Space_bin_ns space(typer, bconfig(DIST), NS1_bounds, NS2_bounds, out_bounds,
                     bconfig(BIN_RES));
  // with cartesian type basis
  Base_tensor basis(space, CARTESIAN_BASIS);

  // output stellar domain radii
  if (rank == 0) {
    std::cout << "Bounds:" << std::endl;

    bco_utils::print_bounds("NS1", NS1_bounds);
    bco_utils::print_bounds("NS2", NS2_bounds);
  }

  // get the domains with inner respectively outer adapted boundary
  // for each of the stars
  std::array<const Domain_shell_outer_adapted*, 2> new_outer_adapted{
      dynamic_cast<const Domain_shell_outer_adapted*>(
          space.get_domain(space.ADAPTED1)),
      dynamic_cast<const Domain_shell_outer_adapted*>(
          space.get_domain(space.ADAPTED2))};

  std::array<const Domain_shell_inner_adapted*, 2> new_inner_adapted{
      dynamic_cast<const Domain_shell_inner_adapted*>(
          space.get_domain(space.ADAPTED1 + 1)),
      dynamic_cast<const Domain_shell_inner_adapted*>(
          space.get_domain(space.ADAPTED2 + 1))};

  // updated the radius mapping for both NS
  for (int i = 0; i < 2; ++i) {
    bco_utils::interp_adapted_mapping(new_inner_adapted[i], 1,
                                      old_space_radius);
    bco_utils::interp_adapted_mapping(new_outer_adapted[i], 1,
                                      old_space_radius);
  }

  // get and print center of each star
  double xc1 = bco_utils::get_center(space, space.NS1);
  double xc2 = bco_utils::get_center(space, space.NS2);

  if (rank == 0)
    std::cout << "xc1: " << xc1 << std::endl << "xc2: " << xc2 << std::endl;

  // start to create the new fields
  // initialized to zero globally
  Scalar logh(space);
  logh.annule_hard();

  Scalar conf(space);
  conf.annule_hard();

  Scalar lapse(space);
  lapse.annule_hard();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();

  Scalar phi(space);
  phi.annule_hard();
  phi.std_base();
  //end create new fields

  // start importing the fields from the single star
  int ndom = space.get_nbr_domains();
  for (int dom = 0; dom < ndom; dom++) {
    // get an index in each domain to iterate over all colocation points
    Index new_pos(space.get_domain(dom)->get_nbr_points());

    do {
      // get cartesian coordinates of the current colocation point
      double x = space.get_domain(dom)->get_cart(1)(new_pos);
      double y = space.get_domain(dom)->get_cart(2)(new_pos);
      double z = space.get_domain(dom)->get_cart(3)(new_pos);

      // define a point shifted suitably to the stellar centers in the binary
      Point absol1(3);
      absol1.set(1) = x - xc1;
      absol1.set(2) = y;
      absol1.set(3) = z;

      Point absol2(3);
      absol2.set(1) = x - xc2;
      absol2.set(2) = y;
      absol2.set(3) = z;

      // initialize lapse and conformal factor as product of both shifted single star solutions
      // this gives an acceptable initial guess
      if (dom != ndom - 1) {
        conf.set_domain(dom).set(new_pos) =
            confin1.val_point(absol1) * confin1.val_point(absol2);
        lapse.set_domain(dom).set(new_pos) =
            lapsein1.val_point(absol1) * lapsein1.val_point(absol2);
        logh.set_domain(dom).set(new_pos) =
            loghin1.val_point(absol1) + loghin1.val_point(absol2);
      } else {
        // We have to set the compactified domain manually since the outer collocation point is always
        // at inf which is undefined numerically
        conf.set_domain(dom).set(new_pos) = 1.;
        lapse.set_domain(dom).set(new_pos) = 1.;
        logh.set_domain(dom).annule_hard();
      }
      // loop over all colocation points
    } while (new_pos.inc());
  }  //end importing fields

  // set logh to zero outside of stars explicitly
  logh.set_domain(space.ADAPTED1 + 1).annule_hard();
  for (int i = space.ADAPTED2 + 1; i < ndom; ++i)
    logh.set_domain(i).annule_hard();

  // employ standard spectral expansion, compatible with the given paraties
  conf.std_base();
  lapse.std_base();
  logh.std_base();
  shift.std_base();

  // save everything to a binary file
  if (rank == 0)
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
}
