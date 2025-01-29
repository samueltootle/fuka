/*
 * Copyright 2021
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author:
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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

#include <cmath>
#include <iterator>
#include <sstream>
#include "coord_fields.hpp"
#include "kadath_bin_bh.hpp"
#include "mpi.h"
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"

using namespace Kadath;
using namespace Kadath::FUKA_Config;

// forward declarations
template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics(space_t const& space,
                       syst_t const& syst,
                       int const ite,
                       double const conv,
                       config_t bconfig);

template <typename space_t, typename syst_t>
void print_diagnostics_pre(space_t const& space,
                           syst_t const& syst,
                           int const ite,
                           double const conv);

template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig);

// end forward declarations

void print_stage(std::string stage_name) {
  std::cout << "############################" << std::endl
            << stage_name << std::endl
            << "############################" << std::endl;
}

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  std::string input_filename;
  if (argc < 2) {
    std::cerr << "Usage: ./solve /<path>/<ID base name>.info "
                 "./<output_path>/"
              << endl;
    std::cerr << "e.g. ./solve initbin.info ./out" << endl;
    std::cerr << "Note: <output_path> is optional and defaults to <path>"
              << endl;
    std::_Exit(EXIT_FAILURE);
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // reading from configuration file
  input_filename = argv[1];
  std::string file_path = extract_path(input_filename);
  kadath_config_boost<BIN_INFO> bconfig(input_filename);

  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();

  // file containing KADATH fields must have same name as config file with only the extension being different
  std::string kadath_filename = bconfig.space_filename();
  if (!fs::exists(kadath_filename)) {
    if (rank == 0) {
      std::cerr << "File: " << kadath_filename << " not found.\n\n";
      std::_Exit(EXIT_FAILURE);
    }
  }

  FILE* fin = fopen(kadath_filename.c_str(), "r");
  Space_bin_bh space(fin);
  Base_tensor basis(space, CARTESIAN_BASIS);
  Scalar conf(space, fin);
  Scalar lapse(space, fin);
  Vector shift(space, fin);
  fclose(fin);

  // set shift to zero - mainly for testing
  if (bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  // update output path based on input arg
  if (argc > 2)
    file_path = argv[2];
  bconfig.set_outputdir(file_path);

  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);
  if (rank == 0) {
    std::cout << "Last Stage Enabled: " << last_stage << std::endl;
  }

  double rm = bco_utils::get_radius(space.get_domain(space.BH1 + 2), INNER_BC);
  double rp = bco_utils::get_radius(space.get_domain(space.BH2 + 2), INNER_BC);

  Metric_flat fmet(space, basis);

  int ndom = space.get_nbr_domains();

  // ignore the interior, excised domains
  std::array<int, 4> excluded_doms{space.BH1, space.BH1 + 1, space.BH2,
                                   space.BH2 + 1};

  // setup coordinate fields
  double xc1 = bco_utils::get_center(space, space.BH1);
  double xc2 = bco_utils::get_center(space, space.BH2);

  double xo = bco_utils::get_center(space, ndom - 1);
  CoordFields<Space_bin_bh> cf_generator(space);

  vec_ary_t coord_vectors = default_binary_vector_ary(space);

  update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2);
  // end setup coord fields

  if (rank == 0) {
    std::cout << "=================================" << endl;
    std::cout << bconfig;
    std::cout << "Mass ratio: " << bconfig(Q) << endl;
    std::cout << "Radii: " << rm << " " << rp << endl;
    std::cout << "Pos: " << xc1 << " " << xc2 << endl;
    std::cout << "=================================" << endl;
  }

  // We initially solve without the shift to obtain a better, smooth guess for conf and lapse
  // before slowly introducing the shift.  This is necessary since we are starting
  // the binary from scratch.
  if (stage_enabled[PRE]) {

    // setup radius scalar fields centered on each BH to be used to fix the radii
    // of the excision regions manually
    scalar_ary_t coord_scalars{Scalar(space), Scalar(space)};
    update_fields(cf_generator, coord_vectors, coord_scalars, xo, xc1, xc2);

    // radius field centered on BH1 such that level_minus has a root at the intended radius
    Scalar level_minus(space);
    level_minus = *coord_scalars[R_BCO1] * *coord_scalars[R_BCO1] - rm * rm;
    level_minus.std_base();

    // radius field centered on BH2 such that level_plus has a root at the intended radius
    Scalar level_plus(space);
    level_plus = *coord_scalars[R_BCO2] * *coord_scalars[R_BCO2] - rp * rp;
    level_plus.std_base();

    if (rank == 0)
      print_stage("Laplace lapse and psi preconditioning");

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    // equation constants
    syst.add_cst("n0", bconfig(FIXED_LAPSE, BCO1));

    // constant vector fields
    syst.add_cst("rm", *coord_scalars[R_BCO1]);
    syst.add_cst("rp", *coord_scalars[R_BCO2]);

    syst.add_cst("levm", level_minus);
    syst.add_cst("levp", level_plus);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);

    // variable fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);

    // definitions for equations
    syst.add_def("NP = P*N");

    syst.add_def("eqP = D^i D_i P");
    syst.add_def("eqNP = D^i D_i NP");

    // lapse and conformal factor equations
    space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq(syst, "eqP= 0", "P", "dn(P)");

    // boundary conditions at inf
    space.add_bc_outer(syst, "NP = 1");
    space.add_bc_outer(syst, "P = 1");

    // boundary conditions on BH1
    space.add_bc_sphere_one(syst, "N = n0");
    // FIXME is this comment correct?
    // this assumes a conformally flat background metric
    space.add_bc_sphere_one(syst, "sm^j * D_j P + P / 4 * D^j sm_j = 0");

    // boundary conditions on BH2
    space.add_bc_sphere_two(syst, "N = n0");
    space.add_bc_sphere_two(syst, "sp^j * D_j P + P / 4 * D^j sp_j = 0");

    // we fix the BH area using our "level" fields such that the horizon
    // is always at the roots of lev* (i.e. the input RMID)
    space.add_eq_int_sphere_one(syst, "integ(levm) = 0");
    space.add_eq_int_sphere_two(syst, "integ(levp) = 0");

    // exclude the domains representing the excised regions
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
    }

    bool endloop = false;
    int ite = 1;
    double conv;

    // begin interative solver
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "grav_np_pre_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {

        print_diagnostics_pre(space, syst, ite, conv);

        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
      }
      ite++;
    }  // end iterative solver

    bconfig.set_filename(converged_filename("PRE", bconfig));
    if (rank == 0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
    if (last_stage_idx != PRE)
      bconfig.set_stage(PRE) = false;

    // delete fields since they are not used by other stages.
    for (auto& el : coord_scalars) {
      el.reset();
    }
  }  //end PRE stage

  // For the forced equal mass stages, we restrict ourselves to the irreducible mass of BCO1
  // this allows one to put different masses in the config file from the start without impacting
  // the early stages
  auto equal_mass_stages = [&](const auto stage, std::string stage_text) {
    if (rank == 0 && stage == FIXED_OMEGA)
      print_stage(
          "Fixed orbital frequency and fixed equal masses in corotation");
    if (rank == 0 && stage == COROT_EQUAL)
      print_stage("Corotation with fixed equal masses");

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    // system constants
    syst.add_cst("PI", M_PI);

    // BH fixing parameters
    syst.add_cst("n0", bconfig(FIXED_LAPSE, BCO1));
    syst.add_cst("M", bconfig(MIRR, BCO1));

    // determine whether the orbital frequency is fixed
    if (bconfig.control(FIXED_GOMEGA))
      syst.add_cst("ome", bconfig(GOMEGA));
    else
      syst.add_var("ome", bconfig(GOMEGA));

    // constant coordinate vector fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // FIXME fixed com shift? is this either zero or from a previous solution?
    //       shouldn't this be zero in the equal mass cases?
    // fixed shifting of the coordinate "center of mass"
    syst.add_cst("xaxis", bconfig(COM));

    // variable fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    // begin of the definitions needed for the system of equations
    syst.add_def("NP      = P*N");
    syst.add_def("Ntilde  = N / P^6");

    // definition describing the orbital rotation field
    syst.add_def("Morb^i  = mg^i + xaxis * ey^i");
    // total shift = inertial + rotation contributions
    syst.add_def("B^i     = bet^i + ome * Morb^i");
    // conformal extrinsic curvature
    syst.add_def(
        "A^ij    = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    // ADM linear momenta
    syst.add_def("intPx   = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy   = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz   = A_ij * ez^j * einf^i / 8 / PI");

    // quasi-local spin definition on the apparent horizons
    syst.add_def("intSm   = A_ij * mm^i * sm^j   / 8 / PI");
    syst.add_def("intSp   = A_ij * mp^i * sp^j   / 8 / PI");

    // BH irreducible mass, given by the proper area
    syst.add_def("intMsq  = P^4 / 16. / PI");

    // vacuum XCTS constraint equations
    syst.add_def("eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");
    // end definitions

    // system of equations - see the space for details regarding add_eq
    space.add_eq(syst, "eqNP    = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP     = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions on variable fields at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i = 0");

    // BCs on the apparent horizon of BH1
    space.add_bc_sphere_one(syst, "N = n0");
    space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i");
    space.add_bc_sphere_one(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_one(syst, "integ(intMsq) - M * M = 0 ");

    // BCs on the apparent horizon of BH2
    space.add_bc_sphere_two(syst, "N = n0");
    space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");
    space.add_bc_sphere_two(
        syst,
        "sp^j * D_j P + P / 4 * D^j sp_j + A_ij * sp^i * sp^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_two(syst, "integ(intMsq) - M * M = 0 ");

    // quasi-equallibrium condition at infinity, if binary orbital frequency is not fixed.
    if (!bconfig.control(FIXED_GOMEGA))
      space.add_eq_int_inf(syst, "integ(dn(N) + 2 * dn(P)) = 0");

    // optimally these shouldn't be part of the system
    // exclusion of fields / domains not implemented properly yet
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2, &syst);

    bool endloop = false;
    int ite = 1;
    double conv;

    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bbh_" + stage_text + "_" << ite - 1;
      bconfig.set_filename(ss.str());

      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);

        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
      }

      update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2, &syst);
      ite++;
    }
    std::string stage_name = "COROT";
    if (stage == FIXED_OMEGA)
      stage_name = "FIXED_" + stage_name;

    bconfig.set_filename(converged_filename(stage_name, bconfig));
    if (rank == 0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift);

    if (last_stage_idx != stage)
      bconfig.set_stage(stage) = false;
  };

  if (stage_enabled[FIXED_OMEGA]) {
    // in case fixed orbital frequency isn't already set, set it temporarily.
    if (!bconfig.control(FIXED_GOMEGA)) {
      bconfig.control(FIXED_GOMEGA) = true;
      equal_mass_stages(FIXED_OMEGA, "fixed_corot_equal");
      bconfig.control(FIXED_GOMEGA) = false;
    } else
      equal_mass_stages(FIXED_OMEGA, "fixed_corot_equal");
  }  // end fixed_omega stage

  if (stage_enabled[COROT_EQUAL]) {
    equal_mass_stages(COROT_EQUAL, "corot_equal");
  }  // end corot_equal stage

  if (stage_enabled[TOTAL]) {
    if (rank == 0 && bconfig.control(COROT_BIN))
      print_stage("Total system - COROT");
    else if (rank == 0)
      print_stage("Total system");

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);

    // BH1 fixing parameters
    syst.add_cst("n0m", bconfig(FIXED_LAPSE, BCO1));
    syst.add_cst("Mm", bconfig(MIRR, BCO1));
    syst.add_cst("chim", bconfig(CHI, BCO1));
    syst.add_cst("CMm", bconfig(MCH, BCO1));

    // BH2 fixing parameters
    syst.add_cst("n0p", bconfig(FIXED_LAPSE, BCO2));
    syst.add_cst("Mp", bconfig(MIRR, BCO2));
    syst.add_cst("chip", bconfig(CHI, BCO2));
    syst.add_cst("CMp", bconfig(MCH, BCO2));

    // constant vector fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // variable fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    // determine whether orbital frequency is fixed
    if (bconfig.control(FIXED_GOMEGA))
      syst.add_cst("ome", bconfig(GOMEGA));
    else
      syst.add_var("ome", bconfig(GOMEGA));

    // setup system accordingly for a corotating binary
    if (bconfig.control(COROT_BIN)) {
      bconfig.set(OMEGA, BCO1) = 0.;
      bconfig.set(OMEGA, BCO2) = 0.;
    } else {
      // otherwise we setup the system for arbitrarily rotating BHs
      // based on either fixed rotation frequency or variable
      // rotation frequency fixed by a dimensionless
      // spin parameter CHI.
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst("omesm", bconfig(FIXED_BCOMEGA, BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
      } else {
        syst.add_var("omesm", bconfig(OMEGA, BCO1));
      }
      syst.add_def(space.BH1 + 2, "ExOme^i = omesm * mm^i");

      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst("omesp", bconfig(FIXED_BCOMEGA, BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
      } else {
        syst.add_var("omesp", bconfig(OMEGA, BCO2));
      }
      syst.add_def(space.BH2 + 2, "ExOme^i = omesp * mp^i");
    }

    // allowed for a shift on the coordinate x axis to allow the system to find
    // the "center of mass" by the vaninishing of the ADM linear momentum at infinity.
    syst.add_var("xaxis", bconfig(COM));

    // begin definitions needed for system of equations
    syst.add_def("NP      = P*N");
    syst.add_def("Ntilde  = N / P^6");

    // definition describing the binary rotation field
    syst.add_def("Morb^i  = mg^i + xaxis * ey^i");
    // total shift = inertial + rotation contributions
    syst.add_def("B^i = bet^i + ome * Morb^i");

    // conformal extrinsic curvature
    syst.add_def(
        "A^ij    = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    // ADM linear momenta
    syst.add_def("intPx   = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy   = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz   = A_ij * ez^j * einf^i / 8 / PI");

    // spin definition on the apparent horizons
    syst.add_def("intSm   = A_ij * mm^i * sm^j   / 8 / PI");
    syst.add_def("intSp   = A_ij * mp^i * sp^j   / 8 / PI");

    // BH irreducible mass given by the proper area
    syst.add_def("intMsq  = P^4 / 16. / PI");

    // vacuum XCTS constraint equations
    syst.add_def("eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");
    // end definitions

    // system of equations - see space for details regarding add_eq
    space.add_eq(syst, "eqNP    = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP     = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions on variable fields at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i = 0");

    // quasi-equallibrium condition if binary orbital frequency is not fixed.
    if (!bconfig.control(FIXED_GOMEGA))
      space.add_eq_int_inf(syst, "integ(dn(N) + 2 * dn(P)) = 0");

    space.add_eq_int_inf(syst, "integ(intPy) = 0");

    // BCs on the apparent horizon of BH1
    space.add_bc_sphere_one(syst, "N = n0m");
    space.add_bc_sphere_one(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_one(syst, "integ(intMsq) - Mm * Mm = 0 ");

    // BCs on the apparent horizon of BH2
    space.add_bc_sphere_two(syst, "N = n0p");
    space.add_bc_sphere_two(
        syst,
        "sp^j * D_j P + P / 4 * D^j sp_j + A_ij * sp^i * sp^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_two(syst, "integ(intMsq) - Mp * Mp = 0 ");

    // in the case of corotation, the tangential term is zero by definition
    if (bconfig.control(COROT_BIN)) {
      space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i");
      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");

    } else {
      // otherwise we include equations on the boundary with the addition of the tangential spin component
      // note that these assume a conformally flat background metric
      space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i + ExOme^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        space.add_eq_int_sphere_one(syst,
                                    "integ(intSm) - chim * CMm * CMm = 0 ");
      }
      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + ExOme^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        space.add_eq_int_sphere_two(syst,
                                    "integ(intSp) - chip * CMp * CMp = 0 ");
      }
    }

    // exclude the excised regions
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    if (rank == 0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    // begin iterative solver
    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bbh_spin_mass_fixed_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);
        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(ss, space, bconfig, conf, lapse, shift);
      }

      update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2, &syst);
      ite++;
    }

    // in the case of corotation or fixed local spin frequency, we calculate MCH accordingly
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)) ||
        bconfig.control(COROT_BIN))
      bconfig(MCH, BCO1) =
          bco_utils::syst_mch(syst, space, "intSm", space.BH1 + 2);
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)) ||
        bconfig.control(COROT_BIN))
      bconfig(MCH, BCO2) =
          bco_utils::syst_mch(syst, space, "intSp", space.BH2 + 2);

    bconfig.set_filename(converged_filename("TOTAL", bconfig));
    if (rank == 0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift);

    if (last_stage_idx != TOTAL)
      bconfig.set_stage(TOTAL) = false;
  }  // TOTAL stage

  auto von_Neumann_stages = [&](const auto stage, std::string stage_text) {
    if (stage == ECC_RED) {
      // determine whether to use PN estimates of the orbital frequency and adot
      // in case of the eccentricity stage
      if (std::isnan(bconfig.set(ADOT)) || std::isnan(bconfig.set(ECC_OMEGA)) ||
          bconfig.control(USE_PN)) {

        bco_utils::KadathPNOrbitalParams(bconfig, bconfig(MCH, BCO1),
                                         bconfig(MCH, BCO2));

        if (rank == 0)
          std::cout << "### Using PN estimate for adot and omega! ###"
                    << std::endl;
      }
      bconfig.set(GOMEGA) = bconfig(ECC_OMEGA);
    }
    // setup radial position vector field for adot - only needed for ECC_RED stage
    Vector CART(space, CON, basis);
    CART = cf_generator.cart();

    if (rank == 0 && stage == ECC_RED)
      print_stage("Ecc. Reduction - Fixed Omega, adot, spin, and mass");
    if (rank == 0 && bconfig.control(COROT_BIN))
      print_stage("Total system with v.Neumann BC - COROT");
    else if (rank == 0)
      print_stage("Total system with v.Neumann BC");

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);

    // BH1 fixing parameters
    syst.add_cst("Mm", bconfig(MIRR, BCO1));
    syst.add_cst("chim", bconfig(CHI, BCO1));
    syst.add_cst("CMm", bconfig(MCH, BCO1));

    // BH2 fixing parameters
    syst.add_cst("Mp", bconfig(MIRR, BCO2));
    syst.add_cst("chip", bconfig(CHI, BCO2));
    syst.add_cst("CMp", bconfig(MCH, BCO2));

    // constant vector fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // variable fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    // determine whether orbital frequency is fixed
    if (bconfig.control(FIXED_GOMEGA) || stage == ECC_RED) {
      syst.add_cst("ome", bconfig(GOMEGA));
    } else {
      syst.add_var("ome", bconfig(GOMEGA));
    }

    // allow for a shift of the "center of mass" on the x-axis
    // enforced by a vanishing ADM linear momentum
    syst.add_var("xaxis", bconfig(COM));
    syst.add_var("yaxis", bconfig(COMY));

    // setup parameters accordingly for a corotating binary
    if (bconfig.control(COROT_BIN)) {
      bconfig.set(OMEGA, BCO1) = 0.;
      bconfig.set(OMEGA, BCO2) = 0.;
    } else {
      // otherwise we setup the system for arbitrarily rotating BHs
      // based on either fixed rotation frequency or variable
      // rotation frequency fixed by a dimensionless
      // spin parameter CHI.
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst("omesm", bconfig(FIXED_BCOMEGA, BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
      } else
        syst.add_var("omesm", bconfig(OMEGA, BCO1));

      syst.add_def(space.BH1 + 2, "ExOme^i = omesm * mm^i");

      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst("omesp", bconfig(FIXED_BCOMEGA, BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
      } else
        syst.add_var("omesp", bconfig(OMEGA, BCO2));
      syst.add_def(space.BH2 + 2, "ExOme^i = omesp * mp^i");
    }

    // begin definitions needed for system of equations
    syst.add_def("NP      = P*N");
    syst.add_def("Ntilde  = N / P^6");

    // definition describing the binary rotation
    syst.add_def("Morb^i  = mg^i + xaxis * ey^i + yaxis * ex^i");

    std::string bigB{"B^i = bet^i + ome * Morb^i"};
    if (stage == ECC_RED) {
      syst.add_cst("adot", bconfig(ADOT));
      syst.add_cst("r", CART);
      syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
      // add contribution of ADOT to the total shift definition
      bigB += " + adot * comr^i";
    }

    // total shift = inertial + orbital contributions
    syst.add_def(bigB.c_str());
    // conformal extrinsic curvature
    syst.add_def(
        "A^ij    = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    // ADM linear momenta
    syst.add_def("intPx   = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy   = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz   = A_ij * ez^j * einf^i / 8 / PI");

    // spin definition on the apparent horizons
    syst.add_def("intSm   = A_ij * mm^i * sm^j   / 8 / PI");
    syst.add_def("intSp   = A_ij * mp^i * sp^j   / 8 / PI");

    // BH irreducible mass
    syst.add_def("intMsq  = P^4 / 16. / PI");

    // vacuum XCTS constraint equations
    syst.add_def("eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");
    // end definitions

    // system of equations - see space for details regarding add_eq
    space.add_eq(syst, "eqNP    = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP     = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions on variable fields at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i = 0");

    // quasi-equallibrium condition if binary orbital frequency is not fixed
    if (!bconfig.control(FIXED_GOMEGA) && stage != ECC_RED) {
      space.add_eq_int_inf(syst, "integ(dn(N) + 2 * dn(P)) = 0");
    }

    // minimize ADM linear momenta at infinity, Pz is zero by symmetry
    space.add_eq_int_inf(syst, "integ(intPx) = 0");
    space.add_eq_int_inf(syst, "integ(intPy) = 0");

    // BCs on the apparent horizon of BH1
    space.add_bc_sphere_one(syst, "dn(NP) = 0");
    space.add_bc_sphere_one(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_one(syst, "integ(intMsq) - Mm * Mm = 0 ");

    // BCs on the apparent horizon of BH2
    space.add_bc_sphere_two(syst, "dn(NP)=0");
    space.add_bc_sphere_two(
        syst,
        "sp^j * D_j P + P / 4 * D^j sp_j + A_ij * sp^i * sp^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_two(syst, "integ(intMsq) - Mp * Mp = 0 ");

    // in the case of corotation, the tangential term is zero by definition
    if (bconfig.control(COROT_BIN)) {
      space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i");
      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");

    } else {
      // otherwise we include equations on the boundary with the addition of the tangential spin component
      // note that these assume a conformally flat background metric
      space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i + ExOme^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_sphere_one(syst,
                                    "integ(intSm) - chim * CMm * CMm = 0 ");

      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + ExOme^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_sphere_two(syst,
                                    "integ(intSp) - chip * CMp * CMp = 0 ");
    }

    // exclude the excised regions
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    if (rank == 0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    // begin iterative solver
    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bbh_" + stage_text << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);

        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(ss, space, bconfig, conf, lapse, shift);
      }
      update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2, &syst);
      ite++;
    }

    // in the case of corotation or fixed local spin frequency, we calculate MCH accordingly
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)) ||
        bconfig.control(COROT_BIN))
      bconfig(MCH, BCO1) =
          bco_utils::syst_mch(syst, space, "intSm", space.BH1 + 2);
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)) ||
        bconfig.control(COROT_BIN))
      bconfig(MCH, BCO2) =
          bco_utils::syst_mch(syst, space, "intSp", space.BH2 + 2);

    bconfig.set_filename(converged_filename(stage_text, bconfig));
    if (rank == 0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  };

  if (stage_enabled[TOTAL_BC]) {
    von_Neumann_stages(TOTAL_BC, "TOTAL_BC");
  }  // TOTAL_BC Stage

  if (stage_enabled[ECC_RED]) {
    von_Neumann_stages(ECC_RED, "ECC_RED");
  }

  auto test_stage = [&](const auto stage, std::string stage_text) {
    // setup radial position vector field for adot - only needed for ECC_RED stage
    Vector CART(space, CON, basis);
    CART = cf_generator.cart();

    if (rank == 0 && bconfig.control(COROT_BIN))
      print_stage("Testing with v.Neumann BC - COROT");
    else if (rank == 0)
      print_stage("Testing with v.Neumann BC");

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);

    // BH1 fixing parameters
    syst.add_cst("Mm", bconfig(MIRR, BCO1));
    syst.add_cst("chim", bconfig(CHI, BCO1));
    syst.add_cst("CMm", bconfig(MCH, BCO1));

    // BH2 fixing parameters
    syst.add_cst("Mp", bconfig(MIRR, BCO2));
    syst.add_cst("chip", bconfig(CHI, BCO2));
    syst.add_cst("CMp", bconfig(MCH, BCO2));

    // constant vector fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // variable fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    // determine whether orbital frequency is fixed
    // FIXME COM should likely not be fixed.  This is
    // here for testing superimposed import of highly asymmetric
    // binaries
    if (bconfig.control(FIXED_GOMEGA) || stage == ECC_RED) {
      syst.add_cst("ome", bconfig(GOMEGA));
      syst.add_cst("xaxis", bconfig(COM));
      syst.add_cst("yaxis", bconfig(COMY));
    } else {
      syst.add_var("ome", bconfig(GOMEGA));
      // allow for a shift of the "center of mass" on the x-axis
      // enforced by a vanishing ADM linear momentum
      syst.add_var("xaxis", bconfig(COM));
      syst.add_var("yaxis", bconfig(COMY));
    }

    // setup parameters accordingly for a corotating binary
    if (bconfig.control(COROT_BIN)) {
      bconfig.set(OMEGA, BCO1) = 0.;
      bconfig.set(OMEGA, BCO2) = 0.;
    } else {
      // otherwise we setup the system for arbitrarily rotating BHs
      // based on either fixed rotation frequency or variable
      // rotation frequency fixed by a dimensionless
      // spin parameter CHI.
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst("omesm", bconfig(FIXED_BCOMEGA, BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
      } else
        syst.add_var("omesm", bconfig(OMEGA, BCO1));

      syst.add_def(space.BH1 + 2, "ExOme^i = omesm * mm^i");

      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst("omesp", bconfig(FIXED_BCOMEGA, BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
      } else
        syst.add_var("omesp", bconfig(OMEGA, BCO2));
      syst.add_def(space.BH2 + 2, "ExOme^i = omesp * mp^i");
    }

    // begin definitions needed for system of equations
    syst.add_def("NP      = P*N");
    syst.add_def("Ntilde  = N / P^6");

    // definition describing the binary rotation
    syst.add_def("Morb^i  = mg^i + xaxis * ey^i + yaxis * ex^i");

    std::string bigB{"B^i = bet^i + ome * Morb^i"};
    if (stage == ECC_RED) {
      syst.add_cst("adot", bconfig(ADOT));
      syst.add_cst("r", CART);
      syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
      // add contribution of ADOT to the total shift definition
      bigB += " + adot * comr^i";
    }

    // total shift = inertial + orbital contributions
    syst.add_def(bigB.c_str());
    // conformal extrinsic curvature
    syst.add_def(
        "A^ij    = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    // ADM linear momenta
    syst.add_def("intPx   = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy   = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz   = A_ij * ez^j * einf^i / 8 / PI");

    // spin definition on the apparent horizons
    syst.add_def("intSm   = A_ij * mm^i * sm^j   / 8 / PI");
    syst.add_def("intSp   = A_ij * mp^i * sp^j   / 8 / PI");

    // BH irreducible mass
    syst.add_def("intMsq  = P^4 / 16. / PI");

    // vacuum XCTS constraint equations
    syst.add_def("eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");
    // end definitions

    // system of equations - see space for details regarding add_eq
    space.add_eq(syst, "eqNP    = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP     = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions on variable fields at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P     = 1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i = 0");

    // quasi-equallibrium condition if binary orbital frequency is not fixed
    if (!bconfig.control(FIXED_GOMEGA) && stage != ECC_RED) {
      // minimize ADM linear momenta at infinity, Pz is zero by symmetry
      space.add_eq_int_inf(syst, "integ(intPx) = 0");
      space.add_eq_int_inf(syst, "integ(intPy) = 0");
      space.add_eq_int_inf(syst, "integ(dn(N) + 2 * dn(P)) = 0");
    }

    // BCs on the apparent horizon of BH1
    space.add_bc_sphere_one(syst, "dn(NP) = 0");
    space.add_bc_sphere_one(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_one(syst, "integ(intMsq) - Mm * Mm = 0 ");

    // BCs on the apparent horizon of BH2
    space.add_bc_sphere_two(syst, "dn(NP)=0");
    space.add_bc_sphere_two(
        syst,
        "sp^j * D_j P + P / 4 * D^j sp_j + A_ij * sp^i * sp^j / P^3 / 4 = 0");
    space.add_eq_int_sphere_two(syst, "integ(intMsq) - Mp * Mp = 0 ");

    // in the case of corotation, the tangential term is zero by definition
    if (bconfig.control(COROT_BIN)) {
      space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i");
      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");

    } else {
      // otherwise we include equations on the boundary with the addition of the tangential spin component
      // note that these assume a conformally flat background metric
      space.add_bc_sphere_one(syst, "B^i = N / P^2 * sm^i + ExOme^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_sphere_one(syst,
                                    "integ(intSm) - chim * CMm * CMm = 0 ");

      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + ExOme^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_sphere_two(syst,
                                    "integ(intSp) - chip * CMp * CMp = 0 ");
    }

    // exclude the excised regions
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    if (rank == 0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    // begin iterative solver
    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bbh_" + stage_text << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);

        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(ss, space, bconfig, conf, lapse, shift);
      }
      update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2, &syst);
      ite++;
    }

    // in the case of corotation or fixed local spin frequency, we calculate MCH accordingly
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)) ||
        bconfig.control(COROT_BIN))
      bconfig(MCH, BCO1) =
          bco_utils::syst_mch(syst, space, "intSm", space.BH1 + 2);
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)) ||
        bconfig.control(COROT_BIN))
      bconfig(MCH, BCO2) =
          bco_utils::syst_mch(syst, space, "intSp", space.BH2 + 2);

    bconfig.set_filename(converged_filename(stage_text, bconfig));
    if (rank == 0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  };
  if (stage_enabled[TESTING]) {
    test_stage(TESTING, "TESTING");
  }

  // clean-up vector fields

  MPI_Finalize();
  return EXIT_SUCCESS;
}

// Print useful diagnostics during iterative solver to monitor the process
template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics(space_t const& space,
                       syst_t const& syst,
                       int const ite,
                       double const conv,
                       config_t bconfig) {

  std::ios_base::fmtflags f(std::cout.flags());
  int ndom = space.get_nbr_domains();

  // local spin angular momentum BH1
  Val_domain integSm(syst.give_val_def("intSm")()(space.BH1 + 2));
  double Sm = space.get_domain(space.BH1 + 2)->integ(integSm, INNER_BC);

  // local spin angular momentum BH2
  Val_domain integSp(syst.give_val_def("intSp")()(space.BH2 + 2));
  double Sp = space.get_domain(space.BH2 + 2)->integ(integSp, INNER_BC);

  // irreducible mass BH1
  Val_domain integMmsq(syst.give_val_def("intMsq")()(space.BH1 + 2));
  double Mirrmsq = space.get_domain(space.BH1 + 2)->integ(integMmsq, INNER_BC);
  double Mirrm = std::sqrt(Mirrmsq);
  // MCH
  double Mchm = std::sqrt(Mirrmsq + Sm * Sm / 4. / Mirrmsq);

  // irreducible mass BH2
  Val_domain integMpsq(syst.give_val_def("intMsq")()(space.BH2 + 2));
  double Mirrpsq = space.get_domain(space.BH2 + 2)->integ(integMpsq, INNER_BC);
  double Mirrp = std::sqrt(Mirrpsq);
  // MCH
  double Mchp = std::sqrt(Mirrpsq + Sp * Sp / 4. / Mirrpsq);

  // ADM Linear Momenta
  Val_domain integPx(syst.give_val_def("intPx")()(ndom - 1));
  double Px = space.get_domain(ndom - 1)->integ(integPx, OUTER_BC);

  Val_domain integPy(syst.give_val_def("intPy")()(ndom - 1));
  double Py = space.get_domain(ndom - 1)->integ(integPy, OUTER_BC);

  Val_domain integPz(syst.give_val_def("intPz")()(ndom - 1));
  double Pz = space.get_domain(ndom - 1)->integ(integPz, OUTER_BC);

  cout << "=======================================" << endl;
#define FORMAT std::setw(12) << std::left << std::showpos
  std::cout << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl;
  std::cout << FORMAT << "Omega: " << bconfig(GOMEGA) << std::endl;
  if (!std::isnan(bconfig.set(ADOT)))
    std::cout << FORMAT << "Adot: " << bconfig(ADOT) << std::endl;
  std::cout << FORMAT << "Axis: " << bconfig(COM) << std::endl;
  std::cout << FORMAT << "P: " << "[" << Px << ", " << Py << ", " << Pz
            << "]\n\n";

  std::cout << FORMAT << "BH1-Mirr: " << Mirrm << std::endl;
  std::cout << FORMAT << "BH1-Mch: " << Mchm << std::endl;
  std::cout << FORMAT << "BH1-R: "
            << bco_utils::get_radius(space.get_domain(space.BH1 + 1), EQUI)
            << std::endl;
  std::cout << FORMAT << "BH1-S: " << Sm << std::endl;
  std::cout << FORMAT << "BH1-Chi: " << Sm / Mchm / Mchm << std::endl;
  std::cout << FORMAT << "BH1-Omega: " << bconfig(OMEGA, BCO1) << "\n\n";

  std::cout << FORMAT << "BH2-Mirr: " << Mirrp << std::endl;
  std::cout << FORMAT << "BH2-Mch: " << Mchp << std::endl;
  std::cout << FORMAT << "BH2-R: "
            << bco_utils::get_radius(space.get_domain(space.BH2 + 1), EQUI)
            << std::endl;
  std::cout << FORMAT << "BH2-S: " << Sp << std::endl;
  std::cout << FORMAT << "BH2-Chi: " << Sp / Mchp / Mchp << std::endl;
  std::cout << FORMAT << "BH2-Omega: " << bconfig(OMEGA, BCO2) << "\n\n";
  std::cout.flags(f);

  std::cout << "=======================================" << std::endl;
}

template <typename space_t, typename syst_t>
void print_diagnostics_pre(space_t const& space,
                           syst_t const& syst,
                           int const ite,
                           double const conv) {
  cout << "=======================================" << endl;
  cout << "Iter " << ite << " \t" << conv << endl;

  std::cout << "R \t"
            << bco_utils::get_radius(space.get_domain(space.BH1 + 1), OUTER_BC)
            << " "
            << bco_utils::get_radius(space.get_domain(space.BH2 + 1), OUTER_BC)
            << std::endl;
  cout << "=======================================" << endl;
}

/**
 * standardized filename for each converged dataset at the end of each stage.
 */
template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig) {
  std::stringstream ss;
  ss << "converged_BBH";
  if (stage != "")
    ss << "_" << stage << ".";
  else
    ss << ".";
  ss << bconfig(DIST) << "." << bconfig(CHI, BCO1) << "." << bconfig(CHI, BCO2)
     << "." << bconfig(MCH, BCO1) + bconfig(MCH, BCO2) << ".q"
     << bconfig(MCH, BCO2) / bconfig(MCH, BCO1) << "." << std::setfill('0')
     << std::setw(2) << bconfig(BIN_RES);
  return ss.str();
}
