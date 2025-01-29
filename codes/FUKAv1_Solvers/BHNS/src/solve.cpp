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
#include <sstream>
#include "Configurator/config_binary.hpp"
#include "EOS/EOS.hh"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath.hpp"
#include "mpi.h"

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

// forward declarations
template <class eos_t, typename config_t>
int BHNS_solver(config_t bconfig, std::string outputdir);

template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics(space_t const& space,
                       syst_t const& syst,
                       int const ite,
                       double const conv,
                       config_t& bconfig);

template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig);

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 2) {
    std::cerr << "Usage: ./solve /<path>/<ID base name>.info "
                 "./<output_path>/"
              << endl;
    std::cerr << "Ex: ./solve initbin.info ./out" << endl;
    std::cerr << "Note: <output_path> is optional and defaults to <path>"
              << endl;
    std::_Exit(EXIT_FAILURE);
  }
  std::string input_filename = argv[1];
  kadath_config_boost<BIN_INFO> bconfig(input_filename);

  std::string outputdir = bconfig.config_outputdir();
  if (argc > 2)
    outputdir = argv[2];

  // setup eos to update central density
  const double h_cut = bconfig.eos<double>(HCUT, BCO1);

  const std::string eos_file = bconfig.eos<std::string>(EOSFILE, BCO1);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);

  int return_status = EXIT_SUCCESS;
  if (eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut);
    return_status = BHNS_solver<eos_t>(bconfig, outputdir);
  } else if (eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS, BCO1) == 0)
                               ? 2000
                               : bconfig.eos<int>(INTERP_PTS, BCO1);

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut, interp_pts);
    return_status = BHNS_solver<eos_t>(bconfig, outputdir);
  } else {
    if (rank == 0)
      std::cerr << eos_type << " is not recognized.\n";
    std::_Exit(EXIT_FAILURE);
  }
  // end eos setup

  MPI_Finalize();

  return return_status;
}

template <class eos_t, typename config_t>
int BHNS_solver(config_t bconfig, std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << bconfig;

  bconfig(Q) = bconfig(MADM, BCO1) / bconfig(MCH, BCO2);

  std::string kadath_filename = bconfig.space_filename();

  FILE* fin = fopen(kadath_filename.c_str(), "r");
  Space_bhns space(fin);
  Scalar conf(space, fin);
  Scalar lapse(space, fin);
  Vector shift(space, fin);
  Scalar logh(space, fin);
  Scalar phi(space, fin);
  fclose(fin);

  if (bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  bconfig.set_outputdir(outputdir);
  auto& stage_enabled = bconfig.return_stages();
  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);
  Base_tensor basis(space, CARTESIAN_BASIS);

  int ndom = space.get_nbr_domains();

  double xc1 = bco_utils::get_center(space, space.NS);
  double xc2 = bco_utils::get_center(space, space.BH);
  double xo = bco_utils::get_center(space, ndom - 1);

  std::array<int, 2> excluded_doms{space.BH, space.BH + 1};

  double loghc = bco_utils::get_boundary_val(space.NS, logh, INNER_BC);

  Metric_flat fmet(space, basis);

  CoordFields<Space_bhns> cfields(space);

  //setup coord fields
  vec_ary_t coord_vectors = default_binary_vector_ary(space);
  update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);
  //end setup coord fields

  if (rank == 0) {
    std::cout << "=================================" << endl;
    std::cout << "BHNS grav input" << endl;
    std::cout << "Distance: " << bconfig(DIST) << endl;
    std::cout << "Omega guess: " << bconfig(GOMEGA) << endl;
    std::cout << "Units: " << bconfig(QPIG) << endl;
    std::cout << "=================================" << endl;
  }

  auto total_stages = [&](auto&& stage) {
    bconfig.set(MIRR, BCO2) =
        bco_utils::mirr_from_mch(bconfig(CHI, BCO2), bconfig(MCH, BCO2));
    if (rank == 0 && stage == TOTAL)
      std::cout << "############################################" << std::endl
                << "Total with mass fixing using\n"
                << "using fixed lapse BC on the BH\n"
                << "############################################" << std::endl;
    else if (rank == 0 && stage == TOTAL_BC)
      std::cout << "############################################" << std::endl
                << "Total syst with Py fixing while fixing the\n"
                << "BH lapse BC with von Neumann condition\n"
                << "############################################" << std::endl;

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    // EOS user defined operators
    Param p;
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
    syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

    // equation constants
    syst.add_cst("4piG", bconfig(QPIG));
    syst.add_cst("PI", M_PI);

    // constant vector fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);

    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // NS fixing parameters
    syst.add_cst("Madm1", bconfig(MADM, BCO1));
    syst.add_cst("Mb1", bconfig(MB, BCO1));
    syst.add_cst("chi1", bconfig(CHI, BCO1));

    // BH fixing parameters
    syst.add_cst("Mirr", bconfig(MIRR, BCO2));
    syst.add_cst("Mch", bconfig(MCH, BCO2));
    syst.add_cst("chi2", bconfig(CHI, BCO2));
    if (stage == TOTAL)
      syst.add_cst("n0", bconfig(FIXED_LAPSE, BCO2));

    // NS free variables
    syst.add_var("qlMadm1", bconfig(QLMADM, BCO1));
    syst.add_var("Hc1", loghc);

    // binary free variables
    syst.add_var("xaxis", bconfig(COM));
    // determine whether orbital frequency is fixed
    if (bconfig.control(FIXED_GOMEGA))
      syst.add_cst("ome", bconfig(GOMEGA));
    else
      syst.add_var("ome", bconfig(GOMEGA));

    // variables fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);
    syst.add_var("H", logh);

    // set omega to zero if we want corotation
    // the equations will simplify based on this
    // control.
    if (bconfig.control(COROT_BIN)) {
      bconfig.set(OMEGA, BCO1) = 0.;
      syst.add_cst("omes1", bconfig(OMEGA, BCO1));
      bconfig.set(OMEGA, BCO2) = 0.;
      syst.add_cst("omes2", bconfig(OMEGA, BCO2));
    } else {
      // with mixed spins, the velocity potential, phi,
      // becomes a variable field
      syst.add_var("phi", phi);

      // it's more meaningful to fix based on chi, however,
      // the capability exists to fix the spin based on
      // a fixed value of omega on the compact object
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst("omes1", bconfig(FIXED_BCOMEGA, BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
      } else {
        syst.add_var("omes1", bconfig(OMEGA, BCO1));
      }

      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst("omes2", bconfig(FIXED_BCOMEGA, BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
      } else {
        syst.add_var("omes2", bconfig(OMEGA, BCO2));
      }

      // for mixed spins, we have additional definitions
      // that are required for both objects that describe
      // the local rotation field
      for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      syst.add_def(space.ADAPTEDBH + 1, "s^i = omes2 * mp^i");
    }

    // convenience definitions
    syst.add_def("h = exp(H)");
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("dHdlnrho = dHdlnrho(h)");
    syst.add_def("delta = h - eps - 1.");
    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    if (stage == TOTAL) {
      syst.add_cst("yaxis", bconfig(COMY));
      syst.add_def("Morb^i = mg^i + xaxis * ey^i");
    } else {
      syst.add_var("yaxis", bconfig(COMY));
      syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");
    }
    syst.add_def("B^i = bet^i + ome * Morb^i");

    syst.add_def(
        "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i");
    syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i");
    syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i");

    syst.add_def(space.ADAPTEDNS + 1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
    syst.add_def(space.ADAPTEDBH + 1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");

    syst.add_def("intMsq  = P^4 / 4. / 4piG");

    for (int d = 0; d < ndom; d++) {
      if (d >= space.ADAPTEDNS + 1) {
        if (!bconfig.control(COROT_BIN))
          syst.add_eq_full(d, "phi= 0");

        syst.add_eq_full(d, "H  = 0");

        syst.add_def(d, "eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst.add_def(d,
                     "eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,
                     "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij "
                     "* D_j Ntilde");

      } else {
        if (bconfig.control(COROT_BIN)) {
          syst.add_def(d, "U^i    = B^i / N");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "firstint = log(h * N / W)");
        } else {
          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "V^i    = N * U^i - B^i");
          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
          syst.add_def(d,
                       "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * "
                       "W * V^i)");
        }

        syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
        syst.add_def(
            d,
            "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare");
        syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

        syst.add_def(d,
                     "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * "
                     "delta + 4piG / 2. * P^5 * Etilde");
        syst.add_def(d,
                     "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta "
                     "* A_ij *A^ij "
                     "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst.add_def(
            d,
            "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
            "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
            "ptilde^i");

        syst.add_def(d, "intMb  = P^6 * rho * W");
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");
      }
    }

    space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq(syst, "eqP= 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

    syst.add_eq_bc(space.ADAPTEDNS, OUTER_BC, "H = 0");
    if (!bconfig.control(COROT_BIN)) {
      syst.add_eq_bc(space.ADAPTEDNS, OUTER_BC, "V^i * D_i H = 0");

      for (int i = space.NS; i < space.ADAPTEDNS; ++i) {
        syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching(i, OUTER_BC, "phi");
        syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot(space.ADAPTEDNS, 2, "eqphi = 0", "phi=0");

      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_outer_NS(syst, "integ(intS1) / Madm1 / Madm1 = chi1");

      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + s^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_BH(syst, "integ(intS2) - chi2 * Mch * Mch = 0 ");
    } else {
      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");
    }

    Index posori1(space.get_domain(space.NS)->get_nbr_points());
    syst.add_eq_val(space.NS, "ex^i * D_i H", posori1);

    syst.add_eq_first_integral(space.NS, space.ADAPTEDNS, "firstint",
                               "H - Hc1");

    space.add_eq_int_volume(syst, space.NS, space.ADAPTEDNS,
                            "integvolume(intM) = qlMadm1");
    space.add_eq_int_volume(syst, space.NS, space.ADAPTEDNS,
                            "integvolume(intMb) = Mb1");

    // determine whether orbital frequency is fixed
    if (!bconfig.control(FIXED_GOMEGA))
      space.add_eq_int_inf(syst, "integ(intPy) = 0");

    if (stage == TOTAL) {
      //For TOTAL we use the fixed lapse condition
      space.add_bc_sphere_two(syst, "N = n0");
    } else {
      //For TOTAL_BC we use the von Neumann condition
      //and constrain Padm^x = 0
      space.add_eq_int_inf(syst, "integ(intPx) = 0");
      space.add_bc_sphere_two(syst, "sp^j * D_j NP = 0");
    }

    space.add_bc_sphere_two(
        syst,
        "sp^j * D_j P + P / 4 * D^j sp_j + A_ij * sp^i * sp^j / P^3 / 4 = 0");
    space.add_eq_int_BH(syst, "integ(intMsq) - Mirr * Mirr = 0 ");

    // excised domains of the BH
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    if (rank == 0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

      std::stringstream ss;
      ss << "total_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);

        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh,
                                  phi);
      }
      ite++;
    }

    // generate stage string for converged filename
    std::string s{"TOTAL"};
    if (stage == TOTAL_BC)
      s += "_BC";
    if (bconfig.control(COROT_BIN))
      s += "_COROT";
    bconfig.set_filename(converged_filename(s, bconfig));
    if (rank == 0) {
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
  };  //end total_stages lambda
  auto rescaling_stages = [&](auto&& stage, std::string&& outputstr) {
    if (stage == ECC_RED) {
      if (rank == 0)
        std::cout << "############################################" << std::endl
                  << "Eccentricity reduction step with fixed omega" << std::endl
                  << "############################################"
                  << std::endl;

      // Condition to set PN estimates of omega and adot
      if (std::isnan(bconfig.set(ADOT)) || std::isnan(bconfig.set(ECC_OMEGA)) ||
          bconfig.control(USE_PN)) {
        bco_utils::KadathPNOrbitalParams(bconfig, bconfig(MADM, BCO1),
                                         bconfig(MCH, BCO2));

        if (rank == 0)
          std::cout << "### Using PN estimate for adot and omega! ###"
                    << std::endl;
      }
      bconfig.set(GOMEGA) = bconfig(ECC_OMEGA);
    } else if (rank == 0)
      std::cout << "############################################" << std::endl
                << "Hydro-rescaling with fixed omega" << std::endl
                << "############################################" << std::endl;

    // rescaling constant for the log specific enthalpy. H
    double H_scale = 0;

    // make a copy of the matter solution since this is rescaled by a constant
    // factor H_scale
    Scalar logh_const(logh);
    logh_const.std_base();

    // setup the "radial" vector field separating the coordinate centers of
    // each compact object to the center-of-mass of the system.  since
    // each object is centered at (+/-x,0,0), we just need the cartesian
    // coordinates shifted by COM and COMY
    Vector CART(space, CON, basis);
    CART = cfields.cart();

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    // EOS user defined operators
    Param p;
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
    syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

    // equation constants
    syst.add_cst("4piG", bconfig(QPIG));
    syst.add_cst("PI", M_PI);

    // constant vector fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);

    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // NS fixing parameters
    syst.add_cst("Madm1", bconfig(MADM, BCO1));
    syst.add_cst("Mb1", bconfig(MB, BCO1));
    syst.add_cst("chi1", bconfig(CHI, BCO1));

    // BH fixing parameters
    syst.add_cst("Mirr", bconfig(MIRR, BCO2));
    syst.add_cst("Mch", bconfig(MCH, BCO2));
    syst.add_cst("chi2", bconfig(CHI, BCO2));

    // binary fixing parameters
    syst.add_cst("ome", bconfig(GOMEGA));

    // constant fields
    syst.add_cst("Hconst", logh_const);

    // NS free variable - rescaling parameter
    syst.add_var("Hscale", H_scale);

    // binary free variables
    syst.add_var("xaxis", bconfig(COM));
    syst.add_var("yaxis", bconfig(COMY));

    // variables fields
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    // set omega to zero if we want corotation
    // the equations will simplify based on this
    // control.
    if (bconfig.control(COROT_BIN)) {
      bconfig.set(OMEGA, BCO1) = 0.;
      syst.add_cst("omes1", bconfig(OMEGA, BCO1));
      bconfig.set(OMEGA, BCO2) = 0.;
      syst.add_cst("omes2", bconfig(OMEGA, BCO2));
    } else {
      // with mixed spins, the velocity potential, phi,
      // becomes a variable field
      syst.add_var("phi", phi);

      // it's more meaningful to fix based on chi, however,
      // the capability exists to fix the spin based on
      // a fixed value of omega on the compact object
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst("omes1", bconfig(FIXED_BCOMEGA, BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
      } else {
        syst.add_var("omes1", bconfig(OMEGA, BCO1));
      }

      if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst("omes2", bconfig(FIXED_BCOMEGA, BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
      } else {
        syst.add_var("omes2", bconfig(OMEGA, BCO2));
      }

      // for mixed spins, we have additional definitions
      // that are required for both objects that describe
      // the local rotation field
      for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      syst.add_def(space.ADAPTEDBH + 1, "s^i = omes2 * mp^i");
    }

    // definition of H
    // - rescaled in the star
    // - constant everywhere else
    for (int d = 0; d < ndom; ++d) {
      if (d <= space.ADAPTEDNS && d >= space.NS)
        syst.add_def(d, "H  = Hconst * (1. + Hscale)");
      else
        syst.add_def(d, "H  = Hconst");
    }

    // convenience definitions
    syst.add_def("h = exp(H)");
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("dHdlnrho = dHdlnrho(h)");
    syst.add_def("delta = h - eps - 1.");
    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    // definitions of the orbital velocity field and its
    // contribution to the shift
    syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");

    std::string bigB{"B^i = bet^i + ome * Morb^i"};
    if (stage == ECC_RED) {
      syst.add_cst("adot", bconfig(ADOT));
      syst.add_cst("r", CART);
      syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
      // add contribution of ADOT to the total shift definition
      bigB += " + adot * comr^i";
    }

    // full shift vector including inertial + orbital contributions
    syst.add_def(bigB.c_str());

    syst.add_def(
        "A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i");
    syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i");
    syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i");

    syst.add_def(ndom - 1, "COMx  = -3 * P^4 * einf^i * ex_i / 8. / PI");
    syst.add_def(ndom - 1, "COMy  =  3 * P^4 * einf^i * ey_i / 8. / PI");
    syst.add_def(ndom - 1, "COMz  =  3 * P^4 * einf^i * ez_i / 8. / PI");

    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");

    syst.add_def(space.ADAPTEDNS + 1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
    syst.add_def(space.ADAPTEDBH + 1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");

    syst.add_def("intMsq  = P^4 / 4. / 4piG");

    for (int d = 0; d < ndom; d++) {
      // vacuum constraint equations definitions
      if (d >= space.ADAPTEDNS + 1) {
        if (!bconfig.control(COROT_BIN))
          syst.add_eq_full(d, "phi= 0");

        syst.add_def(d, "eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst.add_def(d,
                     "eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,
                     "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij "
                     "* D_j Ntilde");

      } else {
        //matter constraint equations definitions
        if (bconfig.control(COROT_BIN)) {
          syst.add_def(d, "U^i    = B^i / N");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "firstint = log(h * N / W)");
        } else {
          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "V^i    = N * U^i - B^i");
          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
          syst.add_def(d,
                       "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * "
                       "W * V^i)");
        }

        syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
        syst.add_def(
            d,
            "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare");
        syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

        syst.add_def(d,
                     "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * "
                     "delta + 4piG / 2. * P^5 * Etilde");
        syst.add_def(d,
                     "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta "
                     "* A_ij *A^ij "
                     "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst.add_def(
            d,
            "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
            "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
            "ptilde^i");

        syst.add_def(d, "intMb  = P^6 * rho * W");
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");
      }
    }

    space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq(syst, "eqP= 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

    syst.add_eq_bc(space.ADAPTEDNS, OUTER_BC, "H = 0");

    // determine equations to add based on corotation
    // or mixed spin binaries based on chi or fixed omega
    if (!bconfig.control(COROT_BIN)) {
      syst.add_eq_bc(space.ADAPTEDNS, OUTER_BC, "V^i * D_i H = 0");

      for (int i = space.NS; i < space.ADAPTEDNS; ++i) {
        syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching(i, OUTER_BC, "phi");
        syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot(space.ADAPTEDNS, 2, "eqphi = 0", "phi=0");

      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_outer_NS(syst, "integ(intS1) / Madm1 / Madm1 = chi1");

      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + s^i");
      if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_BH(syst, "integ(intS2) - chi2 * Mch * Mch = 0 ");
    } else {
      space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");
    }

    space.add_eq_int_volume(syst, space.NS, space.ADAPTEDNS,
                            "integvolume(intMb) = Mb1");

    space.add_eq_int_inf(syst, "integ(intPx) = 0");
    space.add_eq_int_inf(syst, "integ(intPy) = 0");

    space.add_bc_sphere_two(syst, "sp^j * D_j NP = 0");
    space.add_bc_sphere_two(
        syst,
        "sp^j * D_j P + P / 4 * D^j sp_j + A_ij * sp^i * sp^j / P^3 / 4 = 0");

    space.add_eq_int_BH(syst, "integ(intMsq) - Mirr * Mirr = 0 ");

    // excised domains of the BH
    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    if (rank == 0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      // overwrite logh with scaled version
      logh = syst.give_val_def("H");

      update_fields(cfields, coord_vectors, {}, xo, xc1, xc2, &syst);

      // manual update to CART field
      CART = cfields.cart();
      for (int d = 0; d < ndom; ++d)
        update_field(syst, d, "CART^i", CART);
      syst.sec_member();

      std::stringstream ss;
      ss << "rescaling_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);

        std::cout << H_scale << " " << std::endl;

        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh,
                                  phi);
      }

      ite++;
    }
    bconfig.set_filename(converged_filename(outputstr, bconfig));
    if (rank == 0) {
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
  };

  // solve using fixed lapse BC for the BH, but without Py fixing
  if (stage_enabled[TOTAL]) {
    total_stages(TOTAL);
    if (last_stage_idx != TOTAL)
      bconfig.set_stage(TOTAL) = false;
  }

  if (stage_enabled[TOTAL_BC] && bconfig.control(FIXED_GOMEGA)) {
    rescaling_stages(TOTAL_BC, "TOTAL_BC_FIXED_OMEGA");
    bconfig.control(FIXED_GOMEGA) = false;
  }

  // solve using von Neumann lapse BC for the BH and use Py fixing
  if (stage_enabled[TOTAL_BC])
    total_stages(TOTAL_BC);

  // solve using von Neumann lapse BC with fixed orbital omega GOMEGA
  // and radial infall velocity ADOT.
  // Since omega is fixed, the matter solution is only rescaled.
  // Changing resolution and/or properties of the compact objects
  // (e.g. mass, spin) requires rerunning TOTAL_BC before ECC_RED!
  if (stage_enabled[ECC_RED]) {
    rescaling_stages(ECC_RED, "ECC_RED");
  }  //end Ecc_red stage

  bconfig.set_filename(converged_filename("", bconfig));
  if (rank == 0) {
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
  }
  return EXIT_SUCCESS;
}  //end BHNS solver

template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics(space_t const& space,
                       syst_t const& syst,
                       int const ite,
                       double const conv,
                       config_t& bconfig) {

  std::ios_base::fmtflags f(std::cout.flags());
  int ndom = space.get_nbr_domains();
  double baryonic_mass1 = 0.;
  double mass1 = 0.;

  for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
    baryonic_mass1 += syst.give_val_def("intMb")()(d).integ_volume();
    mass1 += syst.give_val_def("intM")()(d).integ_volume();
  }
  auto [rmin, rmax] = bco_utils::get_rmin_rmax(space, space.ADAPTEDNS);
  double r_bh = bco_utils::get_radius(space.get_domain(space.ADAPTEDBH), EQUI);

  Val_domain integS1(syst.give_val_def("intS1")()(space.ADAPTEDNS + 1));
  double S1 = space.get_domain(space.ADAPTEDNS + 1)->integ(integS1, OUTER_BC);
  double chi1 = S1 / bconfig(MADM, BCO1) / bconfig(MADM, BCO1);

  Val_domain integS2(syst.give_val_def("intS2")()(space.ADAPTEDBH + 1));
  double S2 = space.get_domain(space.ADAPTEDBH + 1)->integ(integS2, OUTER_BC);

  Val_domain integPx(syst.give_val_def("intPx")()(ndom - 1));
  double Px = space.get_domain(ndom - 1)->integ(integPx, OUTER_BC);

  Val_domain integPy(syst.give_val_def("intPy")()(ndom - 1));
  double Py = space.get_domain(ndom - 1)->integ(integPy, OUTER_BC);

  Val_domain integPz(syst.give_val_def("intPz")()(ndom - 1));
  double Pz = space.get_domain(ndom - 1)->integ(integPz, OUTER_BC);

  double mirrsq =
      space.get_domain(space.BH + 2)
          ->integ(syst.give_val_def("intMsq")()(space.BH + 2), INNER_BC);
  double Mirr = std::sqrt(mirrsq);
  double Mch = std::sqrt(mirrsq + S2 * S2 / 4. / mirrsq);
#define FORMAT std::setw(13) << std::left << std::showpos
  std::cout << "=======================================" << endl;
  std::cout << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << "\n";
  std::cout << FORMAT << "Omega: " << bconfig(GOMEGA) << endl;
  std::cout << FORMAT << "Axis: " << bconfig(COM) << endl;
  std::cout << FORMAT << "Padm: " << "[" << Px << ", " << Py << ", " << Pz
            << "]\n\n";

  std::cout << FORMAT << "NS-Mb: " << baryonic_mass1 << std::endl;
  std::cout << FORMAT << "NS-Madm_ql: " << mass1 << std::endl;
  std::cout << FORMAT << "NS-R: " << rmin << " " << rmax << std::endl
            << FORMAT << "NS-S: " << S1 << std::endl;
  std::cout << FORMAT << "NS-Chi: " << chi1 << std::endl;
  std::cout << FORMAT << "NS-Omega: " << bconfig(OMEGA, BCO1) << std::endl
            << std::endl;

  std::cout << FORMAT << "BH-Mirr: " << Mirr << std::endl;
  std::cout << FORMAT << "BH-Mch: " << Mch << std::endl;
  std::cout << FORMAT << "BH-R: " << r_bh << std::endl;
  std::cout << FORMAT << "BH-S: " << S2 << std::endl;
  std::cout << FORMAT << "BH-Chi: " << S2 / Mch / Mch << std::endl;
  std::cout << FORMAT << "BH-Omega: " << bconfig(OMEGA, BCO2) << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << endl;
}  //end print diagnostics

//standardized filename for each converged dataset at the end of each stage.
template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig) {
  std::stringstream ss;
  ss << "converged_BHNS";
  if (stage != "")
    ss << "_" << stage << ".";
  else
    ss << ".";
  ss << bconfig(DIST) << "." << bconfig(CHI, BCO1) << "." << bconfig(CHI, BCO2)
     << "." << bconfig(MADM, BCO1) + bconfig(MCH, BCO2) << ".q" << bconfig(Q)
     << "." << std::setfill('0') << std::setw(2) << bconfig(BIN_RES);
  return ss.str();
}  //end converged filename
