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
#include <iterator>
#include <sstream>
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath_adapted_bh.hpp"
#include "mpi.h"

//config_file includes
#include "Configurator/config_bco.hpp"
using namespace Kadath;
using namespace Kadath::FUKA_Config;

template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics(space_t const& space,
                       syst_t const& syst,
                       int const ite,
                       double const conv,
                       config_t& bconfig);

template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics_pre(space_t const& space,
                           syst_t const& syst,
                           int const ite,
                           double const conv,
                           config_t& bconfig,
                           double S = 0);

template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig);

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
    std::cerr << "Ex: ./solve initbh.info ./out" << endl;
    std::cerr << "Note: <output_path> is optional and defaults to <path>"
              << endl;
    std::_Exit(EXIT_FAILURE);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /****************Reading from Config file ************************************/
  input_filename = argv[1];
  kadath_config_boost<BCO_BH_INFO> bconfig(input_filename);

  std::array<bool, NUM_BCO_FIELDS> fields = bconfig.return_fields();
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

  if (rank == 0) {
    std::cout << "Last Stage Enabled: " << last_stage << std::endl;
  }
  /****************************************************************************/

  //file containing KADATH fields must have same name as
  //config file with only the extension being different
  std::string kadath_filename = bconfig.space_filename();

  FILE* fin = fopen(kadath_filename.c_str(), "r");
  Space_adapted_bh space(fin);
  Base_tensor basis(space, CARTESIAN_BASIS);
  Scalar conf(space, fin);
  Scalar lapse(space, fin);
  Vector shift(space, fin);
  fclose(fin);
  if (argc > 2)
    bconfig.set_outputdir(argv[2]);

  Metric_flat fmet(space, basis);

  int ndom = space.get_nbr_domains();

  //exclude excision interior from solver
  std::array<int, 2> excluded_doms{0, 1};

  double xo = bco_utils::get_center(space, 0);

  //setup coord fields
  CoordFields<Space_adapted_bh> cf_generator(space);

  vec_ary_t coord_vectors{default_co_vector_ary(space)};

  scalar_ary_t coord_scalars;
  coord_scalars[R_BCO1] = Scalar(space);

  update_fields_co(cf_generator, coord_vectors, coord_scalars, xo);

  //end setup coord fields

  if (rank == 0) {
    std::cout << "=================================" << endl;
    std::cout << bconfig;
    std::cout << "Radius: " << bconfig(RMID) << endl;
    std::cout << "Pos: " << xo << endl;
    std::cout << "=================================" << endl;
  }
  if (stage_enabled[PRE]) {
    if (rank == 0)
      std::cout << "############################" << std::endl
                << "Grav lapse and psi preconditioning" << std::endl
                << "############################" << std::endl;

    Scalar level(space);
    level = (*coord_scalars[R_BCO1]) * (*coord_scalars[R_BCO1]) -
            bconfig(RMID) * bconfig(RMID);
    level.std_base();

    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");
    syst.add_cst("PI", M_PI);

    syst.add_cst("n0", bconfig(FIXED_LAPSE));
    syst.add_cst("lev", level);

    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("rm", *coord_scalars[R_BCO1]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    syst.add_var("P", conf);
    syst.add_var("N", lapse);

    syst.add_def("NP  = P*N");
    syst.add_def("eqP = D^i D_i P");
    syst.add_def("eqNP= D^i D_i NP");
    syst.add_def("intMsq= P^4 / 16. / PI");

    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");

    space.add_bc_inf(syst, "NP = 1");
    space.add_bc_inf(syst, "P = 1");

    space.add_bc_bh(syst, "N = n0");
    space.add_bc_bh(syst, "sm^j * D_j P + P / 4 * D^j sm_j = 0");

    space.add_eq_int_bh(syst, "integ(lev) = 0");

    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
    }
    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "grav_np_pre_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics_pre(space, syst, ite, conv, bconfig);
        std::cout << "=======================================" << "\n\n";
        if (bconfig.control(CHECKPOINT)) {
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
        }
      }
      ite++;
    }
    /**
     * Set irreducible mass based on fixed radius
     * and Mch based on Mirr and chi
     */
    if (bconfig.control(USE_FIXED_R)) {
      double mirrsq =
          space.get_domain(2)->integ(syst.give_val_def("intMsq")()(2),
                                     INNER_BC);
      bconfig.set(MIRR) = std::sqrt(mirrsq);
      bconfig.set(MCH) =
          std::sqrt(2. / (1 + (1 - bconfig(CHI) * bconfig(CHI)))) *
          bconfig(MIRR);
    }
    if (last_stage_idx != PRE)
      bconfig.set_stage(PRE) = false;

    bconfig.set_filename(converged_filename("pre", bconfig));
    if (rank == 0) {
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
    }
  }  //PRE stage

  if (stage_enabled[TOTAL]) {
    if (rank == 0)
      std::cout << "###################################" << std::endl
                << "Solving using Fixed Lapse condition" << std::endl
                << "###################################" << std::endl;

    //We always fix based on MCH
    bconfig.set(MIRR) = bco_utils::mirr_from_mch(bconfig(CHI), bconfig(MCH));
    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);

    syst.add_cst("n0", bconfig(FIXED_LAPSE));
    syst.add_cst("M", bconfig(MIRR));
    syst.add_cst("CM", bconfig(MCH));
    syst.add_cst("chi", bconfig(CHI));
    syst.add_var("ome", bconfig(OMEGA));

    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    syst.add_def(
        "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    syst.add_def("intPx = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz = A_ij * ez^j * einf^i / 8 / PI");

    syst.add_def("intS  = A_ij * mg^i * sm^j / 8. / PI");
    syst.add_def("intMsq= P^4 / 16. / PI");

    syst.add_def("eqP   = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP  = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");

    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP  = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    space.add_bc_inf(syst, "NP  = 1");
    space.add_bc_inf(syst, "P   = 1");
    space.add_bc_inf(syst, "bet^i = 0");

    space.add_bc_bh(syst, "N = n0");
    space.add_bc_bh(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_bc_bh(syst, "bet^i = N / P^2 * sm^i + ome * mg^i");

    space.add_eq_int_bh(syst, "integ(intS) - chi * CM * CM = 0 ");
    space.add_eq_int_bh(syst, "integ(intMsq) - M * M = 0 ");

    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }

    bool endloop = false;
    int ite = 1;
    double conv;
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bh_fixed_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);
        if (bconfig.control(CHECKPOINT)) {
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
        }
      }

      update_fields_co(cf_generator, coord_vectors, {}, xo, &syst);
      ite++;
    }
    {
      bconfig.set(RMID) = bco_utils::get_radius(space.get_domain(1), OUTER_BC);
      bconfig.set_filename(converged_filename("total", bconfig));
      if (rank == 0)
        bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
      if (last_stage_idx != FIXED_OMEGA)
        bconfig.set_stage(FIXED_OMEGA) = false;
    }
  }  //end TOTAL Stage

  if (stage_enabled[TOTAL_BC]) {
    if (rank == 0)
      std::cout << "############################" << std::endl
                << "Total system with BC" << std::endl
                << "############################" << std::endl;
    bconfig.set(MIRR) = bco_utils::mirr_from_mch(bconfig(CHI), bconfig(MCH));
    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);

    syst.add_cst("M", bconfig(MIRR));
    syst.add_cst("chi", bconfig(CHI));
    syst.add_cst("CM", bconfig(MCH));
    syst.add_cst("xboost", bconfig(BVELX));
    syst.add_cst("yboost", bconfig(BVELY));

    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("rm", *coord_scalars[R_BCO1]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    syst.add_var("ome", bconfig(OMEGA));
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");
    syst.add_def(ndom - 1, "Mk   =  dr(N) / 4 / PI");

    syst.add_def(
        "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    syst.add_def("intPx = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz = A_ij * ez^j * einf^i / 8 / PI");

    syst.add_def("intS = A_ij * mg^i * sm^j / 8. / PI");
    syst.add_def("intMsq = P^4 / 16. / PI");

    syst.add_def("eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");

    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    space.add_bc_inf(syst, "NP = 1");
    space.add_bc_inf(syst, "P = 1");
    space.add_bc_inf(syst, "bet^i = 0");

    space.add_bc_bh(syst, "dn(NP) = 0");
    space.add_bc_bh(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_bc_bh(syst, "bet^i = N / P^2 * sm^i + ome * mg^i");

    space.add_eq_int_bh(syst, "integ(intS) - chi * CM * CM = 0 ");
    space.add_eq_int_bh(syst, "integ(intMsq) - M * M = 0 ");

    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }
    bool endloop = false;
    int ite = 1;
    double conv;
    if (rank == 0)
      print_diagnostics(space, syst, ite, conv, bconfig);
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bh_total_bc_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);
        if (bconfig.control(CHECKPOINT)) {
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
        }
      }

      update_fields_co(cf_generator, coord_vectors, {}, xo, &syst);
      ite++;
    }
    {
      bconfig.set(RMID) = bco_utils::get_radius(space.get_domain(1), OUTER_BC);
      ;
      bconfig(FIXED_LAPSE) = bco_utils::get_boundary_val(2, lapse, INNER_BC);
      bconfig.set_filename(converged_filename("total_bc", bconfig));
      if (rank == 0)
        bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
    }
  }  //end TOTAL_BC Stage
  if (stage_enabled[TESTING]) {
    if (rank == 0)
      std::cout << "############################" << std::endl
                << "TESTING" << std::endl
                << "############################" << std::endl;
    bconfig.set(MIRR) = bco_utils::mirr_from_mch(bconfig(CHI), bconfig(MCH));
    coord_vectors[BCO1_ROT] = Vector(space, CON, basis);
    // MCH = 6, MNS = 1.5, separation 60km
    //const double dist = 40.6229;
    //const double com = 9.14015;

    // MCH1 = 0.91, MCH2 = 0.091, dist = 10
    const double dist = 10;
    const double com = 4.0646536239186215;

    // coord transformation to center or orbital rotation
    const double x_Op = 0.909091;
    double omega_boost = 0.0276269;
    update_fields(cf_generator, coord_vectors, {}, x_Op, xo, 0.);
    System_of_eqs syst(space, 0, ndom - 1);

    fmet.set_system(syst, "f");

    syst.add_cst("PI", M_PI);

    syst.add_cst("M", bconfig(MIRR));
    syst.add_cst("chi", bconfig(CHI));
    syst.add_cst("CM", bconfig(MCH));
    syst.add_cst("yboost", omega_boost);  //omega boost really - placeholder

    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("rm", *coord_scalars[R_BCO1]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    syst.add_var("ome", bconfig(OMEGA));
    syst.add_var("P", conf);
    syst.add_var("N", lapse);
    syst.add_var("bet", shift);

    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");
    syst.add_def(ndom - 1, "Mk   =  dr(N) / 4 / PI");
    syst.add_def("B^i = bet^i + yboost * mg^i");

    syst.add_def(
        "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    syst.add_def("intPx = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz = A_ij * ez^j * einf^i / 8 / PI");

    syst.add_def("intS = A_ij * mm^i * sm^j / 8. / PI");
    syst.add_def("intMsq = P^4 / 16. / PI");

    syst.add_def("eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");

    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    space.add_bc_inf(syst, "NP = 1");
    space.add_bc_inf(syst, "P = 1");
    space.add_bc_inf(syst, "bet^i = 0");

    space.add_bc_bh(syst, "dn(NP) = 0");
    space.add_bc_bh(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_bc_bh(syst, "B^i = N / P^2 * sm^i + ome * mm^i");

    space.add_eq_int_bh(syst, "integ(intS) - chi * CM * CM = 0 ");
    space.add_eq_int_bh(syst, "integ(intMsq) - M * M = 0 ");

    for (int i : excluded_doms) {
      syst.add_eq_full(i, "N = 0");
      syst.add_eq_full(i, "P = 0");
      syst.add_eq_full(i, "bet^i = 0");
    }
    bool endloop = false;
    int ite = 1;
    double conv;
    if (rank == 0)
      print_diagnostics(space, syst, ite, conv, bconfig);
    while (!endloop) {
      endloop = syst.do_newton(1e-8, conv);

      std::stringstream ss;
      ss << "bh_testing_" << ite - 1;
      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(space, syst, ite, conv, bconfig);
        if (bconfig.control(CHECKPOINT)) {
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
        }
      }

      update_fields(cf_generator, coord_vectors, {}, x_Op, xo, 0., &syst);
      ite++;
    }
    {
      bconfig.set(RMID) = bco_utils::get_radius(space.get_domain(1), OUTER_BC);
      ;
      bconfig(FIXED_LAPSE) = bco_utils::get_boundary_val(2, lapse, INNER_BC);
      bconfig.set_filename(converged_filename("testing", bconfig));
      if (rank == 0)
        bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
    }
  }  //end TESTING Stage

  bconfig.set_filename(converged_filename("", bconfig));
  if (rank == 0) {
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}

// runtime diagnostics specific for rotating solutions
template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics(space_t const& space,
                       syst_t const& syst,
                       int const ite,
                       double const conv,
                       config_t& bconfig) {

  int ndom = space.get_nbr_domains();

  Val_domain integS(syst.give_val_def("intS")()(2));
  double S = space.get_domain(2)->integ(integS, INNER_BC);
  print_diagnostics_pre(space, syst, ite, conv, bconfig, S);

  Val_domain integMsq(syst.give_val_def("intMsq")()(2));
  double Mirrsq = space.get_domain(2)->integ(integMsq, INNER_BC);
  double Mirr = std::sqrt(Mirrsq);
  double Mch = std::sqrt(Mirrsq + S * S / 4. / Mirrsq);

  std::cout << "Omega   " << bconfig(OMEGA) << std::endl;
  std::cout << "S       " << S << std::endl;
  std::cout << "Chi     " << S / Mch / Mch << " [" << bconfig(CHI) << "]"
            << std::endl;
  std::cout << "=======================================" << "\n\n";
}

// diagnostics at runtime for pre-conditioning stage
template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics_pre(space_t const& space,
                           syst_t const& syst,
                           int const ite,
                           double const conv,
                           config_t& bconfig,
                           double S) {

  int ndom = space.get_nbr_domains();
  double r = bco_utils::get_radius(space.get_domain(1), OUTER_BC);

  Val_domain integMsq(syst.give_val_def("intMsq")()(2));
  double Mirrsq = space.get_domain(2)->integ(integMsq, INNER_BC);
  double Mirr = std::sqrt(Mirrsq);
  double Mch = std::sqrt(Mirrsq + S * S / 4. / Mirrsq);

  std::cout << "=======================================" << std::endl;
  std::cout << "Iter    " << ite << " \tConv: " << conv << std::endl;
  std::cout << "R       " << r << std::endl;
  std::cout << "Mirr    " << Mirr << " [" << bconfig(MIRR) << "]" << std::endl;
  std::cout << "Mch     " << Mch << " [" << bconfig(MCH) << "]" << std::endl;
}

//standardized filename for each converged dataset at the end of each stage.
template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig) {
  std::stringstream ss;
  ss << "converged_BH";
  if (stage != "")
    ss << "_" << stage << ".";
  else
    ss << ".";
  ss << bconfig(MCH) << "." << bconfig(CHI) << "." << std::setfill('0')
     << std::setw(2) << bconfig(BCO_RES);
  return ss.str();
}
