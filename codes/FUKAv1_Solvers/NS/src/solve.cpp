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
#include <sstream>
#include "Configurator/config_bco.hpp"
#include "EOS/EOS.hh"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath_adapted.hpp"
#include "mpi.h"

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

// forward declarations
template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics_norot(space_t const& space,
                             syst_t const& syst,
                             config_t& bconfig,
                             int const ite,
                             double const conv);

template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics_rot(space_t const& space,
                           syst_t const& syst,
                           config_t& bconfig,
                           int const ite,
                           double const conv);

template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig);

template <class eos, typename config_t>
int NS_solver_3d(config_t bconfig, std::string outputdir);

// end forward declarations

int main(int argc, char** argv) {
  // initialize MPI
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    std::cerr << "Error starting MPI" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  // expecting cmd arguments
  if (argc < 2) {
    std::cerr << "Usage: ./solve /<path>/<ID base name>.info "
                 "./<output_path>/"
              << endl;
    std::cerr << "Ex: ./solve init-3d.info ./out" << endl;
    std::cerr << "Note: <output_path> is optional and defaults to <path>"
              << endl;
    std::_Exit(EXIT_FAILURE);
  }

  // load configuration file
  std::string ifilename{argv[1]};
  kadath_config_boost<BCO_NS_INFO> bconfig(ifilename);

  // specify output directory
  std::string outputdir = bconfig.config_outputdir();
  if (argc > 2)
    outputdir = argv[2];

  // load and setup the EOS
  const double h_cut = bconfig.eos<double>(HCUT);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

  if (eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut);
    return NS_solver_3d<eos_t>(bconfig, outputdir);
  } else if (eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0)
                               ? 2000
                               : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut, interp_pts);
    return NS_solver_3d<eos_t>(bconfig, outputdir);
  } else {
    std::cerr << "Unknown EOSTYPE." << endl;
    std::_Exit(EXIT_FAILURE);
  }

  MPI_Finalize();

  return EXIT_FAILURE;
}  // end main()

template <class eos_t, typename config_t>
int NS_solver_3d(config_t bconfig, std::string outputdir) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence threshold
  // FIXME should this be part of the config?
  double conv_thres = 1e-8;
  // physical origin of the system
  double xo = 0.;
  // logarithm of the central enthalpy, a variable in the system of equations
  double loghc = std::log(bconfig(HC));

  // specify which stages are enabled, pulled from the configuration file
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
  FILE* ff1 = fopen(spacein.c_str(), "r");
  Space_spheric_adapted space(ff1);

  // load the fields defined on the space
  Scalar conf(space, ff1);
  Scalar lapse(space, ff1);
  Vector shift(space, ff1);
  Scalar logh(space, ff1);

  fclose(ff1);

  // should the given shift be dropped?
  if (bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  bconfig.set_outputdir(outputdir);

  // last stage of the enabled stages
  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

  if (rank == 0) {
    std::cout << "Last Stage Enabled: " << last_stage << std::endl;
  }

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

  // the basis is of cartesian type in this case
  Base_tensor basis(space, CARTESIAN_BASIS);
  // the background metric is assumed to be flat, i.e. conformal flatness
  Metric_flat fmet(space, basis);

  // initialize fields defined by the coordinates
  CoordFields<Space_spheric_adapted> cf_generator(space);
  vec_ary_t coord_vectors{default_co_vector_ary(space)};

  // field representing the radius of the domain containing the surface
  scalar_ary_t coord_scalars;
  coord_scalars[R_BCO1] = Scalar(space);

  update_fields_co(cf_generator, coord_vectors, coord_scalars, xo);

  // initialize all fields correctly, based on the current domain setup
  update_fields_co(cf_generator, coord_vectors, coord_scalars, xo);

  // a level function, defining a root at a given fixed radius
  // helper construction to force the system to attain a fixed radius
  // instead resolving the correct surface
  Scalar level(space);
  level = (*coord_scalars[R_BCO1]) * (*coord_scalars[R_BCO1]) -
          bconfig(RMID) * bconfig(RMID);
  level.std_base();

  //end initialze coordinate fields

  if (rank == 0) {
    std::cout << "=================================" << endl;
    std::cout << "Single star solver 3D" << endl;
    std::cout << bconfig;
    std::cout << "=================================" << endl;
  }
  auto norot_solvers = [&](bool fixed) {
    if (fixed) {
      if (rank == 0)
        std::cout << "############################" << std::endl
                  << "TOV with a fixed radius" << std::endl
                  << "############################" << std::endl;
    } else {
      if (rank == 0)
        std::cout << "############################" << std::endl
                  << "TOV with a resolved surface" << std::endl
                  << "############################" << std::endl;
    }

    // setup a system of equations
    System_of_eqs syst(space, 0, ndom - 1);

    // call the (flat) conformal metric "f"
    fmet.set_system(syst, "f");

    // define numerical constants
    syst.add_cst("4piG", bconfig(BCO_QPIG));
    syst.add_cst("Mb", bconfig(MB));
    syst.add_cst("Madm", bconfig(MADM));

    // in case of a fixed radius solve the TOV with the given fixed central enthalpy
    // in case of a resolved surface, solve for the central enthalpy
    if (fixed) {
      syst.add_cst("Hc", loghc);
    } else {
      syst.add_var("Hc", loghc);
    }

    // include the coordinate fields
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);
    syst.add_cst("einf", *coord_vectors[S_INF]);
    syst.add_cst("lev", level);

    // the basic fields, conformal factor, lapse and (log) enthalpy
    syst.add_var("P", conf);
    syst.add_var("N", lapse);

    syst.add_var("H", logh);

    // define common combinations of conformal factor and lapse
    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    // define quantity to be integrated at infinity
    // two (in this case) equivalent definitions of ADM mass
    // as well as the Komar mass
    syst.add_def(ndom - 1, "intMadmalt = -dr(P) / 4piG * 2");
    syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
    syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");

    // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
    syst.add_def("h = exp(H)");

    // define the EOS operators
    Param p;
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);

    // define rest-mass density, internal energy and pressure through the enthalpy
    syst.add_def("rho = rho(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("press = press(h)");

    // definition to rescale the equations
    // delta = p / rho
    syst.add_def("delta = h - eps - 1.");

    for (int d = 0; d < ndom; d++) {
      switch (d) {
        // in the star the constraint equations are sourced by the matter
        case 0:
        case 1:
          // sources
          syst.add_def(d, "Etilde = press * h - press * delta");
          syst.add_def(d, "Stilde = 3 * press * delta");

          // constraint equations
          syst.add_def(d,
                       "eqP    = delta * D^i D_i P + 4piG / 2. * P^5 * Etilde");
          syst.add_def(d,
                       "eqNP   = delta * D^i D_i NP - 4piG / 2. * N * P^5 * "
                       "(Etilde + 2. * Stilde)");

          // definition for the baryonic mass integral
          syst.add_def(d, "intMb = P^6 * rho");
          // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
          syst.add_def(d, "firstint = H + log(N)");

          break;
        // outside the matter is absent and the sources are zero
        default:
          syst.add_eq_full(d, "H = 0");

          syst.add_def(d, "eqP = D^i D_i P");
          syst.add_def(d, "eqNP = D^i D_i NP");
          break;
      }
    }

    // add the constraint equations and demand continuity their normal derivative across domain boundaries
    space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");

    // boundary conditions at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");

    // if the radius of the stellar surface domain is fixed
    // use the helper construction, i.e. a level function with a root defining the radius
    if (fixed) {
      syst.add_eq_bc(1, OUTER_BC, "lev = 0");
    }
    // if the surface is resolved, define it to be where the matter vanishes
    else {
      syst.add_eq_bc(1, OUTER_BC, "H = 0");
    }

    // first integral in the innermost domains with non-zero matter content
    // and condition on the central value, either fixed directly or by the
    // integral below
    syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");

    // if surface is resolved, fix the central enthalpy by one of these integrals
    if (!fixed) {
      // fix the baryonic mass
      if (bconfig.control(MB_FIXING)) {
        space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
        if (rank == 0)
          std::cout << "Fixing Mb..." << std::endl;
      }
      // fix the ADM mass
      else {
        space.add_eq_int_inf(syst, "integ(intMadm) = Madm");
        if (rank == 0)
          std::cout << "Fixing Madm..." << std::endl;
      }
    }

    // print the variation of the surface radius over the whole star
    if (rank == 0) {
      auto rs = bco_utils::get_rmin_rmax(space, 1);
      std::cout << "[Rmin, Rmax] : [" << rs[0] << ", " << rs[1] << "]\n";
    }

    // parameters for the solver loop
    bool endloop = false;
    int ite = 1;
    double conv;

    // solve until convergence is achieved
    while (!endloop) {
      // do exactly one newton step, given the system above
      endloop = syst.do_newton(conv_thres, conv);

      // output files at this iteration and print diagnostics
      std::stringstream ss;
      ss << "norot_3d_";
      if (fixed) {
        ss << "fixed";
      } else {
        ss << "norot_bc";
      }
      ss << "_" << ite - 1;

      bconfig.set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics_norot(space, syst, bconfig, ite, conv);
        std::cout << std::endl;
        if (bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
      }

      // update all coordinate fields, in case the domain extents have changed
      update_fields_co(cf_generator, coord_vectors, {}, xo, &syst);
      ite++;
    }
    // update derived quantities in the config
    bconfig.set(HC) = std::exp(loghc);
    bconfig.set(NC) = EOS<eos_t, DENSITY>::get(bconfig(HC));

    // update the final ADM or baryonic mass depending on which one was fixed
    if (!fixed) {
      Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
      // we update the quasi-local measure of M_ADM with the true value
      // at infinity here
      bconfig.set(QLMADM) =
          space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
      if (bconfig.control(MB_FIXING))
        bconfig.set(MADM) =
            space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
      else
        bconfig.set(MB) = syst.give_val_def("intMb")()(0).integ_volume() +
                          syst.give_val_def("intMb")()(1).integ_volume();
    }
  };

  if (stage_enabled[PRE]) {
    norot_solvers(true);
    bconfig.set_filename(converged_filename("pre", bconfig));
    if (rank == 0) {
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
    }
    if (last_stage_idx != PRE)
      bconfig.set_stage(PRE) = false;
  }  // end pre stage

  if (stage_enabled[NOROT_BC]) {
    norot_solvers(false);
    bconfig.set_filename(converged_filename("norot_bc", bconfig));
    if (rank == 0) {
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
    }
    if (last_stage_idx != NOROT_BC && std::fabs(bconfig(CHI)) > 0)
      bconfig.set_stage(NOROT_BC) = false;
  }  // end norot stage

  // add this field only for the rotating cases
  // surface normal on the sphere outside of the star (not the deformed surface!)
  coord_vectors[S_BCO1] = Vector(space, COV, basis);
  update_fields_co(cf_generator, coord_vectors, {}, xo);

  if (stage_enabled[TOTAL_BC]) {
    if (rank == 0) {
      std::string line{"Rotating - with fixed MADM"};
      if (bconfig.control(MB_FIXING)) {
        line = "Rotating - with fixed Baryonic Mass";
      }
      std::cout << "############################" << std::endl
                << line << std::endl
                << "############################" << std::endl;
    }
    {
      // setup the system of equations
      System_of_eqs syst(space, 0, ndom - 1);

      // conformally flat background metric
      fmet.set_system(syst, "f");

      // constants
      syst.add_cst("4piG", bconfig(BCO_QPIG));

      // the constant coordinate vector fields
      syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
      syst.add_cst("ex", *coord_vectors[EX]);
      syst.add_cst("ey", *coord_vectors[EY]);
      syst.add_cst("ez", *coord_vectors[EZ]);
      syst.add_cst("einf", *coord_vectors[S_INF]);
      syst.add_cst("sm", *coord_vectors[S_BCO1]);

      // central enthalpy is a variable of the system in this case
      syst.add_var("Hc", loghc);

      // angular frequency parameter, coordinate dependent
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA))) {
        bconfig.set(OMEGA) = bconfig(FIXED_BCOMEGA);
        syst.add_cst("ome", bconfig(OMEGA));
        syst.add_var("chi", bconfig(CHI));
      } else {
        syst.add_cst("chi", bconfig(CHI));
        syst.add_var("ome", bconfig(OMEGA));
      }

      // Determine whether we fix the star by its baryonic mass
      // or ADM mass, in each case the other one has to be calculated alongside
      if (bconfig.control(MB_FIXING)) {
        syst.add_cst("Mb", bconfig(MB));
        syst.add_var("Madm", bconfig(MADM));
      } else {
        syst.add_var("Mb", bconfig(MB));
        syst.add_cst("Madm", bconfig(MADM));
      }

      // the basic fields: conformal factor, lapse and shift
      syst.add_var("P", conf);
      syst.add_var("N", lapse);
      syst.add_var("bet", shift);

      // matter, represented by the logarithm of the ethalpy
      syst.add_var("H", logh);

      // common definitions for the equations
      syst.add_def("NP = P*N");
      syst.add_def("Ntilde = N / P^6");

      // rotation vector field, including the shift
      syst.add_def("omega^i = bet^i + ome * mg^i");

      // the extrinsic curvature
      syst.add_def(
          "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / "
          "2. / Ntilde");

      // definition for integrals at infinity:
      // - ADM angular momentum
      // - two (in this case equivalent) definitions for the ADM mass
      // - Komar mass
      syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
      syst.add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
      syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P * 2 / 4piG");
      syst.add_def(ndom - 1,
                   "intMk = (einf^i * D_i N - A_ij * einf^i * bet^j) / 4piG");

      // define the EOS operators
      Param p;
      syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
      syst.add_ope("pres", &EOS<eos_t, PRESSURE>::action, &p);
      syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);

      for (int d = 0; d < ndom; d++) {
        switch (d) {
          case 0:
          case 1:
            // the enthalpy
            syst.add_def(d, "h = exp(H)");

            // density, internal energy and pressure as function of the enthalpy
            syst.add_def(d, "rho = rho(h)");
            syst.add_def(d, "eps = eps(h)");
            syst.add_def(d, "press = pres(h)");
            // p/rho to rescale the equations
            syst.add_def(d, "delta = h - eps - 1.");

            // definitions for the fluid 3-velocity
            // and its Lorentz factor
            syst.add_def(d, "U^i = omega^i / N");
            syst.add_def(d, "Usquare = P^4 * U_i * U^i");
            syst.add_def(d, "Wsquare = 1. / (1. - Usquare)");
            syst.add_def(d, "W = sqrt(Wsquare)");

            // rescaled sources and constraint equations
            syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
            syst.add_def(d,
                         "Stilde = 3 * press * delta + (Etilde + press * "
                         "delta) * Usquare");
            syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

            syst.add_def(d,
                         "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * "
                         "delta + 4piG / 2. * P^5 * Etilde");
            syst.add_def(d,
                         "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * "
                         "delta * A_ij *A^ij "
                         "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
            syst.add_def(
                d,
                "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
                "ptilde^i");

            // integrant of the baryonic mass integral
            syst.add_def(d, "intMb = P^6 * rho(h) * W");
            // the first integral of the Euler equation in case of a axisymmetric rotating star
            syst.add_def(d, "firstint = H + log(N) - log(W)");

            break;
          default:
            // outside the star the matter is absent and the sources are zero
            syst.add_eq_full(d, "H = 0");

            syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst.add_def(
                d, "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
            syst.add_def(d,
                         "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                         "A^ij * D_j Ntilde");
            break;
        }
      }

      // adding the constraint equations to the system
      // and demand continuity of the normal derivative across the domain boundaries
      space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
      space.add_eq(syst, "eqP= 0", "P", "dn(P)");
      space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

      // boundary conditions at infinity
      syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
      syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
      syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

      // boundary condition defining the surface and thus the adapted domain boundary
      syst.add_eq_bc(1, OUTER_BC, "H = 0");

      // first integral together with the condition at the center
      // Hc is fixed by the baryonic mass integral below
      syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");

      // baryonic mass integral of the domains containing matter
      space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");

      // ADM mass integral at infinity to calculate the ADM mass along the rest
      space.add_eq_int_inf(syst, "integ(intMadm) = Madm");

      // ADM angular momentum integral at infinity to fix the spin angular momentum of the star
      space.add_eq_int_inf(syst, "integ(intJ) - chi * Madm * Madm = 0");

      // solver stepping parameters
      bool endloop = false;
      int ite = 1;
      double conv;

      // repeat until convergence is achieved
      while (!endloop) {
        // do exactly one newton step
        endloop = syst.do_newton(conv_thres, conv);

        // output the data and diagnostics at this particular step
        std::stringstream ss;
        ss << "rot_3d_total_bc" << ite - 1;
        bconfig.set_filename(ss.str());

        if (rank == 0) {
          print_diagnostics_rot(space, syst, bconfig, ite, conv);
          if (bconfig.control(CHECKPOINT))
            bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
        }

        // recalculate all static coordinate vector fields,
        // in case the domain boundaries have changed
        update_fields_co(cf_generator, coord_vectors, {}, xo, &syst);
        // count the steps
        ite++;
      }
    }
  }  // TOTAL_BC

  if (stage_enabled[TESTING]) {
    Scalar phi(space);
    phi.annule_hard();
    phi.std_base();
    if (rank == 0) {
      std::string line{"BOOST TESTING"};
      if (bconfig.control(MB_FIXING)) {
        line = "Rotating - with fixed Baryonic Mass";
      }
      std::cout << "############################" << std::endl
                << line << std::endl
                << "############################" << std::endl;
    }
    {

      double H_scale1 = 0;
      // "background", unscaled logarithmic enthalpy from the previous step
      Scalar logh_const(logh);
      logh_const.std_base();

      // setup the system of equations
      System_of_eqs syst(space, 0, ndom - 1);

      // conformally flat background metric
      fmet.set_system(syst, "f");

      // constants
      syst.add_cst("4piG", bconfig(BCO_QPIG));

      // the constant coordinate vector fields
      syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
      syst.add_cst("ex", *coord_vectors[EX]);
      syst.add_cst("ey", *coord_vectors[EY]);
      syst.add_cst("ez", *coord_vectors[EZ]);
      syst.add_cst("einf", *coord_vectors[S_INF]);
      syst.add_cst("sm", *coord_vectors[S_BCO1]);

      // angular frequency parameter, coordinate dependent
      if (!std::isnan(bconfig.set(FIXED_BCOMEGA))) {
        bconfig.set(OMEGA) = bconfig(FIXED_BCOMEGA);
        syst.add_cst("ome", bconfig(OMEGA));
        syst.add_var("chi", bconfig(CHI));
      } else {
        syst.add_cst("chi", bconfig(CHI));
        syst.add_var("ome", bconfig(OMEGA));
      }
      syst.add_cst("yboost", bconfig(BVELY));

      syst.add_cst("Mb", bconfig(MB));
      syst.add_cst("Madm", bconfig(MADM));

      syst.add_var("Hscale1", H_scale1);
      syst.add_var("qlMadm", bconfig(QLMADM));

      // the basic fields: conformal factor, lapse and shift
      syst.add_var("P", conf);
      syst.add_var("N", lapse);
      syst.add_var("bet", shift);
      syst.add_var("phi", phi);

      // matter, represented by the logarithm of the ethalpy
      syst.add_cst("Hconst", logh_const);

      for (int d = 0; d < ndom; d++) {
        // the enthalpy is equal to the constant part everywhere
        // outside of the stars, i.e. zero
        if (d > 1)
          syst.add_def(d, "H  = Hconst");
        // inside the stars, it is the constant part
        // and a small correction factor to fix the correct
        // baryonic mass, due to slightly changing
        // fluid velocities
        else
          syst.add_def(d, "H  = Hconst * (1. + Hscale1)");
      }

      // common definitions for the equations
      syst.add_def("NP = P*N");
      syst.add_def("Ntilde = N / P^6");

      // rotation vector field, including the shift
      syst.add_def("omega^i = bet^i + ome * mg^i + yboost * ey^i");

      // the extrinsic curvature
      syst.add_def(
          "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / "
          "2. / Ntilde");

      // define the irrotational + spinning parts of the fluid velocity
      for (int d = 0; d <= 1; ++d) {
        syst.add_def(d, "s^i  = ome * mg^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }

      // definition for integrals at infinity:
      // - ADM angular momentum
      // - two (in this case equivalent) definitions for the ADM mass
      // - Komar mass
      syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
      syst.add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
      syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P * 2 / 4piG");
      syst.add_def(ndom - 1,
                   "intMk = (einf^i * D_i N - A_ij * einf^i * bet^j) / 4piG");

      // define the EOS operators
      Param p;
      syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
      syst.add_ope("pres", &EOS<eos_t, PRESSURE>::action, &p);
      syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
      syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

      for (int d = 0; d < ndom; d++) {
        switch (d) {
          case 0:
          case 1:
            // the enthalpy
            syst.add_def(d, "h = exp(H)");

            // density, internal energy and pressure as function of the enthalpy
            syst.add_def(d, "rho = rho(h)");
            syst.add_def(d, "eps = eps(h)");
            syst.add_def(d, "press = pres(h)");
            syst.add_def(d, "dHdlnrho = dHdlnrho(h)");
            // p/rho to rescale the equations
            syst.add_def(d, "delta = h - eps - 1.");

            // definitions for the fluid 3-velocity
            // and its Lorentz factor
            syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
            syst.add_def(d, "W      = sqrt(Wsquare)");
            syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
            syst.add_def(d, "Usquare= P^4 * U_i * U^i");
            syst.add_def(d, "V^i    = N * U^i - omega^i");
            syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
            syst.add_def(d,
                         "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 "
                         "* W * V^i)");

            // rescaled sources and constraint equations
            syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
            syst.add_def(d,
                         "Stilde = 3 * press * delta + (Etilde + press * "
                         "delta) * Usquare");
            syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

            syst.add_def(d,
                         "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * "
                         "delta + 4piG / 2. * P^5 * Etilde");
            syst.add_def(d,
                         "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * "
                         "delta * A_ij *A^ij "
                         "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
            syst.add_def(
                d,
                "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
                "ptilde^i");

            // integrant of the baryonic mass integral
            syst.add_def(d, "intMb = P^6 * rho(h) * W");
            // the first integral of the Euler equation in case of a axisymmetric rotating star
            //          syst.add_def(d, "firstint = H + log(N) - log(W)");

            break;
          default:
            // outside the star the matter is absent and the sources are zero
            syst.add_eq_full(d, "phi = 0");
            //syst.add_eq_full(d, "H = 0");

            syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst.add_def(
                d, "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
            syst.add_def(d,
                         "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                         "A^ij * D_j Ntilde");
            break;
        }
      }

      // adding the constraint equations to the system
      // and demand continuity of the normal derivative across the domain boundaries
      space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
      space.add_eq(syst, "eqP= 0", "P", "dn(P)");
      space.add_eq(syst, "eqbet^i=  0", "bet^i", "dn(bet^i)");

      // boundary conditions at infinity
      syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
      syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
      syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=yboost * ey^i");

      // boundary condition defining the surface and thus the adapted domain boundary
      syst.add_eq_bc(1, OUTER_BC, "H = 0");

      // baryonic mass integral of the domains containing matter
      space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");

      // ADM mass integral at infinity to calculate the ADM mass along the rest
      space.add_eq_int_inf(syst, "integ(intMadm) = qlMadm");

      // ADM angular momentum integral at infinity to fix the spin angular momentum of the star
      space.add_eq_int_inf(syst, "integ(intJ) - chi * Madm * Madm = 0");

      syst.add_eq_bc(1, OUTER_BC, "V^i * D_i H = 0");
      // in case of the stellar domains
      syst.add_eq_vel_pot(0, 2, "eqphi = 0", "phi=0");
      syst.add_eq_matching(0, OUTER_BC, "phi");
      syst.add_eq_matching(0, OUTER_BC, "dn(phi)");
      syst.add_eq_vel_pot(1, 2, "eqphi = 0", "phi=0");
      // solver stepping parameters
      bool endloop = false;
      int ite = 1;
      double conv;

      // repeat until convergence is achieved
      while (!endloop) {
        // do exactly one newton step
        endloop = syst.do_newton(conv_thres, conv);
        logh = syst.give_val_def("H");

        // output the data and diagnostics at this particular step
        std::stringstream ss;
        ss << "rot_3d_testing" << ite - 1;
        bconfig.set_filename(ss.str());

        if (rank == 0) {
          print_diagnostics_rot(space, syst, bconfig, ite, conv);
          if (bconfig.control(CHECKPOINT))
            bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh,
                                    phi);
        }

        // recalculate all static coordinate vector fields,
        // in case the domain boundaries have changed
        update_fields_co(cf_generator, coord_vectors, {}, xo, &syst);
        // count the steps
        ite++;
      }
    }
    // generate output filename
    bconfig.set_filename(converged_filename("TESTING", bconfig));

    // save data and configuration to files
    if (rank == 0) {
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
  }  // TESTING

  bconfig.set(QLMADM) = bconfig(MADM);

  // update derived quantities in the config
  bconfig.set(HC) = std::exp(loghc);
  bconfig.set(NC) = EOS<eos_t, DENSITY>::get(bconfig(HC));

  // generate output filename
  bconfig.set_filename(converged_filename("", bconfig));

  // save data and configuration to files
  if (rank == 0) {
    bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
  }

  // finalize MPI
  MPI_Finalize();
  return EXIT_SUCCESS;
}  // end ns solver 3d

// diagnostics at runtime
#define FORMAT std::setw(10) << std::left << std::showpos

template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics_norot(space_t const& space,
                             syst_t const& syst,
                             config_t& bconfig,
                             int const ite,
                             double const conv) {

  // total number of domains
  int ndom = space.get_nbr_domains();

  // compute the baryonic mass at volume integral from the given integrant
  double baryonic_mass = syst.give_val_def("intMb")()(0).integ_volume() +
                         syst.give_val_def("intMb")()(1).integ_volume();

  // compute the ADM mass as surface integral at infinity
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // get the maximum and minimum coordinate radius along the surface,
  // i.e. the adapted domain boundary
  auto rs = bco_utils::get_rmin_rmax(space, 1);

  // alternative, equivalent ADM mass integral
  Val_domain integMadmalt(syst.give_val_def("intMadmalt")()(ndom - 1));
  double Madmalt = space.get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

  // output to standard output
  std::ios_base::fmtflags f(std::cout.flags());
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl
            << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Madm_ql: " << Madmalt << " ["
            << std::abs(Madm - Madmalt) / Madm << "]" << std::endl
            << FORMAT << "Mk: " << Mk << " [" << std::abs(Madm - Mk) / Madm
            << "]" << std::endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n";
  std::cout.flags(f);
}  // end print diagnostics norot

// runtime diagnostics specific for rotating solutions
template <typename space_t, typename syst_t, typename config_t>
void print_diagnostics_rot(space_t const& space,
                           syst_t const& syst,
                           config_t& bconfig,
                           int const ite,
                           double const conv) {

  // print all the diagnostics as in the non-rotating case first
  print_diagnostics_norot(space, syst, bconfig, ite, conv);

  // total number of domains
  int ndom = space.get_nbr_domains();

  // compute the ADM angular momentum as surface integral at infinity
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

  // compute the ADM mass as surface integral at infinity
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // output the dimensionless spin and angular frequency parameter
  std::ios_base::fmtflags f(std::cout.flags());
  std::cout << FORMAT << "Chi: " << J / Madm / Madm << std::endl
            << FORMAT << "Omega: " << bconfig(OMEGA) << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << "\n\n";
}  // end print diagnostics rot

//standardized filename for each converged dataset at the end of each stage.
template <typename config_t>
std::string converged_filename(const std::string& stage, config_t bconfig) {
  std::stringstream ss;
  ss << "converged_NS";
  if (stage != "")
    ss << "_" << stage << ".";
  else
    ss << ".";
  ss << bconfig(MADM) << "." << bconfig(CHI) << "." << std::setfill('0')
     << std::setw(2) << bconfig(BCO_RES);
  return ss.str();
}
