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

#include "EOS/EOS.hh"
#include "EOS/FUKA_EOS_Utilities.hh"
#include "bco_utilities.hpp"
#include "kadath_bin_ns.hpp"

// config_file includes
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"

#include <iostream>
#include <numeric>

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;
using namespace Kadath::FUKA_EOS;

// conversion from solar mass to km
double M2km = 1.4769994423016508;

// convenient type definitions
using vec_d = std::vector<double>;
using ary_d = std::array<double, 2>;
using ary_i = std::array<int, 2>;

// forward declaration
template <class eos_t>
struct bns_reader {
  template <class config_t>
  void operator()(config_t& bconfig) {
    if (std::isnan(bconfig.set(MADM, BCO1)) ||
        std::isnan(bconfig.set(MADM, BCO2))) {
      std::cerr << "Missing \"Madm\" in config file\n"
                   "Setting to quasi-local \"Madm\"! \n";
      bconfig.set(MADM, BCO1) = bconfig(QLMADM, BCO1);
      bconfig.set(MADM, BCO2) = bconfig(QLMADM, BCO2);
    }
    if (std::isnan(bconfig.set(COMY)))
      bconfig.set(COMY) = 0.;

    // read domain decomposition and fields from binary file
    std::string in_spacefile = bconfig.space_filename();
    FILE* fich = fopen(in_spacefile.c_str(), "r");

    Space_bin_ns space(fich);
    Scalar conf(space, fich);
    Scalar lapse(space, fich);
    Vector shift(space, fich);
    Scalar logh(space, fich);
    Scalar phi(space, fich);
    fclose(fich);

    // define basic fields
    Base_tensor basis(shift.get_basis());
    Metric_flat fmet(space, basis);

    // coordinate dependent fields
    CoordFields<Space_bin_ns> cfields(space);

    // get center of each star
    const std::array<double, 2> x_nuc{bco_utils::get_center(space, space.NS1),
                                      bco_utils::get_center(space, space.NS2)};

    // coordinate origin
    double xo = 0;

    // coordinate dependent vector fields
    vec_ary_t coord_vectors = default_binary_vector_ary(space);

    // vector field for expansion factor
    Vector CART(space, CON, basis);
    CART = cfields.cart();

    // initialize all coordinate dependent fields with the given centers
    update_fields(cfields, coord_vectors, {}, xo, x_nuc[0], x_nuc[1]);

    // compute the location of the maximum densities
    std::array<double, 2> x_max;
    {
      System_of_eqs syst_H(space);
      fmet.set_system(syst_H, "f");

      syst_H.add_cst("H", logh);
      syst_H.add_cst("ex", *coord_vectors[EX]);

      syst_H.add_def("dH = (ex^i * D_i H) / H");
      syst_H.add_def("dH2 = ex^i * D_i dH");

      Scalar dHdx = syst_H.give_val_def("dH");
      Scalar dHdx2 = syst_H.give_val_def("dH2");

      for (int j : {0, 1}) {
        double xc = x_nuc[j];
        double err = 1;

        Point absol(3);
        absol.set(1) = xc;
        absol.set(2) = 0;
        absol.set(3) = 0;

        while (err > 1e-14) {
          xc = xc - dHdx.val_point(absol) / dHdx2.val_point(absol);
          absol.set(1) = xc;
          err = std::abs(dHdx.val_point(absol));
        }
        x_max[j] = xc;
      }
    }

    // compute the "center of mass" shifted maximum density locations
    double x_max_com1 = x_max[0] + bconfig(COM);
    double x_max_com2 = x_max[1] + bconfig(COM);

    // setup a system of equations to compute derived quantities
    // using the spectral representation of the fields
    int ndom = space.get_nbr_domains();
    System_of_eqs syst(space);

    // conformally flat background metric
    fmet.set_system(syst, "f");

    // setup the eos operators
    Param p;
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
    syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

    // constants in the system
    syst.add_cst("4piG", bconfig(QPIG));

    syst.add_cst("Madm1", bconfig(MADM, BCO1));
    syst.add_cst("Mb1", bconfig(MB, BCO1));
    syst.add_cst("chi1", bconfig(CHI, BCO1));

    syst.add_cst("Madm2", bconfig(MADM, BCO2));
    syst.add_cst("Mb2", bconfig(MB, BCO2));
    syst.add_cst("chi2", bconfig(CHI, BCO2));

    // coordinate fields
    // be aware that these names are hardcoded in the coord_fields.hpp!
    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);

    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // "center of mass" and orbital angular frequency parameter
    syst.add_cst("xaxis", bconfig(COM));
    syst.add_cst("yaxis", bconfig(COMY));
    syst.add_cst("ome", bconfig(GOMEGA));

    // the actual fields representing the solution
    syst.add_cst("P", conf);
    syst.add_cst("N", lapse);
    syst.add_cst("bet", shift);
    syst.add_cst("H", logh);

    // different definitions for the velocity
    // between the corotation and the irrotational / spinning cases
    if (bconfig.control(COROT_BIN)) {
      syst.add_cst("omes1", bconfig(OMEGA, BCO1));
      syst.add_cst("omes2", bconfig(OMEGA, BCO2));
    } else {
      syst.add_cst("phi", phi);
      syst.add_cst("omes1", bconfig(OMEGA, BCO1));
      syst.add_cst("omes2", bconfig(OMEGA, BCO2));

      for (int d = space.NS1; d <= space.ADAPTED1; ++d) {
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      for (int d = space.NS2; d <= space.ADAPTED2; ++d) {
        syst.add_def(d, "s^i  = omes2 * mp^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
    }

    // common definitions
    syst.add_def("NP     = P*N");
    syst.add_def("Ntilde = N / P^6");

    // definitions of derived matter quantities
    syst.add_def("h = exp(H)");
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("dHdlnrho = dHdlnrho(h)");
    syst.add_def("delta = h - eps - 1.");

    // orbital vector field and resulting "total" shift
    // check for ADOT so we don't get errors.
    std::string eccstr{};
    if (!std::isnan(bconfig.set(ADOT))) {
      syst.add_cst("adot", bconfig(ADOT));
      syst.add_cst("r", CART);
      syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
      eccstr += " + adot * comr^i";
    }
    std::string omegastr{"omega^i = bet^i + ome * Morb^i" + eccstr};
    syst.add_def("Morb^i = mg^i + xaxis * ey^i");
    syst.add_def(omegastr.c_str());

    // normalized derivative of the enthalpy along binary axis
    syst.add_def("dH = (ex^i * D_i H) / H");

    // conformal extrinsic curvature
    syst.add_def(
        "A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    // integrant of the ADM linear momentum at infinity
    syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i");
    syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i");
    syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i");

    // definitions from https://arxiv.org/abs/1506.01689
    syst.add_def(ndom - 1, "intCOMx = 3. / 2. / 4piG * P^4 * ex_i * einf^i");
    syst.add_def(ndom - 1, "intCOMy = 3. / 2. / 4piG * P^4 * ey_i * einf^i");
    syst.add_def(ndom - 1, "intCOMz = 3. / 2. / 4piG * P^4 * ez_i * einf^i");

    // quasi-local spin angular momentum surface integral integrants for both stars
    syst.add_def(space.ADAPTED1 + 1, "intS = A_ij * mm^i * sm^j / 2. / 4piG");
    syst.add_def(space.ADAPTED2 + 1, "intS = A_ij * mp^i * sp^j / 2. / 4piG");

    for (int d = 0; d < ndom; d++) {
      if ((d >= space.ADAPTED2 + 1) || d == space.ADAPTED1 + 1) {
      } else {
        if (bconfig.control(COROT_BIN)) {
          syst.add_def(d, "U^i    = omega^i / N");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "firstint = log(h * N / W)");
        } else {
          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "V^i    = N * U^i - omega^i");
          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
        }

        syst.add_def(d, "intMb  = P^6 * rho * W");
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");

        syst.add_def(d, "intH  = P^6 * H * W");
      }
    }

    // ADM and Komar integrants at infinity
    syst.add_def(ndom - 1, "intJ = multr(A_ij * Morb^j * einf^i) / 2. / 4piG");
    syst.add_def(ndom - 1, "Madm = -dr(P) * 2 / 4piG");
    syst.add_def(ndom - 1, "Mk   =  dr(N) / 4piG");
    // integrant to compute proper area over a spherical surface
    syst.add_def("intMsq = P^4");

    // BNS component quantities
    ary_i nuc_doms{space.NS1, space.NS2};
    ary_i adapted_doms{space.ADAPTED1, space.ADAPTED2};
    ary_i shells{};
    ary_d areal_r;
    ary_d central_logh;
    ary_d central_rho;
    ary_d central_P;
    ary_d central_dHdx;
    ary_d central_euler;
    ary_d H_int;
    ary_d ql_madm;
    ary_d ql_spin;
    ary_d madm;
    ary_d mbs;
    ary_d xcom;

    std::vector<ary_d> r_extrema;
    std::vector<vec_d> mb_distro;

    // compute everything for both stars
    for (int i : {0, 1}) {
      int dom = adapted_doms[i];
      shells[i] = dom - nuc_doms[i] - 1;

      // compute area and areal radius from that
      double A = space.get_domain(dom + 1)->integ(syst.give_val_def("intMsq")()(
                                                      dom + 1),
                                                  INNER_BC);
      areal_r[i] = sqrt(A / 4. / acos(-1.));

      // point defining the center of the star
      Point ns_c(3);
      ns_c.set(1) = x_nuc[i];
      ns_c.set(2) = 0;
      ns_c.set(3) = 0;

      // central matter quantities
      central_logh[i] = logh.val_point(ns_c);
      central_rho[i] = EOS<eos_t, DENSITY>::get(std::exp(central_logh[i]));
      central_P[i] = EOS<eos_t, PRESSURE>::get(std::exp(central_logh[i]));
      central_dHdx[i] = syst.give_val_def("dH")().val_point(ns_c);
      central_euler[i] = syst.give_val_def("firstint")().val_point(ns_c);

      // volume integrated quantities, log enthalpy, quasi-local ADM mass and baryonic mass
      H_int[i] = 0;
      vec_d ql_m{};
      vec_d baryonic_mass{};

      for (int d = nuc_doms[i]; d <= adapted_doms[i]; ++d) {
        H_int[i] += syst.give_val_def("intH")()(d).integ_volume();
        ql_m.push_back(syst.give_val_def("intM")()(d).integ_volume());
        baryonic_mass.push_back(syst.give_val_def("intMb")()(d).integ_volume());
      }

      mb_distro.push_back(baryonic_mass);
      mbs[i] = std::accumulate(baryonic_mass.begin(), baryonic_mass.end(), 0.);

      ql_madm[i] = std::accumulate(ql_m.begin(), ql_m.end(), 0.);
      // get quasi-local spin as surface integral outside of the star
      ql_spin[i] =
          space.get_domain(dom + 1)->integ(syst.give_val_def("intS")()(dom + 1),
                                           OUTER_BC);

      // the correct ADM mass (at infinite separation), is given from the single star solution
      // and thus fixed
      madm[i] = bconfig(MADM, i);

      // minimal and maximal (coordinate) radius across the adapted surface
      r_extrema.push_back(bco_utils::get_rmin_rmax(space, dom));

      // center-of-mass shifted center
      xcom[i] = x_nuc[i] + bconfig(COM);
    }

    // binary quantities
    double adm_inf =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("Madm")()(ndom - 1),
                                          OUTER_BC);
    double komar =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("Mk")()(ndom - 1),
                                          OUTER_BC);
    double Jinf =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("intJ")()(ndom - 1),
                                          OUTER_BC);
    double Px = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPx")()(
                                                      ndom - 1),
                                                  OUTER_BC);
    double Py = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPy")()(
                                                      ndom - 1),
                                                  OUTER_BC);
    double Pz = space.get_domain(ndom - 1)->integ(syst.give_val_def("intPz")()(
                                                      ndom - 1),
                                                  OUTER_BC);

    // center of mass defined like in https://arxiv.org/abs/1506.01689
    double COMx = space.get_domain(ndom - 1)->integ(syst.give_val_def(
                                                        "intCOMx")()(ndom - 1),
                                                    OUTER_BC) /
                  adm_inf;
    double COMy = space.get_domain(ndom - 1)->integ(syst.give_val_def(
                                                        "intCOMy")()(ndom - 1),
                                                    OUTER_BC) /
                  adm_inf;
    double COMz = space.get_domain(ndom - 1)->integ(syst.give_val_def(
                                                        "intCOMz")()(ndom - 1),
                                                    OUTER_BC) /
                  adm_inf;

    // ADM mass at infinite separation
    double Minf = std::accumulate(madm.begin(), madm.end(), 0.);
    // binding energy
    double e_bind = adm_inf - Minf;
    // end binary quantities

#ifdef FORMAT
#undef FORMAT
#endif

#define FORMAT1                                                     \
  std::setw(25) << std::right << std::setprecision(5) << std::fixed \
                << std::showpos
#define FORMAT                                                           \
  std::setw(25) << std::right << std::setprecision(5) << std::scientific \
                << std::showpos

    // helper lambda to print shell radii
    auto print_shells = [&](int dom_min, int dom_max) {
      int cnt = 1;
      for (int i = dom_min; i < dom_max; ++i) {
        std::string shell{"SHELL" + std::to_string(cnt) + " = "};
        std::cout << FORMAT1 << shell << "["
                  << bco_utils::get_radius(space.get_domain(i), INNER_BC)
                  << ", " << bco_utils::get_radius(space.get_domain(i), EQUI)
                  << "]" << std::endl;
        cnt++;
      }
    };

    // helper to print contents of a collection of values
    auto print_vec = [&](auto& vec) {
      for (auto& e : vec) {
        std::cout << e << ",";
      }
    };

    auto res_r = space.get_domain(0)->get_nbr_points()(0);
    auto res_t = space.get_domain(0)->get_nbr_points()(1);
    auto res_p = space.get_domain(0)->get_nbr_points()(2);

    // string to identify the stars
    std::array<std::string, 2> ns_str{" NS_MINUS ", " NS_PLUS "};
    ary_i bounds{space.NS2 - 1, space.OUTER - 1};
    std::string header(22, '#');
    for (int i = 0; i <= 1; ++i) {
      std::cout << header + ns_str[i] + header + "\n";
      std::cout << FORMAT1 << "Center_COM = " << "(" << xcom[i] << ", 0, 0)\n"
                << FORMAT1 << "Coord R_IN = "
                << bco_utils::get_radius(space.get_domain(nuc_doms[i]), EQUI)
                << '\n'
                << FORMAT1 << "Coord R = " << "[" << r_extrema[i][0] << ","
                << r_extrema[i][1] << "] (" << "[" << r_extrema[i][0] * M2km
                << "," << r_extrema[i][1] * M2km << "] km)\n";
      print_shells(adapted_doms[i] + 1, bounds[i]);
      std::cout
          << FORMAT1 << "Coord R_OUT = "
          << bco_utils::get_radius(space.get_domain(adapted_doms[i] + 1), EQUI)
          << "\n"
          << FORMAT1 << "Areal R = " << areal_r[i] << " [" << areal_r[i] * M2km
          << "km]\n"
          // baryonic mass and fractions per stellar domain covering the star
          << FORMAT1 << "Baryonic Mass = " << mbs[i] << " (";
      print_vec(mb_distro[i]);
      std::cout
          << ")\n"
          << FORMAT1 << "Isolated ADM Mass = " << bconfig(MADM, i) << "\n"
          << FORMAT1 << "Quasi-local Madm = " << ql_madm[i]
          << " Diff:" << std::fabs(1. - ql_madm[i] / bconfig(MADM, i))
          << std::endl
          // quasi-local spin angular momentum
          << FORMAT1 << "Quasi-local S = " << ql_spin[i]
          << std::endl
          // dimensionless spin (constant, given by the imported single star!)
          // angular frequency paramter of the star, describing the magnitude of the spin component of the velocity field
          << FORMAT1
          << "Chi = " << ql_spin[i] / bconfig(MADM, i) / bconfig(MADM, i)
          << " [" << bconfig(CHI, i) << "]\n"
          << FORMAT1 << "Omega = " << bconfig(OMEGA, i)
          << std::endl
          // location of the maximal density (and compared to the nucleus' center)
          << FORMAT1 << "x(max(Density)) = " << x_max[i] << " ("
          << (x_nuc[i] - x_max[i]) / x_nuc[i] << ")\n"
          << FORMAT << "Central Density = " << central_rho[i] << std::endl
          << FORMAT << "Central log(h) = " << central_logh[i] << std::endl
          << FORMAT << "Central Pressure = " << central_P[i] << std::endl
          << FORMAT << "Central dlog(h)/dx = " << central_dHdx[i]
          << std::endl
          // central values of the Euler constant, log enthalpy and its derivative
          << FORMAT1 << "Central Euler Constant = " << central_euler[i]
          << std::endl
          // integrated log enthalpy
          << FORMAT1 << "Integrated log(h) = " << H_int[i] << "\n\n";
    }

    // boundaries of additional outer shells
    auto outer_shells = space.get_n_shells_outer();
    if (outer_shells > 0) {
      std::cout << FORMAT1 << "Outer shell bounds\n";
      print_shells(ndom - 1 - outer_shells, ndom - 1);
    }
    auto M1 = bconfig(MADM, BCO1);
    auto M2 = bconfig(MADM, BCO2);
    if (M2 > M1)
      std::swap(M1, M2);
    // mass ratio, ratio of the ADM masses at infinity
    std::cout
        << header + " Binary " + header + '\n'
        << FORMAT1 << std::fixed << "RES = " << "[" << res_r << "," << res_t
        << "," << res_p << "]\n"
        << FORMAT1 << "Q = " << M2 / M1
        << std::endl
        // distance between the stellar centers, defined at the setup stage
        << FORMAT1 << std::setprecision(2) << "Separation = " << bconfig(DIST)
        << " [" << bconfig(DIST) / Minf << "] (" << bconfig(DIST) * M2km
        << "km)"
        << std::endl
        // orbital angular frequency parameter
        << FORMAT1 << "Orbital Omega = " << bconfig(GOMEGA)
        << std::endl
        // Komar and ADM mass of the binary
        << FORMAT1 << "Komar mass = " << komar << std::endl
        << FORMAT1 << "Adm mass = " << adm_inf
        << ", Diff: " << fabs(adm_inf - komar) / (adm_inf + komar) << std::endl
        << FORMAT1 << "Total mass (Minf) = " << Minf
        << std::endl
        // ADM angular momentum of the binary
        << FORMAT1 << "Adm moment. = " << Jinf
        << std::endl
        // binding energy, defined by the gravitational mass difference at finite separation
        << FORMAT << "Binding energy = " << e_bind
        << std::endl
        // dimensionless orbital frequency
        << FORMAT << "Minf * Ome = " << Minf * bconfig(GOMEGA)
        << std::endl
        // dimensionless binding energy
        << FORMAT << "E_b / Minf = " << e_bind / Minf
        << std::endl
        // ADM linear momentum
        << FORMAT << "Px = " << Px << std::endl
        << FORMAT << "Py = " << Py << std::endl
        << FORMAT << "Pz = " << Pz
        << std::endl
        // "center of mass" defined by a vanishing ADM momentum at infinity
        // With an analytical estimate of the center of mass estimate from Osokine+
        << FORMAT1 << "COMx = " << bconfig(COM) << ", A-COMx = " << COMx
        << std::endl
        << FORMAT1 << "COMy = " << bconfig(COMY) << ", A-COMy = " << COMy
        << std::endl
        << FORMAT1 << "A-COMz = " << COMz << std::endl;
  }
};

int main(int argc, char** argv) {
  // usage information
  if (argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<str: BNS ID basename>.info "
              << std::endl;
    std::cerr << "e.g. ./reader converged.TOTAL_BC.9.info" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  // read the configuration
  std::string in_filename = argv[1];
  kadath_config_boost<BIN_INFO> bconfig(in_filename);

  // setup the eos
  EOS_initialize::init(bconfig, NODES::BCO1);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);
  EOS_Function_Dispatcher::dispatch<bns_reader>(bconfig, eos_type, bconfig);

  return EXIT_SUCCESS;
}