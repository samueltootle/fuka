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
 */
#include "EOS/EOS.hh"
#include "EOS/FUKA_EOS_Utilities.hh"
#include "bco_utilities.hpp"
#include "kadath.hpp"

// config_file includes
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"

#include <iostream>
#include <numeric>
using namespace Kadath;
using namespace Kadath::FUKA_Config;
using namespace Kadath::FUKA_EOS;

double M2km = 1.4769994423016508;

template <class eos_t>
struct reader_output {
  template <class config_t>
  void operator()(config_t bconfig) {
    std::string in_spacefile = bconfig.space_filename();

    FILE* fich = fopen(in_spacefile.c_str(), "r");
    Space_bhns space(fich);
    Scalar conf(space, fich);
    Scalar lapse(space, fich);
    Vector shift(space, fich);
    Scalar logh(space, fich);
    Scalar phi(space, fich);
    fclose(fich);

    Base_tensor basis(shift.get_basis());
    Metric_flat fmet(space, basis);
    CoordFields<Space_bhns> cfields(space);
    int ndom = space.get_nbr_domains();

    double xc1 = bco_utils::get_center(space, space.NS);
    double xc2 = bco_utils::get_center(space, space.BH);
    double xo = bco_utils::get_center(space, ndom - 1);

    // setup coord fields
    vec_ary_t coord_vectors = default_binary_vector_ary(space);
    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);
    // end setup coord fields

    System_of_eqs syst(space);
    fmet.set_system(syst, "f");

    // setup operators for the EOS
    Param p;
    set_eos_ope_struct<eos_t>()(syst, p);

    syst.add_cst("4piG", bconfig(QPIG));
    syst.add_cst("PI", M_PI);

    syst.add_cst("Mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("m1", *coord_vectors[BCO1_ROT]);
    syst.add_cst("m2", *coord_vectors[BCO2_ROT]);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);

    syst.add_cst("s1", *coord_vectors[S_BCO1]);
    syst.add_cst("s2", *coord_vectors[S_BCO2]);
    syst.add_cst("esurf", *coord_vectors[S_INF]);

    syst.add_cst("xaxis", bconfig(COM));
    syst.add_cst("yaxis", bconfig(COMY));
    syst.add_cst("ome", bconfig(GOMEGA));
    syst.add_cst("omes1", bconfig(OMEGA, BCO1));

    for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
      syst.add_def(d, "s^i  = omes1 * m1^i");
    }

    syst.add_cst("P", conf);
    syst.add_cst("N", lapse);
    syst.add_cst("bet", shift);
    syst.add_cst("H", logh);
    syst.add_cst("phi", phi);

    syst.add_def("NP     = P*N");
    syst.add_def("Ntilde = N / P^6");

    syst.add_def("Morb^i = Mg^i + xaxis * ey^i + yaxis * ex^i");
    syst.add_def("B^i = bet^i + ome * Morb^i");

    syst.add_def("dH = (ex^i * D_i H) / H");

    syst.add_def(
        "A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
        "Ntilde");

    syst.add_def("intPx = A_i^j * ex_j * esurf^i / 8 / PI");
    syst.add_def("intPy = A_i^j * ey_j * esurf^i / 8 / PI");
    syst.add_def("intPz = A_i^j * ez_j * esurf^i / 8 / PI");
    syst.add_def(space.BH + 2, "h^ij = f^ij - s2^i * s2^j");

    syst.add_def(space.ADAPTEDNS + 1, "intS1 = A_ij * m1^i * s1^j / 8. / PI");
    syst.add_def(space.ADAPTEDBH + 1, "intS2 = A_ij * m2^i * s2^j / 8. / PI");
    syst.add_def("intAsq = P^4 / 4 / PI");
    syst.add_def("intMsq = intAsq / 4");
    syst.add_def(space.ADAPTEDBH + 1, "intMsq = P^4 / 16. / PI");
    syst.add_def(space.ADAPTEDNS + 1, "intMsq = P^4 / 4.  / PI");

    syst.add_def(ndom - 1, "COMx  = -3 * P^4 *esurf^i * ex_i / 8. / PI");
    syst.add_def(ndom - 1, "COMy  = 3 * P^4 * esurf^i * ey_i / 8. / PI");
    syst.add_def(ndom - 1, "COMz  = 3 * P^4 * esurf^i * ez_i / 8. / PI");
    syst.add_def("dtgamma = D_k bet^k + 6 / P * B^k * D_k P");

    for (int d = 0; d < ndom; d++) {
      if (d >= space.ADAPTEDNS + 1) {
        syst.add_def(d, "eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst.add_def(d,
                     "eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,
                     "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij "
                     "* D_j Ntilde");

      } else {
        syst.add_def(d, "h      = exp(H)");

        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");

        syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
        syst.add_def(d, "W      = sqrt(Wsquare)");

        syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
        syst.add_def(d, "V^i    = N * U^i - B^i");

        syst.add_def(d, "Usquare= P^4 * U_i * U^i");

        syst.add_def(d, "E      = rho(h) * h * Wsquare - press(h)");
        syst.add_def(d, "S      = 3 * press(h) + (E + press(h)) * Usquare");
        syst.add_def(d, "p^i    = rho(h) * h * Wsquare * U^i");

        syst.add_def(d, "intMb  = P^6 * rho(h) * W");
        syst.add_def(d, "intM   = - D_i D^i P / 2. / PI");

        syst.add_def(d, "firstint = h * N / W + eta_i * V^i");

        syst.add_def(d, "eqphi  = D_i (P^6 * W * rho(h) * V^i)");

        syst.add_def(
            d,
            "eqP    = D^i D_i P + A_ij * A^ij / P^7 / 8 + 4piG / 2. * P^5 * E");
        syst.add_def(d,
                     "eqNP   = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij *A^ij "
                     "- 4piG / 2. * N * P^5 * (E + 2. * S)");
        syst.add_def(d,
                     "eqbet^i= D_j D^j bet^i + D^i D_j bet^j / 3. "
                     "- 2. * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * p^i");
        syst.add_def(d, "intH  = P^6 * H * W");
      }
    }
    syst.add_def(ndom - 1, "intJ = multr(A_ij * Morb^j * esurf^i) / 8 / PI");
    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");
    syst.add_def(ndom - 1, "Mk   =  dr(N) / 4 / PI");

    // NS Quantities
    double loghc1 = bco_utils::get_boundary_val(space.NS, logh, INNER_BC);
    double pressc1 = EOS<eos_t, PRESSURE>::get(std::exp(loghc1));
    double rhoc1 = EOS<eos_t, DENSITY>::get(std::exp(loghc1));

    std::vector<double> baryonic_mass1{};
    std::vector<double> int_H1{};
    std::vector<double> adm_mass1{};
    for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
      baryonic_mass1.push_back(syst.give_val_def("intMb")()(d).integ_volume());
      int_H1.push_back(syst.give_val_def("intH")()(d).integ_volume());
      adm_mass1.push_back(syst.give_val_def("intM")()(d).integ_volume());
    }
    double MB1 =
        std::accumulate(baryonic_mass1.begin(), baryonic_mass1.end(), 0.);
    double ql_madm1 = std::accumulate(adm_mass1.begin(), adm_mass1.end(), 0.);
    double intH1 = std::accumulate(int_H1.begin(), int_H1.end(), 0.);
    double ql_spinns =
        space.get_domain(space.ADAPTEDNS + 1)
            ->integ(syst.give_val_def("intS1")()(space.ADAPTEDNS + 1),
                    OUTER_BC);

    auto [rmin, rmax] = bco_utils::get_rmin_rmax(space, space.ADAPTEDNS);
    double A = space.get_domain(space.ADAPTEDNS + 1)
                   ->integ(syst.give_val_def("intAsq")()(space.ADAPTEDNS + 1),
                           INNER_BC);
    double areal_rns = sqrt(A);

    double rin_ns = bco_utils::get_radius(space.get_domain(space.NS), EQUI);
    double rout_ns = bco_utils::get_radius(
        space.get_domain(space.ADAPTEDNS + 1 + bconfig(NSHELLS, BCO1)), EQUI);

    auto dHdx = syst.give_val_def("dH")();
    double dHdx1 = bco_utils::get_boundary_val(space.NS, dHdx, INNER_BC);
    double euler1 =
        bco_utils::get_boundary_val(space.NS, syst.give_val_def("firstint")(),
                                    INNER_BC);
    double NS_x_com = xc1 + bconfig(COM);
    double NS_py =
        space.get_domain(space.ADAPTEDNS + 1)
            ->integ(syst.give_val_def("intPy")()(space.ADAPTEDNS + 1),
                    OUTER_BC);
    double NS_px =
        space.get_domain(space.ADAPTEDNS + 1)
            ->integ(syst.give_val_def("intPx")()(space.ADAPTEDNS + 1),
                    OUTER_BC);

    auto npts = space.get_domain(space.ADAPTEDNS)->get_nbr_points();
    Index pos_eq(npts);
    pos_eq.set(0) = npts(0) - 1;  /// Set to outer radius
    pos_eq.set(1) = npts(1) - 1;  /// Set theta to be on the xy plane.

    Index pos_pole(npts);
    pos_pole.set(0) = npts(0) - 1;  /// Set to outer radius

    Val_domain logh_dr(logh(space.ADAPTEDNS).der_r());
    double mass_shedding_parameter = logh_dr(pos_eq) / logh_dr(pos_pole);

    double conf_eq = conf(space.ADAPTEDNS)(pos_eq);
    double Circumferential_R =
        conf_eq * conf_eq * space.get_domain(space.ADAPTEDNS)->get_radius()(pos_eq);
    double& CR = Circumferential_R;
    // END NS Quantities

    // BH Quantities
    double ql_spinbh =
        space.get_domain(space.ADAPTEDBH + 1)
            ->integ(syst.give_val_def("intS2")()(space.ADAPTEDBH + 1),
                    OUTER_BC);

    double mirrsq =
        space.get_domain(space.BH + 2)
            ->integ(syst.give_val_def("intMsq")()(space.BH + 2), INNER_BC);
    double mirr = sqrt(mirrsq);
    double mch = sqrt(mirrsq + ql_spinbh * ql_spinbh / 4. / mirrsq);
    double rin_bh = bco_utils::get_radius(space.get_domain(space.BH), EQUI);
    double r_bh = bco_utils::get_radius(space.get_domain(space.BH + 1), EQUI);
    double rout_bh =
        bco_utils::get_radius(space.get_domain(space.OUTER - 1), EQUI);
    A = space.get_domain(space.ADAPTEDBH + 1)
            ->integ(syst.give_val_def("intAsq")()(space.ADAPTEDBH + 1),
                    INNER_BC);
    double areal_rbh = sqrt(A);
    auto [lapsemin, lapsemax] =
        bco_utils::get_field_min_max(lapse, space.BH + 2, INNER_BC);
    auto [confmin, confmax] =
        bco_utils::get_field_min_max(conf, space.BH + 2, INNER_BC);
    double BH_x_com = xc2 + bconfig(COM);
    double BH_py =
        space.get_domain(space.ADAPTEDBH + 1)
            ->integ(syst.give_val_def("intPy")()(space.ADAPTEDBH + 1),
                    INNER_BC);
    double BH_px =
        space.get_domain(space.ADAPTEDBH + 1)
            ->integ(syst.give_val_def("intPx")()(space.ADAPTEDBH + 1),
                    INNER_BC);
    // END BH Quantities

    // binary quantities
    double adm_inf =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("Madm")()(ndom - 1),
                                          OUTER_BC);
    double komar =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("Mk")()(ndom - 1),
                                          OUTER_BC);
    double e_diff = fabs(2 * (adm_inf - komar) / (adm_inf + komar));
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

    double& Madm1 = bconfig(MADM, BCO1);

    double Minf = Madm1 + mch;
    double e_bind = adm_inf - Minf;

    // center of mass defined like in https://arxiv.org/abs/1506.01689
    double COMx =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("COMx")()(ndom - 1),
                                          OUTER_BC) /
        adm_inf;
    double COMy =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("COMy")()(ndom - 1),
                                          OUTER_BC) /
        adm_inf;
    double COMz =
        space.get_domain(ndom - 1)->integ(syst.give_val_def("COMz")()(ndom - 1),
                                          OUTER_BC) /
        adm_inf;
    // END binary quantities

#ifdef FORMAT
#undef FORMAT
#endif
#define FORMAT1                                                     \
  std::setw(25) << std::right << std::setprecision(5) << std::fixed \
                << std::showpos
#define FORMAT                                                           \
  std::setw(25) << std::right << std::setprecision(5) << std::scientific \
                << std::showpos
    auto print_shells = [&](int dom_min, int dom_max) {
      int cnt = 1;
      for (int i = dom_min; i < dom_max; ++i) {
        std::string shell{"SHELL" + std::to_string(cnt) + " = "};
        Index pos(
            space.get_domain(i + 1)->get_radius().get_conf().get_dimensions());
        std::cout << FORMAT1 << shell
                  << space.get_domain(i + 1)->get_radius()(pos) << std::endl;
        cnt++;
      }
    };
    auto print_shell_mb = [&](auto& vec) {
      for (auto& e : vec) {
        std::cout << e << ",";
      }
    };

    int nshells1 = space.ADAPTEDNS - space.NS - 1;
    int nshells2 = space.OUTER - (space.BH + 2) - 1;
    auto res_r = space.get_domain(0)->get_nbr_points()(0);
    auto res_t = space.get_domain(0)->get_nbr_points()(1);
    auto res_p = space.get_domain(0)->get_nbr_points()(2);
    std::string header(22, '#');
    std::cout << header + " Neutron Star " + header + "\n";
    std::cout << FORMAT1 << "Center_COM = " << "(" << NS_x_com << ", 0, 0)\n"
              << FORMAT1 << "Coord R_IN = " << rin_ns << std::endl;
    // Print inner shells
    print_shells(space.NS + 1, space.ADAPTEDNS);
    std::cout << FORMAT1 << "Coord R = " << "[" << rmin << "," << rmax << "] ("
              << "[" << rmin * M2km << "," << rmax * M2km << "] km)"
              << std::endl;
    // Print outer shells
    print_shells(space.ADAPTEDNS + 1, space.BH - 1);
    std::cout << FORMAT1 << "Coord R_OUT = " << rout_ns << std::endl
              << FORMAT1 << "Areal R = " << areal_rns << " ["
              << areal_rns * M2km << "km]\n"
              << FORMAT << "Circumferential R = " << CR << " [" << CR * M2km
              << "km]\n"
              << FORMAT << "Mass Shedding = " << mass_shedding_parameter << endl
              << FORMAT1 << "NS Mb = " << MB1 << " (";
    print_shell_mb(baryonic_mass1);
    std::cout
        << ")\n"
        << FORMAT1 << "Isolated ADM Mass = " << bconfig(MADM, BCO1) << "\n"
        << FORMAT1 << "Quasi-local Madm = " << ql_madm1
        << " Diff:" << std::fabs(1. - ql_madm1 / bconfig(MADM, BCO1))
        << std::endl
        // quasi-local spin angular momentum
        << FORMAT1 << "Quasi-local S = " << ql_spinns
        << std::endl
        // dimensionless spin (constant, given by the imported single star!)
        << FORMAT1 << "Chi = " << ql_spinns / Madm1 / Madm1 << " ["
        << bconfig(CHI, BCO1)
        << "]\n"
        // angular frequency paramter of the star, describing the magnitude of
        // the spin component of the velocity field
        << FORMAT1 << "Omega = " << bconfig(OMEGA, BCO1) << std::endl
        << FORMAT1 << "Local P_y = " << NS_py << std::endl
        << FORMAT1 << "Local P_x = " << NS_px << std::endl
        << FORMAT << "Central Density = " << rhoc1 << std::endl
        << FORMAT << "Central log(h) = " << loghc1 << std::endl
        << FORMAT << "Central Pressure = " << pressc1 << std::endl
        << FORMAT << "Central dlog(h)/dx = " << dHdx1 << std::endl
        << FORMAT << "Central Euler Constant = " << euler1 << std::endl
        << FORMAT1 << "Integrated log(h) = " << intH1 << "\n\n";

    std::cout << header + " Black Hole " + header + "\n"
              << FORMAT1 << "Center_COM = " << "(" << BH_x_com << ", 0, 0)\n"
              << FORMAT1 << "Coord R_IN = " << rin_bh << std::endl
              << FORMAT1 << "Coord R = " << r_bh << " [" << r_bh * M2km
              << "km]\n";
    print_shells(space.BH + 2, space.OUTER - 1);
    std::cout << FORMAT1 << "Coord R_OUT = " << rout_bh << std::endl
              << FORMAT1 << "Areal R = " << areal_rbh << " ["
              << areal_rbh * M2km << "km]\n"
              << FORMAT1 << " LAPSE = [" << lapsemin << ", " << lapsemax
              << "]\n"
              << FORMAT1 << " PSI = [" << confmin << ", " << confmax << "]\n"
              << FORMAT1 << "Mirr = " << mirr << "\n"
              << FORMAT1 << "Mch = " << mch << " [" << bconfig(MCH, BCO2)
              << "]\n"
              << FORMAT1 << "Chi = " << ql_spinbh / (mch * mch) << " ["
              << bconfig(CHI, BCO2) << "]\n"
              << FORMAT1 << "S = " << ql_spinbh << std::endl
              << FORMAT1 << "Local P_y = " << BH_py << std::endl
              << FORMAT1 << "Local P_x = " << BH_px << std::endl
              << FORMAT1 << "Omega = " << bconfig(OMEGA, BCO2) << "\n\n";

    // boundaries of additional outer shells
    auto outer_shells = space.get_n_shells_outer();
    if (outer_shells > 0) {
      std::cout << FORMAT1 << "Outer shell bounds\n";
      print_shells(ndom - 1 - outer_shells, ndom - 1);
    }

    auto M1 = bconfig(MADM, BCO1);
    auto M2 = bconfig(MCH, BCO2);
    if (M2 > M1)
      std::swap(M1, M2);
    auto Mtot = M1 + M2;
    std::cout
        << header + " Binary " + header + "\n"
        << FORMAT1 << std::fixed << "RES = " << "[" << res_r << "," << res_t
        << "," << res_p
        << "]\n"
        // mass ratio, ratio of the ADM masses at infinity
        << FORMAT1 << "Q = " << M2 / M1 << std::endl
        << FORMAT1
        << std::setprecision(2)
        // Separation distance in geometrized units
        << "Separation = " << bconfig(DIST)
        << " ["
        // [Proper separation], (coordinate separation [km])
        << bconfig(DIST) / Mtot << "] (" << bconfig(DIST) * M2km << "km)"
        << std::endl
        // orbital angular frequency parameter
        << FORMAT1 << "Orbital Omega = " << bconfig(GOMEGA)
        << std::endl
        // Komar and ADM mass of the binary
        << FORMAT1 << "Komar mass = " << komar << std::endl
        << FORMAT1 << "Adm mass = " << adm_inf << ", Diff: " << e_diff
        << std::endl
        << FORMAT1 << "Total Mass = " << Minf << " [" << Mtot
        << "]\n"
        // ADM angular momentum of the binary
        << FORMAT1 << "Adm moment. = " << Jinf
        << std::endl
        // binding energy, defined by the gravitational mass difference at finite separation
        << FORMAT1 << "Binding energy = " << e_bind
        << std::endl
        // dimensionless orbital frequency
        << FORMAT1 << "Minf * Ome = " << Minf * bconfig(GOMEGA)
        << std::endl
        // dimensionless binding energy
        << FORMAT1 << "E_b / Minf = " << e_bind / Minf
        << std::endl
        // ADM linear momentum
        << FORMAT << "ADM P_x = " << Px << std::endl
        << FORMAT << "ADM P_y = " << Py << std::endl
        << FORMAT << "ADM P_z = " << Pz
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

  if (argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<str: BHNS ID basename>.info "
              << std::endl;
    std::cerr << "Ex: ./reader converged.TOTAL_BC.9.info" << std::endl;
    std::_Exit(EXIT_FAILURE);
  }

  // Name of config.info file
  std::string in_filename = argv[1];
  kadath_config_boost<BIN_INFO> bconfig(in_filename);

  // setup the EOS
  EOS_initialize::init(bconfig, NODES::BCO1);
  const std::string eos_type =
      bconfig.eos<std::string>(EOS_PARAMS::EOSTYPE, NODES::BCO1);
  EOS_Function_Dispatcher::dispatch<reader_output>(bconfig, eos_type, bconfig);

  return EXIT_SUCCESS;
}
