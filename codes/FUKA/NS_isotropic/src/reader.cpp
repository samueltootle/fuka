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
// FUKA includes
#include "Configurator/config_bco.hpp"
#include "EOS/EOS.hh"
#include "EOS/FUKA_EOS_Utilities.hh"
#include "bco_utilities.hpp"

// Kadath includes
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_polar.hpp"

// C++ includes
#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <type_traits>

using namespace Kadath;
using namespace Kadath::FUKA_Config;
using namespace Kadath::FUKA_EOS;

// conversion from solar mass to km
constexpr double M2km = 1.4769994423016508;

template <class eos_t, typename config_t>
void reader_2d_diffrot(config_t bconfig) {
  // load the space (and thus the domain setup)
  auto spacein = bconfig.space_filename();
  FILE* ff1 = fopen(spacein.c_str(), "r");
  Space_polar_adapted space(ff1);
  Scalar ome(space);
  ome = bconfig(BCO_PARAMS::OMEGA);
  ome.std_base();

  // load the fields defined on the space
  Scalar lap_Aterm(space, ff1);
  Scalar nu(space, ff1);
  Scalar logh(space, ff1);
  Scalar lap_Bterm(space, ff1);
  Scalar lap_wterm(space, ff1);
  if (bconfig.set_field(BCO_FIELDS::DIFF_OMEGA))
    ome = Scalar(space, ff1);
  fclose(ff1);
  lap_wterm.affect_parameters();
  lap_wterm.set_parameters()->set_m_quant() = 1;
  lap_wterm.std_base();

  // central values of the matter fields
  double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);
  double hc = std::exp(loghc);
  double nc = EOS<eos_t, DENSITY>::get(hc);
  double pc = EOS<eos_t, PRESSURE>::get(hc);

  // minimal and maximal radius of the adapted surface domain
  auto [rmin, rmax] = bco_utils::get_rmin_rmax(space, 1);
  // inner radius of the nucleus
  double rin1 = bco_utils::get_radius(space.get_domain(0), OUTER_BC);

  int ndom = space.get_nbr_domains();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  syst.add_cst("H", logh);
  syst.add_cst("nu", nu);
  syst.add_cst("lapAterm", lap_Aterm);
  syst.add_cst("lapBterm", lap_Bterm);
  syst.add_cst("lapwterm", lap_wterm);
  if (bconfig.set_field(BCO_FIELDS::DIFF_OMEGA))
    syst.add_cst("Omega", ome);
  else
    syst.add_cst("Omega", bconfig(BCO_PARAMS::OMEGA));

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  set_eos_ope_struct<eos_t>()(syst, p);

  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");
  syst.add_def("delta = h - eps - 1.");
  syst.add_def("edens = rho * (1 + eps)");

  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_def("B = (divrsint(lapBterm) + 1) / N");
  syst.add_def("w = divrsint(lapwterm)");
  syst.add_def("Brsint = multrsint(B)");
  Scalar Brsint(syst.give_val_def("Brsint")());
  Scalar psi(log(Brsint));
  psi.std_base();

  // syst.add_def("psi = log(Brsint)");
  syst.add_cst("psi", psi);

  syst.add_def("diffAB = B^2 - A^2");

  // eq. 4.21 - Full ADM surface integral at spatial infinity.  This does not give accurate
  // values for some reason.
  syst.add_def(ndom - 1,
               "intMadmFULL = - (dr(A^2 + B^2) + divr(B^2 - A^2)) / 4 / 4piG ");

  // If we assume that at infinity A = B = 1, we can obtain two "equivalent expressions"
  // However, we find that intMadmA does not give as accurate of results as intMadmB
  // intMadmB, however, gives very accurate results as compared to Mkomar and MADM that
  // has been computed for the same configuration using the full 3D code.
  syst.add_def(ndom - 1, "intMadmA = - (dr(A)) / 4piG ");
  syst.add_def(ndom - 1, "intMadmB = - (dr(B)) / 4piG ");

  // eq. 4.14 Gourgoulhon (4.15 can also be used)
  syst.add_def(ndom - 1,
               "intMk = B * (dr(N) - multrsint(multrsint(B^2) / 2 / N * w * "
               "dr(w)))  / 4piG");

  // eq 4.40 evaluated at spatial infinity, Gourgoulhon
  syst.add_def(ndom - 1, "intJ = -multrsint(multrsint(dr(w))) / 4 / 4piG");

  // (3.31) Upper Phi component of U vector (r, theta = 0)
  syst.add_def("UphiU = (Omega - w) / N");
  // (3.32)
  syst.add_def("U = multrsint(B * UphiU)");
  syst.add_def("Usq = U*U");
  syst.add_def("Wsq = 1 / (1 - Usq)");
  syst.add_def("W = sqrt(Wsq)");
  // Lower Phi component of U vector = u . \xi
  // (3.85)
  syst.add_def("UphiL = multrsint(multrsint(B^2 * UphiU))");
  syst.add_def("j = Wsq / N * UphiL");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
      // in the star the constraint equations are sourced by the matter
      case 0:
      case 1:
        // sources rescaled by p/rho as in Papenfort2021
        syst.add_def(d, "E = Wsq * (press * h - Wsq * press * delta / Wsq)");
        syst.add_def(d, "Srrtt = press * delta");
        syst.add_def(d, "Spp = delta * press * (1 + Usq) + E * Usq");
        syst.add_def(d, "S = 2 * Srrtt + Spp");

        // Phi component of the pressure
        syst.add_def(d, "pphi = multrsint(B * (E + Srrtt) * U)");

        // Volume integral for Angular momentum, eq 4.38 Gourgoulhon
        syst.add_def(d, "intJV = pphi * A^2 * B * 4piG / 2");

        // constraint equations
        syst.add_def(d,
                     "DDA = -scal(grad(nu), grad(nu)) + 2 * 4piG * A^2 * Spp");
        syst.add_def(d, "intDDA = - lap(A) * multrsint(A^2) * multr(B)");

        // definition for the baryonic mass integral, eq 4.5
        syst.add_def(d, "intMb = W * rho * A^2 * B * 4piG / 2");

        // Momentum volume integral over the star.
        // P^6 doesn't yield reasonable results as expected.
        syst.add_def(d, "vintJ = rho * j * A^2 * B * 4piG / 2");
        syst.add_def(d, "intT = Omega * vintJ");
        syst.add_def(d, "inteps  = intMb * eps(h)");

        // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
        syst.add_def(d, "firstint = H + log(N) - 0.5 * log(Wsq)");

        syst.add_def(d,
                     "GRV2 = divrsint(2 * 4piG * A^2 * Spp"
                     " + delta * 3 * multrsint(multrsint(B^2)) / 4 / N^2 * "
                     "scal(grad(w), grad(w))"
                     " - delta * scal(grad(nu), grad(nu))) / 2");

        // syst.add_def(d, "lam3 = 4piG * (3 * press + (edens + press) * Usq * Wsq) * A^2 * B");
        break;
      // outside the matter is absent and the sources are zero
      default:
        // syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(lapAterm))") ;
        syst.add_def(d, "DDA = -delta * scal(grad(nu), grad(nu))");
        syst.add_def(d, "intDDA = - lap2(A) * 2 / 4piG");
        syst.add_def(
            d,
            "GRV2 = divrsint(3 * multrsint(multrsint(B^2)) / 4 / N^2 * "
            "scal(grad(w), grad(w)) - scal(grad(nu), grad(nu))) / 2");
        break;
    }
  }
  // syst.add_def(1, "OmegaK = w + dr(w) / 2 / psi "
  //         "+ sqrt(N^2 * dr(nu) / dr(psi) /  Brsint^2 + (dr(w) / dr(psi) / 2)^2)");
  // syst.add_def(1, "OmegaK = w + dr(w) / 2 / psi");
  // syst.add_def(0, "OmegaK = w + dr(w)");
  syst.sec_member();

  double VMadm = 0;
  Scalar intDDA(syst.give_val_def("intDDA")());
  intDDA.coef_i();

  double VJadm = 0;
  Scalar intJV(syst.give_val_def("intJV")());
  intJV.coef_i();

  // P / rho
  Scalar P_o_rho(syst.give_val_def("delta")());
  P_o_rho.coef_i();

  double baryonic_mass = 0;
  Scalar intMb(syst.give_val_def("intMb")());
  intMb.coef_i();

  double GRV2 = 0;
  Scalar intGRV2(syst.give_val_def("GRV2")());
  intGRV2.coef_i();

  double T_integral = 0.0;
  Scalar Tint = syst.give_val_def("intT")();

  double J_vintegral = 0.0;
  Scalar Jvol = syst.give_val_def("vintJ")();

  double eps_integral = 0.0;
  Scalar epsint = syst.give_val_def("inteps")();

  // double LAM3=0;
  // Scalar intLAM3(syst.give_val_def("lam3")());
  // intLAM3.coef_i();

  for (int i = 0; i < 2; ++i) {
    VMadm += intDDA(i).integ_volume();
    // LAM3 += intLAM3(i).integ_volume();

    // To obtain the rescaled angular momentum
    // We need to compute pphi / (P / rho)
    // However, since the surface is defined by P = rho = 0, this produces
    // NaNs.  Therefore, we compute this manually here such
    // that we can assert that P/rho on the boundary is zero
    Val_domain J(intJV(i));
    Val_domain Porho(P_o_rho(i));
    Val_domain GRV2_vd(intGRV2(i));

    Index pos(space.get_domain(i)->get_nbr_points());
    Val_domain J_o_Porho(J);
    do {
      double j = J(pos);
      double porho = Porho(pos);
      double grv2_val = GRV2_vd(pos);
      if (std::fabs(porho) <= 1e-15) {
        J_o_Porho.set(pos) = 0.;
        // GRV2_vd.set(pos) = 0.;
      } else {
        J_o_Porho.set(pos) = j / porho;
        // GRV2_vd.set(pos) = grv2_val / porho;
      }
    } while (pos.inc());
    VJadm += J_o_Porho.integ_volume();
    baryonic_mass += intMb(i).integ_volume();
    T_integral += Tint(i).integ_volume();
    J_vintegral += Jvol(i).integ_volume();
    eps_integral += epsint(i).integ_volume();
    double tmp = GRV2_vd.integ_volume();
    cout << i << ": " << tmp << '\n';
    GRV2 += tmp;
  }

  intGRV2.coef_i();
  for (int i = 2; i < ndom; ++i) {
    auto tmp = intGRV2(i).integ_volume();
    cout << i << ": " << tmp << '\n';
    GRV2 += tmp;
  }

  cout << "GRV2: " << GRV2 << '\n';
  // cout << "LAM3: " << LAM3 << '\n';

  Val_domain integMadm(syst.give_val_def("intMadmFULL")()(ndom - 1));
  double MadmFULL = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
  Val_domain integMadmA(syst.give_val_def("intMadmA")()(ndom - 1));
  double MadmA = space.get_domain(ndom - 1)->integ(integMadmA, OUTER_BC);
  Val_domain integMadmB(syst.give_val_def("intMadmB")()(ndom - 1));
  double MadmB = space.get_domain(ndom - 1)->integ(integMadmB, OUTER_BC);

  // Komar mass at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // ADM angular momentum at infinity
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

  T_integral *= 0.5;
  double W_be = T_integral + eps_integral + baryonic_mass - MadmB;
  //std::cout << "T_integral: " << T_integral << std::endl;
  //std::cout << "W_be: " << W_be << std::endl;
  //std::cout << "Beta: " << T_integral / W_be << std::endl;

  auto npts = space.get_domain(1)->get_nbr_points();

  Index pos_eq(npts);
  pos_eq.set(0) = npts(0) - 1;  /// Set to outer radius
  pos_eq.set(1) = npts(1) - 1;  /// Set theta to be on the xy plane.

  Index pos_pole(npts);
  pos_pole.set(0) = npts(0) - 1;  /// Set to outer radius

  Scalar logh_dr(logh.der_r());
  double mass_shedding_parameter = logh_dr(1)(pos_eq) / logh_dr(1)(pos_pole);

  auto B(syst.give_val_def("B")());
  auto N(syst.give_val_def("N")()(1));
  auto r(space.get_domain(1)->get_radius());
  auto [Bmin, Bmax] = bco_utils::get_field_min_max(B, 2, INNER_BC);
  double CR = B(1)(pos_eq) * r(pos_eq);

#ifdef FORMAT
#undef FORMAT
#endif
#define FORMAT                                                      \
  std::setw(25) << std::right << std::setprecision(5) << std::fixed \
                << std::showpos
  auto print_shells = [&](int dom_min, int dom_max) {
    int cnt = 1;
    for (int i = dom_min; i < dom_max; ++i) {
      std::string shell{"SHELL" + std::to_string(cnt) + " = "};
      std::cout << FORMAT << shell
                << bco_utils::get_radius(space.get_domain(i), OUTER_BC)
                << std::endl;
      cnt++;
    }
  };

  auto res_r = space.get_domain(0)->get_nbr_points()(0);
  auto res_t = space.get_domain(0)->get_nbr_points()(1);

  // output to stdout
  std::cout << FORMAT << "RES = " << "[" << res_r << "," << res_t << "]\n"
            << FORMAT << "Coord R_IN = " << rin1 << std::endl
            << FORMAT << "Coord R = " << "[" << rmin << ", " << rmax << "]\n";
  std::cout << FORMAT << "Coord R_OUT = "
            << bco_utils::get_radius(space.get_domain(2), OUTER_BC) << "\n";
  print_shells(3, ndom - 1);
  cout << endl;

  std::cout << FORMAT << "Circumferential R = " << CR << " [" << CR * M2km
            << "km]\n"
            << FORMAT << "Mass Shedding = " << mass_shedding_parameter << "\n"
            << FORMAT << "Baryonic Mass = " << baryonic_mass << std::endl
            << FORMAT << "ADM Mass = " << MadmB << " [" << MadmA << ", "
            << MadmFULL << "]\n"
            << FORMAT << "ADM Momentum = " << J << " [" << VJadm << "]\n"
            << FORMAT << "Chi = " << J / MadmB / MadmB << " [" << bconfig(CHI)
            << "]\n"
            << FORMAT << "Omega = " << bconfig(OMEGA) << std::endl
            << FORMAT << std::scientific << "Central Density = " << nc
            << std::endl
            << FORMAT << std::scientific << "Central h = " << hc << std::endl
            << FORMAT << std::scientific << "Central log(h) = " << loghc
            << std::endl
            << FORMAT << std::scientific << "Central Pressure = " << pc
            << "\n"
            << FORMAT << "Beta: " << T_integral / W_be << "\n\n";
  // << FORMAT << std::scientific << "Central dlog(h)/dx = " << central_dHdx << std::endl
  // << FORMAT << std::scientific << "Central Euler Constant = "<< central_euler << std::endl
  // << FORMAT << "Integrated log(h) = "    << H_integral << "\n\n";

  std::cout << FORMAT << "Mk = " << Mk << std::scientific
            << ", Diff: " << 2. * fabs(MadmB - Mk) / (MadmB + Mk) << " ["
            << 2. * fabs(MadmA - Mk) / (MadmA + Mk) << ", "
            << 2. * fabs(MadmFULL - Mk) / (MadmFULL + Mk) << "]\n";
  // << FORMAT << "Px = "   << Px   << std::endl
  // << FORMAT << "Py = "   << Py   << std::endl
  // << FORMAT << "Pz = "   << Pz   << std::endl;
}

template <class eos_t, typename config_t>
void reader_2d_norot(config_t bconfig) {

  // load the space (and thus the domain setup)
  auto spacein = bconfig.space_filename();
  FILE* ff1 = fopen(spacein.c_str(), "r");
  Space_polar_adapted space(ff1);

  // load the fields defined on the space
  Scalar lap_Aterm(space, ff1);
  Scalar nu(space, ff1);
  Scalar logh(space, ff1);
  Scalar lap_Bterm(space, ff1);
  fclose(ff1);

  // central values of the matter fields
  double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);
  double hc = std::exp(loghc);
  double nc = EOS<eos_t, DENSITY>::get(hc);
  double pc = EOS<eos_t, PRESSURE>::get(hc);

  // minimal and maximal radius of the adapted surface domain
  auto [rmin, rmax] = bco_utils::get_rmin_rmax(space, 1);
  // inner radius of the nucleus
  double rin1 = bco_utils::get_radius(space.get_domain(0), OUTER_BC);

  int ndom = space.get_nbr_domains();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  syst.add_cst("H", logh);
  syst.add_cst("nu", nu);
  syst.add_cst("lapAterm", lap_Aterm);
  syst.add_cst("lapBterm", lap_Bterm);

  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_def("B = (divrsint(lapBterm) + 1) / N");
  syst.add_def("Brsint = multrsint(B)");
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
  syst.add_def("delta = h - eps - 1.");

  syst.add_def(ndom - 1,
               "intMadmFULL = - (dr(A^2 + B^2) + divr(B^2 - A^2)) / 4 / 4piG ");

  // If we assume that at infinity A = B = 1, we can obtain two "equivalent expressions"
  // However, we find that intMadmA does not give as accurate of results as intMadmB
  // intMadmB, however, gives very accurate results as compared to Mkomar and MADM that
  // has been computed for the same configuration using the full 3D code.
  syst.add_def(ndom - 1, "intMadmA = - (dr(A)) / 4piG ");
  syst.add_def(ndom - 1, "intMadmB = - (dr(B)) / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
      // in the star the constraint equations are sourced by the matter
      case 0:
      case 1:
        // sources
        syst.add_def(d, "E = press * h - press * delta");
        syst.add_def(d, "S = delta * 3 * press");
        syst.add_def(d, "Spp = press * delta");

        // constraint equations
        syst.add_def(d,
                     "eqnu = delta * ( lap(nu) + scal(grad(nu), "
                     "grad(lapAterm)) ) - 4piG * A^2 * (E + S)");
        syst.add_def(d,
                     "eqlapAterm = delta * ( lap2(lapAterm) + scal(grad(nu), "
                     "grad(nu)) ) - 2 * 4piG * A^2 * Spp");
        // Extra...
        // syst.add_def(d, "eqNA = dr(drNA) + 3 * divr(drNA) - 4 * 4piG * NA * A^2 * press") ;

        // // definition for the baryonic mass integral
        syst.add_def(d, "intMb = rho * A^2 * B * 4piG / 2");
        syst.add_def(d, "intDDA = - lap2(A) * multrsint(A^2) * multr(B)");

        // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
        syst.add_def(d, "firstint = H + log(N)");

        syst.add_def(d,
                     "GRV2 = divrsint(2 * 4piG * A^2 * Spp"
                     " - delta * scal(grad(nu), grad(nu))) / 2");

        break;
      // outside the matter is absent and the sources are zero
      default:
        syst.add_eq_full(d, "H = 0");
        syst.add_def(d, "DDA = -delta * scal(grad(nu), grad(nu))");
        syst.add_def(d, "intDDA = - lap2(A) * 2 / 4piG");
        syst.add_def(d, "GRV2 = divrsint(scal(grad(nu), grad(nu))) / 2");

        syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(lapAterm))");
        syst.add_def(d,
                     "eqlapAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))");
        break;
    }
  }

  double baryonic_mass = 0;
  Scalar intMb(syst.give_val_def("intMb")());
  intMb.coef_i();

  double VMadm = 0;
  Scalar intDDA(syst.give_val_def("intDDA")());
  intDDA.coef_i();

  for (int i = 0; i < 2; ++i) {
    VMadm += intDDA(i).integ_volume();
    baryonic_mass += intMb(i).integ_volume();
  }
  cout << "VMadm: " << VMadm << endl;

  Val_domain integMadm(syst.give_val_def("intMadmFULL")()(ndom - 1));
  double MadmFULL = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
  Val_domain integMadmA(syst.give_val_def("intMadmA")()(ndom - 1));
  double MadmA = space.get_domain(ndom - 1)->integ(integMadmA, OUTER_BC);
  Val_domain integMadmB(syst.give_val_def("intMadmB")()(ndom - 1));
  double MadmB = space.get_domain(ndom - 1)->integ(integMadmB, OUTER_BC);

  // Komar mass at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  auto npts = space.get_domain(1)->get_nbr_points();

  Index pos_eq(npts);
  pos_eq.set(0) = npts(0) - 1;  /// Set to outer radius
  pos_eq.set(1) = npts(1) - 1;  /// Set theta to be on the xy plane.

  Index pos_pole(npts);
  pos_pole.set(0) = npts(0) - 1;  /// Set to outer radius

  auto B(syst.give_val_def("A")()(1));
  auto r(space.get_domain(1)->get_radius());
  double CR = B(pos_eq) * r(pos_eq);

#ifdef FORMAT
#undef FORMAT
#endif
#define FORMAT                                                      \
  std::setw(25) << std::right << std::setprecision(5) << std::fixed \
                << std::showpos
  auto print_shells = [&](int dom_min, int dom_max) {
    int cnt = 1;
    for (int i = dom_min; i < dom_max; ++i) {
      std::string shell{"SHELL" + std::to_string(cnt) + " = "};
      std::cout << FORMAT << shell
                << bco_utils::get_radius(space.get_domain(i), OUTER_BC)
                << std::endl;
      cnt++;
    }
  };

  auto res_r = space.get_domain(0)->get_nbr_points()(0);
  auto res_t = space.get_domain(0)->get_nbr_points()(1);

  // output to stdout
  std::cout << FORMAT << "RES = " << "[" << res_r << "," << res_t << "]\n"
            << FORMAT << "Coord R_IN = " << rin1 << std::endl
            << FORMAT << "Coord R = " << "[" << rmin << ", " << rmax << "]\n";
  std::cout << FORMAT << "Coord R_OUT = "
            << bco_utils::get_radius(space.get_domain(2), OUTER_BC) << "\n";
  print_shells(3, ndom - 1);
  cout << endl;

  std::cout << FORMAT << "Circumferential R = " << CR << " [" << CR * M2km
            << "km]\n"
            << FORMAT << "Baryonic Mass = " << baryonic_mass << std::endl;
  std::cout << FORMAT << "ADM Mass = " << MadmB << " [" << MadmA << ", "
            << MadmFULL << "]\n"
            << FORMAT << std::scientific << "Central Density = " << nc
            << std::endl
            << FORMAT << std::scientific << "Central h = " << hc << std::endl
            << FORMAT << std::scientific << "Central log(h) = " << loghc
            << std::endl
            << FORMAT << std::scientific << "Central Pressure = " << pc
            << "\n\n";
  // << FORMAT << std::scientific << "Central dlog(h)/dx = " << central_dHdx << std::endl
  // << FORMAT << std::scientific << "Central Euler Constant = "<< central_euler << std::endl
  // << FORMAT << "Integrated log(h) = "    << H_integral << "\n\n";

  std::cout << FORMAT << "Mk = " << Mk << std::scientific
            << ", Diff: " << 2. * fabs(MadmB - Mk) / (MadmB + Mk) << std::endl;
  // << FORMAT << "Px = "   << Px   << std::endl
  // << FORMAT << "Py = "   << Py   << std::endl
  // << FORMAT << "Pz = "   << Pz   << std::endl;
}

template <class eos_t>
struct reader_select {
  template <typename config_t>
  void operator()(config_t& bconfig) {
    if (bconfig.set_field(BCO_FIELDS::LAP_BTERM) &&
        bconfig.set_field(BCO_FIELDS::LAP_WTERM)) {
      reader_2d_diffrot<eos_t>(bconfig);
    } else {
      reader_2d_norot<eos_t>(bconfig);
    }
  }
};

int main(int argc, char** argv) {
  // expecting a configuration file on execution
  if (argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  // load the configuration
  std::string ifilename{argv[1]};
  kadath_config_boost<BCO_NS_INFO> bconfig(ifilename);

  // setup the EOS
  EOS_initialize::init(bconfig);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);
  EOS_Function_Dispatcher::dispatch<reader_select>(bconfig, eos_type, bconfig);

  return EXIT_SUCCESS;
}
