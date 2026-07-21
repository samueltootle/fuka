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
#include <iostream>
#include <numeric>
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath_bin_bh.hpp"
#include "mpi.h"
using namespace Kadath;
using namespace Kadath::FUKA_Config;
using vec_d = std::vector<double>;
using ary_d = std::array<double, 2>;
using ary_i = std::array<int, 2>;

double M2km = 1.4769994423016508;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./reader /<path>/<str: BBH ID basename>.info "
                  << std::endl;
        std::cerr << "Ex: ./reader converged.TOTAL_BC.9.info" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    int output = 0;
    if (argv[2] != 0x0)
        output = std::stoi(argv[2]);

    //Name of config.info file
    std::string in_filename = argv[1];
    kadath_config_boost<BIN_INFO> bconfig(in_filename);

    std::string in_spacefile = bconfig.space_filename();

    FILE* fich = fopen(in_spacefile.c_str(), "r");
    Space_bin_bh space(fich);
    Scalar conf(space, fich);
    Scalar lapse(space, fich);
    Vector shift(space, fich);
    fclose(fich);
    if (std::isnan(bconfig.set(ADOT)))
        bconfig.set(ADOT) = 0.;

    Base_tensor basis(shift.get_basis());
    Metric_flat fmet(space, basis);

    const std::array<double, 2> x_nuc{bco_utils::get_center(space, space.BH1),
                                      bco_utils::get_center(space, space.BH2)};
    double xo = (x_nuc[0] + x_nuc[1]) / 2.;

    //setup coord fields
    CoordFields<Space_bin_bh> cf_generator(space);
    vec_ary_t coord_vectors = default_binary_vector_ary(space);

    update_fields(cf_generator, coord_vectors, {}, xo, x_nuc[0], x_nuc[1]);
    Vector CART(space, CON, basis);
    CART = cf_generator.cart();
    //end setup coord fields

    int ndom = space.get_nbr_domains();
    System_of_eqs syst(space);
    fmet.set_system(syst, "f");
    syst.add_cst("PI", M_PI);
    syst.add_cst("ome", bconfig(GOMEGA));
    syst.add_cst("P", conf);
    syst.add_cst("N", lapse);
    syst.add_cst("bet", shift);
    syst.add_cst("xaxis", bconfig(COM));
    syst.add_cst("yaxis", bconfig(COMY));
    syst.add_cst("adot", bconfig(ADOT));

    syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("r", CART);

    syst.add_cst("ex", *coord_vectors[EX]);
    syst.add_cst("ey", *coord_vectors[EY]);
    syst.add_cst("ez", *coord_vectors[EZ]);

    syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
    syst.add_cst("mp", *coord_vectors[BCO2_ROT]);

    syst.add_cst("sm", *coord_vectors[S_BCO1]);
    syst.add_cst("sp", *coord_vectors[S_BCO2]);

    syst.add_cst("einf", *coord_vectors[S_INF]);

    syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");
    syst.add_def("Ntilde = N / P^6");
    syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
    syst.add_def("B^i = bet^i + ome * Morb^i + adot * comr^i");

    syst.add_def(
        "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) /2. / "
        "Ntilde");

    syst.add_def("intPx  = A_ij * ex^j * einf^i / 8 / PI");
    syst.add_def("intPy  = A_ij * ey^j * einf^i / 8 / PI");
    syst.add_def("intPz  = A_ij * ez^j * einf^i / 8 / PI");
    syst.add_def(space.BH1 + 2, "intS  = A_ij * mm^i * sm^j   / 8 / PI");
    syst.add_def(space.BH2 + 2, "intS  = A_ij * mp^i * sp^j   / 8 / PI");
    syst.add_def("intA = P^4");
    syst.add_def("intMsq = intA / 16. / PI");

    syst.add_def(ndom - 1, "intJ = multr(A_ij * Morb^i * einf^j) / 8 / PI");
    syst.add_def(ndom - 1, "Madm = -dr(P) / 2 / PI");
    syst.add_def(ndom - 1, "Mk   =  dr(N) / 4 / PI");

    //diagnostic from https://arxiv.org/abs/1506.01689
    syst.add_def(ndom - 1, "COMx  = -3 * P^4 * einf^i * ex_i / 8. / PI");
    syst.add_def(ndom - 1, "COMy  = 3 * P^4 * einf^i * ey_i / 8. / PI");
    syst.add_def(ndom - 1, "COMz  = 3 * P^4 * einf^i * ez_i / 8. / PI");

    syst.add_def("Morb^i  = mg^i + xaxis * ey^i");

    syst.add_def("dtgamma = D_k bet^k + 6 / P * B^k * D_k P");

    ary_i nuc_doms{space.BH1, space.BH2};
    ary_d mirrs;
    ary_d areal_rs;
    ary_d mchs;
    ary_d chis;
    ary_d spins;
    ary_d xcom;
    for (int i : {0, 1}) {
        int dom = nuc_doms[i];

        auto A = space.get_domain(dom + 2)->integ(syst.give_val_def("intA")()(
                                                      dom + 2),
                                                  INNER_BC);
        areal_rs[i] = sqrt(A / 4. / acos(-1.));

        auto mirrsq = space.get_domain(dom + 2)->integ(syst.give_val_def(
                                                           "intMsq")()(dom + 2),
                                                       INNER_BC);
        mirrs[i] = sqrt(mirrsq);
        spins[i] = space.get_domain(dom + 2)->integ(syst.give_val_def("intS")()(
                                                        dom + 2),
                                                    INNER_BC);
        mchs[i] = std::sqrt(mirrsq + spins[i] * spins[i] / 4. / mirrsq);
        chis[i] = spins[i] / mchs[i] / mchs[i];
        // center-of-mass shifted center
        xcom[i] = x_nuc[i] + bconfig(COM);
    }
    ary_i shells;
    shells[0] = space.BH2 - 3;
    shells[1] = space.OUTER - 6 - shells[0];

    double mirr = std::accumulate(mirrs.begin(), mirrs.end(), 0.);
    double mch = std::accumulate(mchs.begin(), mchs.end(), 0.);

    //quantities measured at inf
    auto outer_dom = space.get_domain(ndom - 1);
    double dtgamma_inf =
        outer_dom->integ(syst.give_val_def("dtgamma")()(ndom - 1), OUTER_BC);
    double adm_inf =
        outer_dom->integ(syst.give_val_def("Madm")()(ndom - 1), OUTER_BC);
    double COMx =
        outer_dom->integ(syst.give_val_def("COMx")()(ndom - 1), OUTER_BC) /
        adm_inf;
    double COMy =
        outer_dom->integ(syst.give_val_def("COMy")()(ndom - 1), OUTER_BC) /
        adm_inf;
    double COMz =
        outer_dom->integ(syst.give_val_def("COMz")()(ndom - 1), OUTER_BC) /
        adm_inf;
    double komar =
        outer_dom->integ(syst.give_val_def("Mk")()(ndom - 1), OUTER_BC);

    // J at infinity :
    Val_domain Jinf_vd = syst.give_val_def("intJ")()(ndom - 1);
    double Jinf =
        outer_dom->integ(syst.give_val_def("intJ")()(ndom - 1), OUTER_BC);
    double mu = mchs[0] * mchs[1] / mch;
    double e_bind = adm_inf - mch;

    // Linear Momenta
    double Px =
        outer_dom->integ(syst.give_val_def("intPx")()(ndom - 1), OUTER_BC);
    double Py =
        outer_dom->integ(syst.give_val_def("intPy")()(ndom - 1), OUTER_BC);
    double Pz =
        outer_dom->integ(syst.give_val_def("intPz")()(ndom - 1), OUTER_BC);

    auto res_r = space.get_domain(0)->get_nbr_points()(0);
    auto res_t = space.get_domain(0)->get_nbr_points()(1);
    auto res_p = space.get_domain(0)->get_nbr_points()(2);

    //print useful data without headers in a row.
    if (output == 2) {
        std::cout << bconfig(BIN_RES) << "," << bconfig(DIST) << ",";
        std::cout << std::setprecision(15) << std::scientific << mirrs[0] << ","
                  << bco_utils::get_radius(space.get_domain(nuc_doms[0] + 1),
                                           EQUI)
                  << "," << mirrs[1] << ","
                  << bco_utils::get_radius(space.get_domain(nuc_doms[1] + 1),
                                           EQUI)
                  << "," << mu << "," << adm_inf << "," << komar << ","
                  << fabs(adm_inf - komar) / (adm_inf + komar) << "," << mirr
                  << "," << Jinf << "," << bconfig(GOMEGA) << "," << e_bind
                  << "," << bconfig(OMEGA, BCO1) << "," << bconfig(OMEGA, BCO2)
                  << "," << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                  << "," << bconfig(MCH, BCO1) << "," << bconfig(MCH, BCO2)
                  << "," << mch * bconfig(GOMEGA) << "," << e_bind / mu << ","
                  << Jinf / mu / mirr << "," << bconfig(DIST) / mch << ","
                  << bconfig(MCH, BCO2) / bconfig(MCH, BCO1) << std::endl;
        std::_Exit(EXIT_SUCCESS);
    }

#define FORMAT1                                                       \
    std::setw(25) << std::right << std::setprecision(5) << std::fixed \
                  << std::showpos
    auto print_shells = [&](int dom_min, int dom_max) {
        int cnt = 1;
        for (int i = dom_min; i < dom_max; ++i) {
            std::string shell{"SHELL" + std::to_string(cnt) + " = "};
            std::cout << FORMAT1 << shell
                      << bco_utils::get_radius(space.get_domain(i), OUTER_BC)
                      << std::endl;
            cnt++;
        }
    };
    std::array<std::string, 2> bh_str{" BH_MINUS ", " BH_PLUS "};
    ary_i bounds{space.BH2 - 1, space.OUTER - 1};
    auto M1 = bconfig(MCH, BCO1);
    auto M2 = bconfig(MCH, BCO2);
    if (M2 > M1)
        std::swap(M1, M2);

#ifdef FORMAT
#undef FORMAT
#endif

#define FORMAT                                                             \
    std::setw(25) << std::right << std::setprecision(5) << std::scientific \
                  << std::showpos
    std::string header(22, '#');
    for (int i = 0; i <= 1; ++i) {
        auto [lapsemin, lapsemax] =
            bco_utils::get_field_min_max(lapse, nuc_doms[i] + 2, INNER_BC);
        auto [confmin, confmax] =
            bco_utils::get_field_min_max(conf, nuc_doms[i] + 2, INNER_BC);
        std::string idx = std::to_string(i + 1);

        std::cout << header + bh_str[i] + header + "\n";
        std::cout << FORMAT1 << "Center_COM = " << "(" << xcom[i] << ", 0, 0)\n"
                  << FORMAT1 << "Coord R_IN = "
                  << bco_utils::get_radius(space.get_domain(nuc_doms[i]), EQUI)
                  << '\n'
                  << FORMAT1 << "Coord R = "
                  << bco_utils::get_radius(space.get_domain(nuc_doms[i] + 1),
                                           EQUI)
                  << '\n';
        print_shells(nuc_doms[i] + 2, bounds[i]);
        std::cout << FORMAT1 << "Coord R_OUT = "
                  << bco_utils::get_radius(space.get_domain(bounds[i]), EQUI)
                  << "\n"
                  << FORMAT1 << "Areal R = " << areal_rs[i] << '\n'
                  << FORMAT1 << " LAPSE = " << "[" << lapsemin << ", "
                  << lapsemax << "]\n"
                  << FORMAT1 << " PSI = " << "[" << confmin << ", " << confmax
                  << "]\n"
                  << FORMAT1 << "Mirr = " << mirrs[i] << "[" << bconfig(MIRR, i)
                  << "]\n"
                  << FORMAT1 << "Mch = " << mchs[i] << "[" << bconfig(MCH, i)
                  << "]\n"
                  << FORMAT1 << "Chi = " << chis[i] << "[" << bconfig(CHI, i)
                  << "]\n"
                  << FORMAT << "S = " << spins[i] << std::endl
                  << FORMAT << "Omega = " << bconfig(OMEGA, i) << "\n";
    }
    std::cout << header + " Binary " + header + '\n'
              << FORMAT1 << std::fixed << "RES = " << "[" << res_r << ","
              << res_t << "," << res_p << "]\n"
              << FORMAT1 << "Q = " << M2 / M1 << std::endl
              << FORMAT1 << std::setprecision(2)
              << "Separation = " << bconfig(DIST) << " [" << bconfig(DIST) / mch
              << "] (" << bconfig(DIST) * M2km << "km)" << std::endl
              << FORMAT << "Orbital Omega = " << bconfig(GOMEGA) << std::endl
              << FORMAT << "Komar mass = " << komar << std::endl
              << FORMAT << "Adm mass = " << adm_inf
              << ", Diff: " << 2. * fabs(adm_inf - komar) / (adm_inf + komar)
              << std::endl
              << FORMAT1 << "Total Mirr = " << mirr << std::endl
              << FORMAT1 << "Total Mch = " << mch << std::endl
              << FORMAT << "Adm moment. = " << Jinf << std::endl
              << FORMAT << "Binding energy = " << e_bind << std::endl
              << FORMAT << "M * Ome = " << mch * bconfig(GOMEGA) << std::endl
              << FORMAT << "E_b / mu = " << e_bind / mu << std::endl
              << FORMAT << "Px = " << Px << '\n'
              << FORMAT << "Py = " << Py << '\n'
              << FORMAT << "Pz = " << Pz << '\n'
              << FORMAT << "COMx = " << bconfig(COM) << ", A-COMx = " << COMx
              << std::endl
              << FORMAT << "COMy = " << bconfig(COMY) << ", A-COMy = " << COMy
              << std::endl
              << FORMAT << "A-COMz = " << COMz << std::endl;

    return EXIT_SUCCESS;
}
