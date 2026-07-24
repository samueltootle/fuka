/*
 * Copyright 2022
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
#pragma once
#include "bco_utilities.hpp"
#include "mpi.h"

/**
 * \addtogroup BH_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template <typename config_t, typename space_t>
int bh_3d_xcts_solver<config_t, space_t>::von_Neumann_stage(
    std::string stage_text) {
    int exit_status = EXIT_SUCCESS;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto current_file = bconfig.config_filename_abs();
    if (!bconfig.control(RESOLVE) && solution_exists("TOTAL_BC")) {
        if (bconfig.config_filename_abs() == current_file) {
            if (rank == 0)
                std::cout << "Solved previously: "
                          << bconfig.config_filename_abs() << std::endl;
            return EXIT_SUCCESS;
        } else {
            return RELOAD_FILE;
        }
    }

    if (rank == 0)
        std::cout << "############################" << std::endl
                  << "Total system with von Neumann BC" << std::endl
                  << "############################" << std::endl;
    bconfig.set(MIRR) =
        Kadath::bco_utils::mirr_from_mch(bconfig(CHI), bconfig(MCH));

    double xo = 0.;

    System_of_eqs syst(space, 0, ndom - 1);
    syst_init(syst);

    // in the event we want a linear boost
    syst.add_cst("xboost", bconfig(BVELX));
    syst.add_cst("yboost", bconfig(BVELY));

    syst.add_def("intS  = A_ij * mg^i * sm^j / 8. / PI");

    // definition of constraint equations - conformally flat/maximal slice only
    // K = 0 -> R = 0; \partial_t(metric) = 0
    syst.add_def("eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");

    // Add equations and continuity constraints to the numerical space
    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    space.add_bc_inf(syst, "NP = 1");
    space.add_bc_inf(syst, "P = 1");
    space.add_bc_inf(syst, "bet^i = xboost * ex^i + yboost * ey^i");

    // excision boundary conditions
    space.add_bc_bh(syst, "dn(NP) = 0");
    space.add_bc_bh(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_bc_bh(syst, "bet^i = N / P^2 * sm^i + ome * mg^i");

    // equations fixing the radius of the BH based on (MIRR)
    // and the spin angular momentum based on fixed dimensionless spin (CHI)
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
        print_diagnostics(syst, ite, conv);
    while (!endloop) {
        endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

        std::stringstream ss;
        ss << "bh_total_bc_" << ite - 1;
        bconfig.set_filename(ss.str());
        if (rank == 0) {
            print_diagnostics(syst, ite, conv);
            if (bconfig.control(CHECKPOINT)) {
                checkpoint();
            }
        }

        update_fields_co(cfields, coord_vectors, {}, xo, &syst);
        ite++;
        check_max_iter_exceeded(rank, ite, conv);
    }

    bconfig.set(RMID) =
        Kadath::bco_utils::get_radius(space.get_domain(1), OUTER_BC);
    bconfig.set(FIXED_LAPSE) =
        Kadath::bco_utils::get_boundary_val(2, lapse, INNER_BC);
    bconfig.set_filename(converged_filename("TOTAL_BC"));
    if (rank == 0)
        checkpoint();

    return exit_status;
}

template <typename config_t, typename space_t>
int bh_3d_xcts_solver<config_t, space_t>::fixed_lapse_stage() {
    int exit_status = EXIT_SUCCESS;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto current_file = bconfig.config_filename_abs();
    if (!bconfig.control(RESOLVE) && solution_exists("TOTAL")) {
        if (bconfig.config_filename_abs() == current_file) {
            if (rank == 0)
                std::cout << "Solved previously: "
                          << bconfig.config_filename_abs() << std::endl;
            return EXIT_SUCCESS;
        } else {
            return RELOAD_FILE;
        }
    }

    if (rank == 0)
        std::cout << "############################" << std::endl
                  << "Total system with Fixed Lapse BC" << std::endl
                  << "############################" << std::endl;
    bconfig.set(MIRR) =
        Kadath::bco_utils::mirr_from_mch(bconfig(CHI), bconfig(MCH));

    double xo = 0.;

    System_of_eqs syst(space, 0, ndom - 1);
    syst_init(syst);

    // in the event we want a linear boost
    syst.add_cst("xboost", bconfig(BVELX));
    syst.add_cst("yboost", bconfig(BVELY));
    syst.add_cst("n0", bconfig(FIXED_LAPSE));

    syst.add_def("intS  = A_ij * mg^i * sm^j / 8. / PI");

    // definition of constraint equations - conformally flat/maximal slice only
    // K = 0 -> R = 0; \partial_t(metric) = 0
    syst.add_def("eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");

    // Add equations and continuity constraints to the numerical space
    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    space.add_bc_inf(syst, "NP = 1");
    space.add_bc_inf(syst, "P = 1");
    space.add_bc_inf(syst, "bet^i = xboost * ex^i + yboost * ey^i");

    // excision boundary conditions
    space.add_bc_bh(syst, "N = n0");
    space.add_bc_bh(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_bc_bh(syst, "bet^i = N / P^2 * sm^i + ome * mg^i");

    // equations fixing the radius of the BH based on (MIRR)
    // and the spin angular momentum based on fixed dimensionless spin (CHI)
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
        print_diagnostics(syst, ite, conv);
    while (!endloop) {
        endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

        std::stringstream ss;
        ss << "bh_total_" << ite - 1;
        bconfig.set_filename(ss.str());
        if (rank == 0) {
            print_diagnostics(syst, ite, conv);
            if (bconfig.control(CHECKPOINT)) {
                checkpoint();
            }
        }

        update_fields_co(cfields, coord_vectors, {}, xo, &syst);
        ite++;
        check_max_iter_exceeded(rank, ite, conv);
    }

    bconfig.set(RMID) =
        Kadath::bco_utils::get_radius(space.get_domain(1), OUTER_BC);
    bconfig.set_filename(converged_filename("TOTAL"));
    if (rank == 0)
        checkpoint();

    return exit_status;
}

template <typename config_t, typename space_t>
int bh_3d_xcts_solver<config_t, space_t>::binary_boost_stage(
    kadath_config_boost<BIN_INFO>& binconfig,
    const size_t bco) {
    int exit_status = EXIT_SUCCESS;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
        std::cout << "############################" << std::endl
                  << "Binary boost using von Neumann BC" << std::endl
                  << "############################" << std::endl;
    bconfig.set(MIRR) =
        Kadath::bco_utils::mirr_from_mch(bconfig(CHI), bconfig(MCH));

    // generate filename string unique to this binary setup
    std::stringstream ss;
    ss << "BIN_BOOST" << "_" << binconfig(DIST) << "_" << binconfig(GOMEGA);
    auto const boost_converged_filename{ss.str()};

    // check if we have a non-boosted solution
    if (!bconfig.control(RESOLVE) &&
        solution_exists(boost_converged_filename)) {
        if (rank == 0)
            std::cout << "Solved previously: " << bconfig.config_filename_abs()
                      << std::endl;
        return EXIT_SUCCESS;
    }

    double xo = 0.;

    update_fields_co(cfields, coord_vectors, {}, xo);

    System_of_eqs syst(space, 0, ndom - 1);
    syst_init(syst);

    syst.add_def("intS = A_ij * mg^i * sm^j / 8. / PI");

    // boost based on binary config
    syst.add_cst("omega_boost", binconfig(GOMEGA));
    syst.add_def("B^i = bet^i + omega_boost * mg^i");

    // definition of constraint equations - conformally flat/maximal slice only
    // K = 0 -> R = 0; \partial_t(metric) = 0
    syst.add_def("eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
    syst.add_def("eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(
        "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j "
        "Ntilde");

    // Add equations and continuity constraints to the numerical space
    space.add_eq(syst, "eqNP = 0", "N", "dn(N)");
    space.add_eq(syst, "eqP = 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i = 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    space.add_bc_inf(syst, "NP = 1");
    space.add_bc_inf(syst, "P = 1");
    space.add_bc_inf(syst, "bet^i = 0");

    // excision boundary conditions
    space.add_bc_bh(syst, "dn(NP) = 0");
    space.add_bc_bh(
        syst,
        "sm^j * D_j P + P / 4 * D^j sm_j + A_ij * sm^i * sm^j / P^3 / 4 = 0");
    space.add_bc_bh(syst, "B^i = N / P^2 * sm^i + ome * mg^i");

    // equations fixing the radius of the BH based on (MIRR)
    // and the spin angular momentum based on fixed dimensionless spin (CHI)
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
        print_diagnostics(syst, ite, conv);
    while (!endloop) {
        endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

        std::stringstream ss;
        ss << "bh_bin_boost_" << ite - 1;
        bconfig.set_filename(ss.str());
        if (rank == 0) {
            print_diagnostics(syst, ite, conv);
            if (bconfig.control(CHECKPOINT)) {
                checkpoint();
            }
        }

        update_fields_co(cfields, coord_vectors, {}, xo, &syst);
        ite++;
        check_max_iter_exceeded(rank, ite, conv);
    }

    bconfig.set(RMID) =
        Kadath::bco_utils::get_radius(space.get_domain(1), OUTER_BC);
    bconfig.set_filename(converged_filename(boost_converged_filename));
    if (rank == 0)
        checkpoint();

    return exit_status;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath