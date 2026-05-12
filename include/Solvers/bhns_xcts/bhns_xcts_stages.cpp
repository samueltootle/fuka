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
#include "bco_utilities.hpp"
#include "mpi.h"

/**
 * \addtogroup BHNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template <class eos_t, typename config_t, typename space_t>
int bhns_xcts_solver<eos_t, config_t, space_t>::hydrostatic_equilibrium_stage(
    const std::string stage_text) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_status = EXIT_SUCCESS;

    auto current_file = bconfig.config_filename_abs();
    /*if(solution_exists(stage_text) && !bconfig.control(RESOLVE)){
    if(bconfig.config_filename_abs() == current_file){
      if(rank == 0)
        std::cout << "Solved previously: " << bconfig.config_filename_abs() <<
  std::endl; return EXIT_SUCCESS; } else { return RELOAD_FILE;
    }
  }*/

    // get central values of the logarithmic enthalpy
    double loghc =
        ::Kadath::bco_utils::get_boundary_val(space.NS, logh, INNER_BC);

    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

    if (rank == 0 && solver_stage == TOTAL)
        std::cout << std::string(42, '#') << std::endl
                  << "TOTAL - Hydrostatic equilibrium stage\n"
                  << "using fixed lapse BC on the BH\n"
                  << std::string(42, '#') << std::endl;
    else if (rank == 0 && solver_stage == TOTAL_BC)
        std::cout << std::string(42, '#') << std::endl
                  << "TOTAL_BC - Hydrostatic equilibrium stage\n"
                  << "using v.Neumann lapse BC\n"
                  << std::string(42, '#') << std::endl;

    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);
    // setup a system of equations
    System_of_eqs syst(space, 0, ndom - 1);

    // the actual solution fields of the equations, i.e.
    // conformal factor, lapse, shift and (logarithmic) enthalpy
    // set this first before running syst_init()
    syst.add_var("H", logh);

    // populate all the boiler-plate constants, variables, and definitions
    syst_init(syst);

    // center of mass on the x-axis connecting both companions
    // and orbital angular frequency parameter
    // in this case, both are fixed by the two central force-balance
    // equations
    if (bconfig.control(FIXED_GOMEGA)) {
        syst.add_cst("ome", bconfig(GOMEGA));
        syst.add_cst("xaxis", bconfig(COM));
        syst.add_cst("yaxis", bconfig(COMY));
    } else {
        syst.add_var("ome", bconfig(GOMEGA));
        syst.add_var("xaxis", bconfig(COM));
        syst.add_var("yaxis", bconfig(COMY));
    }

    // central (logarithmic) enthalpy is a variable, fixed by the baryonic mass
    // integral
    syst.add_var("Hc1", loghc);

    // orbital rotation vector field, corrected by the "center of mass" shift
    syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");

    // full shift vector field, incoporating the orbital part
    syst.add_def("B^i= bet^i + ome * Morb^i");

    // the actual equations, defined differently in the different domains
    for (int d = 0; d < ndom; d++) {
        // if outside the stellar domains, without matter sources
        // resort to the source-free constraint equations
        // and set matter (and velocity potential) to zero
        if (d >= space.ADAPTEDNS + 1) {
            if (!bconfig.control(COROT_BIN))
                syst.add_eq_full(d, "phi= 0");

            syst.add_eq_full(d, "H  = 0");

            syst.add_def(d, "eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst.add_def(
                d,
                "eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
            syst.add_def(
                d,
                "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * "
                "D_j Ntilde");

        }
        // in case of the domains describing the NS interior
        // define all derived matter related quantities
        // and enforce the transformed contraint equations (i.e. multiplied by p /
        // rho)
        else {
            // 3-velocity and first integral of the Euler equation
            // in case of corotation
            if (bconfig.control(COROT_BIN)) {
                syst.add_def(d, "U^i    = B^i / N");
                syst.add_def(d, "Usquare= P^4 * U_i * U^i");
                syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
                syst.add_def(d, "W      = sqrt(Wsquare)");
                syst.add_def(d, "firstint = log(h * N / W)");
            }
            // 3-velocity, (approximate) first integral of the Euler equation
            // and velocity potential equation
            // in case of irrotational or spinning companions
            else {
                syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
                syst.add_def(d, "W      = sqrt(Wsquare)");
                syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
                syst.add_def(d, "Usquare= P^4 * U_i * U^i");
                syst.add_def(d, "V^i    = N * U^i - B^i");
                syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
                syst.add_def(d,
                             "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i "
                             "(P^6 * W * V^i)");
            }
            // transformed source terms
            syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
            syst.add_def(d,
                         "Stilde = 3 * press * delta + (Etilde + press * "
                         "delta) * Usquare");
            syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

            // transformed constraint equations
            syst.add_def(
                d,
                "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta "
                "+ 4piG / 2. * P^5 * Etilde");
            syst.add_def(
                d,
                "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * "
                "A_ij *A^ij "
                "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
            syst.add_def(
                d,
                "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
                "ptilde^i");

            // baryonic mass volume integrant
            syst.add_def(d, "intMb  = P^6 * rho * W");
            // quasi-local ADM mass volume integrant
            syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");
        }
    }
    // add equations to the system and demand continuity
    // along the normal of the domains
    space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq(syst, "eqP= 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

    // boundary conditions defining the boundary of the adapted domains, i.e.
    // vanishing matter
    syst.add_eq_bc(space.ADAPTEDNS, OUTER_BC, "H = 0");

    // in case of irrotational or spinning companions
    // solve also for the velocity potential
    if (!bconfig.control(COROT_BIN)) {
        syst.add_eq_bc(space.ADAPTEDNS, OUTER_BC, "V^i * D_i H = 0");

        for (int i = space.NS; i < space.ADAPTEDNS; ++i) {
            syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
            syst.add_eq_matching(i, OUTER_BC, "phi");
            syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
        }
        syst.add_eq_vel_pot(space.ADAPTEDNS, 2, "eqphi = 0", "phi=0");

        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
            space.add_eq_int_outer_NS(syst,
                                      "integ(intS1) / Madm1 / Madm1 = chi1");

        space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + s^i");
        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
            space.add_eq_int_BH(syst, "integ(intS2) - chi2 * Mch * Mch = 0 ");
    } else {
        space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");
    }

    // force-balance equations at the center of each star
    // to fix the orbital frequency as well as the "center of mass"
    // (the latter can get very inaccurate in extreme configurations
    // with large residual linear momenta, for which the next stage is given)
    Index posori1(space.get_domain(space.NS)->get_nbr_points());
    syst.add_eq_val(space.NS, "ex^i * D_i H", posori1);

    // add the first integral to get (approximate) hydrostatic equilibrium
    syst.add_eq_first_integral(space.NS,
                               space.ADAPTEDNS,
                               "firstint",
                               "H - Hc1");

    // fix the central enthalpy by baryonic mass volume integrals
    space.add_eq_int_volume(syst,
                            space.NS,
                            space.ADAPTEDNS,
                            "integvolume(intM) = qlMadm1");

    // compute a quasi-local approximation of the ADM component masses
    space.add_eq_int_volume(syst,
                            space.NS,
                            space.ADAPTEDNS,
                            "integvolume(intMb) = Mb1");

    if (!bconfig.control(FIXED_GOMEGA)) {
        space.add_eq_int_inf(syst, "integ(intPy) = 0");
        space.add_eq_int_inf(syst, "integ(intPx) = 0");
    }

    if (solver_stage == TOTAL) {
        // For TOTAL we use the fixed lapse condition
        space.add_bc_sphere_two(syst, "N = n0");
    } else {
        // For TOTAL_BC we use the von Neumann condition
        // and constrain Padm^x = 0
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

    // print initial diagnostics
    if (rank == 0)
        print_diagnostics(syst, 0, 0);

    // iterative solver variables
    bool endloop = false;
    int ite = 1;
    double conv;

    // loop until desired convergence is achieved
    while (!endloop) {
        // do exactly one Newton step
        endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

        // recompute all coordinate dependent fields
        // to make sure that they are updated correctly along
        // with the changing adapted domains
        update_fields(cfields, coord_vectors, {}, xo, xc1, xc2, &syst);

        // generate output filename for this iteration
        std::stringstream ss;
        ss << stage_text << ite - 1;
        bconfig.set_filename(ss.str());

        // print diagnostics and output configuration as well as the binary data
        if (rank == 0) {
            print_diagnostics(syst, ite, conv);
            if (bconfig.control(CHECKPOINT))
                checkpoint();
        }

        ite++;
        check_max_iter_exceeded(rank, ite, conv);
    }

    // since the ADM mass at infinite separation is not known
    // a priori for corotating stars,
    // we use the quasi-local measurement to update the ADM masses
    // approximately
    if (bconfig.control(COROT_BIN)) {
        bconfig.set(MADM, BCO1) = bconfig(QLMADM, BCO1);

        std::string ss = converged_filename(stage_text + "_COROT");
        bconfig.set_filename(ss);
    } else {
        std::string ss = converged_filename(stage_text);
        bconfig.set_filename(ss);
    }

    // output final configuration and binary data
    if (rank == 0)
        checkpoint();
    MPI_Barrier(MPI_COMM_WORLD);
    return exit_status;
}

template <class eos_t, typename config_t, typename space_t>
int bhns_xcts_solver<eos_t, config_t, space_t>::hydro_rescaling_stages(
    std::string stage_text) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_status = EXIT_SUCCESS;
    /*if(solution_exists(stage_text)) {
    if(rank == 0)
      std::cout << "Solved previously: " << bconfig.config_filename_abs() <<
  std::endl; return EXIT_SUCCESS;
  }*/

    if (rank == 0)
        if (solver_stage == ECC_RED)
            std::cout << std::string(42, '#') << std::endl
                      << "Eccentricity reduction step with fixed omega, chi = "
                      << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                      << " and q = " << bconfig(Q) << std::endl
                      << std::string(42, '#') << std::endl;
        else
            std::cout << std::string(42, '#') << std::endl
                      << "Hydro Rescaling stage \n"
                      << "using fixed orbital velocity with \n"
                      << "chi = " << bconfig(CHI, BCO1) << ","
                      << bconfig(CHI, BCO2) << " and q = " << bconfig(Q)
                      << std::endl
                      << std::string(42, '#') << std::endl;

    if (solver_stage == ECC_RED) {
        // determine whether to use PN estimates of the orbital frequency and adot
        // in case of the eccentricity stage
        if (std::isnan(bconfig.set(ADOT)) ||
            std::isnan(bconfig.set(ECC_OMEGA)) || bconfig.control(USE_PN)) {
            ::Kadath::bco_utils::KadathPNOrbitalParams(bconfig,
                                                       bconfig(MADM, BCO1),
                                                       bconfig(MCH, BCO2));

            if (rank == 0)
                std::cout << "### Using PN estimate for adot and omega! ###"
                          << std::endl;
            bconfig.set(ECC_OMEGA) = bconfig(GOMEGA);
        } else {
            bconfig.set(GOMEGA) = bconfig(ECC_OMEGA);
        }
    }

    // setup background position vector field - only needed for ECC_RED stage
    Vector CART(space, CON, basis);
    CART = cfields.cart();

    // residual scaling factors for the ethalpy
    // used to correct the baryonic mass
    double H_scale = 0;

    // set the shift of the "center of mass" along the y-axis
    // to zero, potentially fixing x-components of the ADM linear
    // momentum at infinity introduced by the eccentricity reduction
    // parameters
    if (std::isnan(bconfig.set(COMY)))
        bconfig.set(COMY) = 0.;

    // "background", unscaled logarithmic enthalpy from the previous step
    Scalar logh_const(logh);
    logh_const.std_base();

    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

    // setup a system of equations
    System_of_eqs syst(space, 0, ndom - 1);

    // constant "background" enthalpy
    syst.add_cst("Hconst", logh_const);

    // residual scaling factors for the ethalpy
    syst.add_var("Hscale", H_scale);

    // definition of H
    // - rescaled in the star
    // - constant everywhere else
    for (int d = 0; d < ndom; ++d) {
        if (d <= space.ADAPTEDNS && d >= space.NS)
            syst.add_def(d, "H  = Hconst * (1. + Hscale)");
        else
            syst.add_def(d, "H  = Hconst");
    }

    // populate all the boiler-plate constants, variables, and definitions
    syst_init(syst);

    bool fixed_com = this->solver_stage != ECC_RED;
    if (fixed_com) {
        // "center of mass" on the x-axis, connecting both stellar centers
        // These are fixed on the initial ID import as integrals at INF
        // cause instabilities in the initial solution
        syst.add_cst("xaxis", bconfig(COM));
        syst.add_cst("yaxis", bconfig(COMY));
    } else {
        // "center of mass" on the x-axis, connecting both stellar centers
        // fixed by the vanishing of the ADM linear momentum at infinity
        syst.add_var("xaxis", bconfig(COM));

        // same on the y-axis in case finite momenta develope by
        // the eccentricity reduction parameters
        syst.add_var("yaxis", bconfig(COMY));
    }

    // no additional force-balance is computed,
    // the matter distribution is fixed modulo the scaling factos above,
    // therefore the orbital frequency is a constant similar to
    // eccentricity reduced ID
    syst.add_cst("ome", bconfig(GOMEGA));

    // orbital rotation vector field, corrected by the "center of mass" shift
    syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");

    std::string bigB{"B^i = bet^i + ome * Morb^i"};
    if (solver_stage == ECC_RED) {
        syst.add_cst("adot", bconfig(ADOT));
        syst.add_cst("r", CART);
        syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
        // add contribution of ADOT to the total shift definition
        bigB += " + adot * comr^i";
    }

    // full shift vector including inertial + orbital contributions
    syst.add_def(bigB.c_str());

    // the actual equations, defined differently in the different domains
    for (int d = 0; d < ndom; d++) {
        // if outside the stellar domains, without matter sources
        // resort to the source-free constraint equations
        // and set matter (and velocity potential) to zero
        if (d >= space.ADAPTEDNS + 1) {
            if (!bconfig.control(COROT_BIN))
                syst.add_eq_full(d, "phi= 0");

            syst.add_def(d, "eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst.add_def(
                d,
                "eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
            syst.add_def(
                d,
                "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * "
                "D_j Ntilde");

        }
        // in case of the domains harboring the two stars
        // define all derived matter related quantities
        // and enforce the transformed contraint equations (i.e. multiplied by p /
        // rho)
        else {
            // 3-velocity and first integral of the Euler equation
            // in case of corotation
            if (bconfig.control(COROT_BIN)) {
                syst.add_def(d, "U^i    = B^i / N");
                syst.add_def(d, "Usquare= P^4 * U_i * U^i");
                syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
                syst.add_def(d, "W      = sqrt(Wsquare)");
                syst.add_def(d, "firstint = log(h * N / W)");
            }
            // 3-velocity, (approximate) first integral of the Euler equation
            // and velocity potential equation
            // in case of irrotational or spinning companions
            else {
                syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
                syst.add_def(d, "W      = sqrt(Wsquare)");
                syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
                syst.add_def(d, "Usquare= P^4 * U_i * U^i");
                syst.add_def(d, "V^i    = N * U^i - B^i");
                syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
                syst.add_def(d,
                             "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i "
                             "(P^6 * W * V^i)");
            }
            // transformed source terms
            syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
            syst.add_def(d,
                         "Stilde = 3 * press * delta + (Etilde + press * "
                         "delta) * Usquare");
            syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

            // transformed constraint equations
            syst.add_def(
                d,
                "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta "
                "+ 4piG / 2. * P^5 * Etilde");
            syst.add_def(
                d,
                "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * "
                "A_ij *A^ij "
                "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
            syst.add_def(
                d,
                "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
                "ptilde^i");

            // baryonic mass volume integrant
            syst.add_def(d, "intMb  = P^6 * rho * W");
            // quasi-local ADM mass volume integrant
            syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");
        }
    }
    // add equations to the system and demand continuity
    // along the normal of the domains
    space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq(syst, "eqP= 0", "P", "dn(P)");
    space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
    syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

    // boundary conditions defining the boundary of the adapted domains, i.e.
    // vanishing matter
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
            space.add_eq_int_outer_NS(syst,
                                      "integ(intS1) / Madm1 / Madm1 = chi1");

        space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i + s^i");
        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
            space.add_eq_int_BH(syst, "integ(intS2) - chi2 * Mch * Mch = 0 ");
    } else {
        space.add_bc_sphere_two(syst, "B^i = N / P^2 * sp^i");
    }

    space.add_eq_int_volume(syst,
                            space.NS,
                            space.ADAPTEDNS,
                            "integvolume(intM) = qlMadm1");
    space.add_eq_int_volume(syst,
                            space.NS,
                            space.ADAPTEDNS,
                            "integvolume(intMb) = Mb1");

    // stage == ECC_RED
    if (!fixed_com) {
        space.add_eq_int_inf(syst, "integ(intPx) = 0");
        space.add_eq_int_inf(syst, "integ(intPy) = 0");
    }

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

    // print initial diagnostics
    if (rank == 0)
        print_diagnostics(syst, 0, 0);

    // iterative solver variables
    bool endloop = false;
    int ite = 1;
    double conv;

    // loop until desired convergence is achieved
    while (!endloop) {
        // do exactly one Newton step
        endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

        // recompute all coordinate dependent fields
        // to make sure that they are updated correctly along
        // with the changing adapted domains
        update_fields(cfields, coord_vectors, {}, xo, xc1, xc2, &syst);

        // generate output filename for this iteration
        std::stringstream ss;
        ss << "total_" << ite - 1;
        bconfig.set_filename(ss.str());

        // overwrite logh with scaled version for later use and output
        logh = syst.give_val_def("H");

        // print diagnostics and output configuration as well as the binary data
        if (rank == 0) {
            print_diagnostics(syst, ite, conv);
            if (bconfig.control(CHECKPOINT))
                checkpoint();
        }

        ite++;
        check_max_iter_exceeded(rank, ite, conv);
    }

    // since the ADM mass at infinite separation is not known
    // a priori for corotating stars,
    // we use the quasi-local measurement to update the ADM masses
    // approximately
    if (bconfig.control(COROT_BIN)) {
        bconfig.set(MADM, BCO1) = bconfig(QLMADM, BCO1);

        std::string ss = converged_filename(stage_text + "_COROT");
        bconfig.set_filename(ss);
    } else {
        std::string ss = converged_filename(stage_text);
        bconfig.set_filename(ss);
    }

    // output final configuration and binary data
    if (rank == 0)
        checkpoint();
    return exit_status;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath