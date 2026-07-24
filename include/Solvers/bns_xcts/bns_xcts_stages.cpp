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
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
namespace bco_u = ::Kadath::bco_utils;

template <class eos_t, typename config_t, typename space_t>
int bns_xcts_solver<eos_t, config_t, space_t>::hydrostatic_equilibrium_stage() {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_status = EXIT_SUCCESS;

    // get central values of the logarithmic enthalpy
    double loghc1 = bco_u::get_boundary_val(space.NS1, logh, INNER_BC);
    double loghc2 = bco_u::get_boundary_val(space.NS2, logh, INNER_BC);

    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

    // first stage, baryonic mass fixing but without explicit corrections to the
    // linear momentum at infinty
    if (rank == 0)
        std::cout << "############################" << std::endl
                  << "Total with mass fixing, chi = " << bconfig(CHI, BCO1)
                  << "," << bconfig(CHI, BCO2) << " and q = " << bconfig(Q)
                  << std::endl
                  << "############################" << std::endl;

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
    syst.add_var("xaxis", bconfig(COM));
    syst.add_var("ome", bconfig(GOMEGA));

    // central (logarithmic) enthalpy is a variable, fixed by the baryonic mass
    // integral
    syst.add_var("Hc1", loghc1);
    syst.add_var("Hc2", loghc2);

    // orbital rotation vector field, corrected by the "center of mass" shift
    syst.add_def("Morb^i = mg^i + xaxis * ey^i");

    // full shift vector field, incoporating the orbital part
    syst.add_def("B^i= bet^i + ome * Morb^i");
    syst_init_Aterms(syst);

    // the actual equations, defined differently in the different domains
    for (int d = 0; d < ndom; d++) {
        // if outside the stellar domains, without matter sources
        // resort to the source-free constraint equations
        // and set matter (and velocity potential) to zero
        if ((d >= space.ADAPTED2 + 1) ||
            ((d >= space.ADAPTED1 + 1) && (d < space.NS2))) {
            if (!bconfig.control(COROT_BIN))
                syst.add_eq_full(d, "phi= 0");

            syst.add_eq_full(d, "H  = 0");

            syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst.add_def(
                d,
                "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
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
    syst.add_eq_bc(space.ADAPTED1, OUTER_BC, "H = 0");
    syst.add_eq_bc(space.ADAPTED2, OUTER_BC, "H = 0");

    // in case of irrotational or spinning companions
    // solve also for the velocity potential
    if (!bconfig.control(COROT_BIN)) {
        // boundary conditions at the stellar surface,
        // i.e. where the velocity potential equation becomes
        // degenerate
        syst.add_eq_bc(space.ADAPTED1, OUTER_BC, "V^i * D_i H = 0");
        syst.add_eq_bc(space.ADAPTED2, OUTER_BC, "V^i * D_i H = 0");

        // add the equation and the matchings to the system
        // in case of the stellar domains
        for (int i = space.NS1; i < space.ADAPTED1; ++i) {
            syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
            syst.add_eq_matching(i, OUTER_BC, "phi");
            syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
        }
        syst.add_eq_vel_pot(space.ADAPTED1, 2, "eqphi = 0", "phi=0");

        for (int i = space.NS2; i < space.ADAPTED2; ++i) {
            syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
            syst.add_eq_matching(i, OUTER_BC, "phi");
            syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
        }
        syst.add_eq_vel_pot(space.ADAPTED2, 2, "eqphi = 0", "phi=0");

        // in case the frequency parameter of the star is not fixed by hand
        // solve for the given dimensionless spin through the quasi-local
        // spin surface integral
        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
            space.add_eq_int_outer_sphere_one(
                syst,
                "integ(intS1) / Madm1 / Madm1 = chi1");

        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
            space.add_eq_int_outer_sphere_two(
                syst,
                "integ(intS2) / Madm2 / Madm2 = chi2");
    }

    // force-balance equations at the center of each star
    // to fix the orbital frequency as well as the "center of mass"
    // (the latter can get very inaccurate in extreme configurations
    // with large residual linear momenta, for which the next stage is given)
    Index posori1(space.get_domain(space.NS1)->get_nbr_points());
    syst.add_eq_val(space.NS1, "ex^i * D_i H", posori1);

    Index posori2(space.get_domain(space.NS2)->get_nbr_points());
    syst.add_eq_val(space.NS2, "ex^i * D_i H", posori2);

    // add the first integral to get (approximate) hydrostatic equilibrium
    syst.add_eq_first_integral(space.NS1,
                               space.ADAPTED1,
                               "firstint",
                               "H - Hc1");
    syst.add_eq_first_integral(space.NS2,
                               space.ADAPTED2,
                               "firstint",
                               "H - Hc2");

    // fix the central enthalpy by baryonic mass volume integrals
    space.add_eq_int_volume(syst,
                            space.NS1,
                            space.ADAPTED1,
                            "integvolume(intMb) = Mb1");
    space.add_eq_int_volume(syst,
                            space.NS2,
                            space.ADAPTED2,
                            "integvolume(intMb) = Mb2");

    // compute a quasi-local approximation of the ADM component masses
    // space.add_eq_int_volume(syst, space.NS1, space.ADAPTED1,
    //                         "integvolume(intM) = qlMadm1");
    // space.add_eq_int_volume(syst, space.NS2, space.ADAPTED2,
    //                         "integvolume(intM) = qlMadm2");

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

        update_config_quantities(syst);

        // generate output filename for this iteration
        std::stringstream ss;
        ss << "total_" << ite - 1;
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
        bconfig.set(MADM, BCO2) = bconfig(QLMADM, BCO2);

        std::string ss = converged_filename("TOTAL_COROT");
        bconfig.set_filename(ss);
    } else {
        std::string ss = converged_filename("TOTAL");
        bconfig.set_filename(ss);
    }

    // output final configuration and binary data
    if (rank == 0)
        checkpoint();
    return exit_status;
}

template <class eos_t, typename config_t, typename space_t>
int bns_xcts_solver<eos_t, config_t, space_t>::hydro_rescaling_stages(
    std::string stage_text) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int exit_status = EXIT_SUCCESS;

    if (rank == 0)
        if (solver_stage == TOTAL_BC) {
            std::cout
                << "############################" << std::endl
                << "Py fixing without first integral and fixed omega, chi = "
                << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                << " and q = " << bconfig(Q) << std::endl
                << "with Py BC fixing" << std::endl
                << "############################" << std::endl;
        } else if (solver_stage == ECC_RED) {
            std::cout << "############################" << std::endl
                      << "Eccentricity reduction step with fixed omega, chi = "
                      << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                      << " and q = " << bconfig(Q) << std::endl
                      << "############################" << std::endl;
        } else if (solver_stage == HEADON) {
            std::cout << "############################" << std::endl
                      << "Head-on collision of two stars from rest."
                      << "q = " << bconfig(Q) << std::endl
                      << "############################" << std::endl;
        }

    if (solver_stage == ECC_RED) {
        // determine whether to use PN estimates of the orbital frequency and adot
        // in case of the eccentricity stage
        if (std::isnan(bconfig.set(ADOT)) ||
            std::isnan(bconfig.set(ECC_OMEGA)) || bconfig.control(USE_PN)) {
            bco_u::KadathPNOrbitalParams(bconfig,
                                         bconfig(MADM, BCO1),
                                         bconfig(MADM, BCO2));

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
    double H_scale1 = 0;
    double H_scale2 = 0;

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
    syst.add_var("Hscale1", H_scale1);
    syst.add_var("Hscale2", H_scale2);

    for (int d = 0; d < ndom; d++)
        // the enthalpy is equal to the constant part everywhere
        // outside of the stars, i.e. zero
        if ((d >= space.ADAPTED2 + 1) ||
            ((d >= space.ADAPTED1 + 1) && (d < space.NS2)))
            syst.add_def(d, "H  = Hconst");

    // inside the stars, it is the constant part
    // and a small correction factor to fix the correct
    // baryonic mass, due to slightly changing
    // fluid velocities
    for (int d = space.NS1; d <= space.ADAPTED1; ++d) {
        syst.add_def(d, "H  = Hconst * (1. + Hscale1)");
    }
    for (int d = space.NS2; d <= space.ADAPTED2; ++d) {
        syst.add_def(d, "H  = Hconst * (1. + Hscale2)");
    }

    // populate all the boiler-plate constants, variables, and definitions
    syst_init(syst);

    if (solver_stage == HEADON) {
        // "center of mass" on the x-axis, connecting both stellar centers
        // fixed by the vanishing of the ADM linear momentum at infinity
        syst.add_cst("xaxis", bconfig(COM));

        // same on the y-axis in case finite momenta develope by
        // the eccentricity reduction parameters
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
    syst_init_Aterms(syst);

    // the actual equations, defined differently in the different domains
    for (int d = 0; d < ndom; d++) {
        // if outside the stellar domains, without matter sources
        // resort to the source-free constraint equations
        // and set matter (and velocity potential) to zero
        if ((d >= space.ADAPTED2 + 1) ||
            ((d >= space.ADAPTED1 + 1) && (d < space.NS2))) {
            if (!bconfig.control(COROT_BIN))
                syst.add_eq_full(d, "phi= 0");

            syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst.add_def(
                d,
                "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
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
    syst.add_eq_bc(space.ADAPTED1, OUTER_BC, "H = 0");
    syst.add_eq_bc(space.ADAPTED2, OUTER_BC, "H = 0");

    // in case of irrotational or spinning companions
    // solve also for the velocity potential
    if (!bconfig.control(COROT_BIN)) {
        // boundary conditions at the stellar surface,
        // i.e. where the velocity potential equation becomes
        // degenerate
        syst.add_eq_bc(space.ADAPTED1, OUTER_BC, "V^i * D_i H = 0");
        syst.add_eq_bc(space.ADAPTED2, OUTER_BC, "V^i * D_i H = 0");

        // add the equation and the matchings to the system
        // in case of the stellar domains
        for (int i = space.NS1; i < space.ADAPTED1; ++i) {
            syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
            syst.add_eq_matching(i, OUTER_BC, "phi");
            syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
        }
        syst.add_eq_vel_pot(space.ADAPTED1, 2, "eqphi = 0", "phi=0");

        for (int i = space.NS2; i < space.ADAPTED2; ++i) {
            syst.add_eq_vel_pot(i, 2, "eqphi = 0", "phi=0");
            syst.add_eq_matching(i, OUTER_BC, "phi");
            syst.add_eq_matching(i, OUTER_BC, "dn(phi)");
        }
        syst.add_eq_vel_pot(space.ADAPTED2, 2, "eqphi = 0", "phi=0");

        // in case the frequency parameter of the star is not fixed by hand
        // solve for the given dimensionless spin through the quasi-local
        // spin surface integral
        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
            space.add_eq_int_outer_sphere_one(
                syst,
                "integ(intS1) / Madm1 / Madm1 = chi1");

        if (std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
            space.add_eq_int_outer_sphere_two(
                syst,
                "integ(intS2) / Madm2 / Madm2 = chi2");
    }

    if (solver_stage != HEADON) {
        // enforcing the vanishing of the x- & y-component of the ADM linear momentum
        // at infinity, fixing the "center of mass" shift on the x-axis
        space.add_eq_int_inf(syst, "integ(intPx) = 0");
        space.add_eq_int_inf(syst, "integ(intPy) = 0");
    }

    // fix the central enthalpy by baryonic mass volume integrals
    space.add_eq_int_volume(syst,
                            space.NS1,
                            space.ADAPTED1,
                            "integvolume(intMb) = Mb1");
    space.add_eq_int_volume(syst,
                            space.NS2,
                            space.ADAPTED2,
                            "integvolume(intMb) = Mb2");

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

        update_config_quantities(syst);

        // overwrite logh with scaled version for later use and output
        logh = syst.give_val_def("H");

        // generate output filename for this iteration
        std::stringstream ss;
        ss << "total_" << ite - 1;
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
        bconfig.set(MADM, BCO2) = bconfig(QLMADM, BCO2);

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