#include "FUKA_Solvers/utilities/format_settings.hpp"
#include "Solvers/fuka_syst/fuka_syst_setup.hpp"
#include "utilities/solver_utilities.hpp"

namespace Kadath::FUKA_Solvers {
// NOROT Routines
template <class eos_t>
NS_XCTS_UNIFORM_ROT<eos_t>::NS_XCTS_UNIFORM_ROT(
    NS_XCTS_BASE::base_config_t* config_,
    ns_sequence const& seq_,
    Parameter_sequence<BCO_PARAMS> const& res_,
    std::string outputdir_,
    int const rank_)
    : NS_XCTS_BASE(config_, seq_, res_, outputdir_, rank_), spinup(nullptr) {
    stagename = "UNIFORM_ROT";
    solver_stage = ::Kadath::FUKA_Config::STAGES::UNIFORM_ROT;
    if (!seq->is_set() && !bconfig->control(CONTROLS::SEQUENCES)) {
        initialize_config_from_fixing_values(*bconfig, *seq);
    }
    load_solution_from_file();
    initialize_EOS(*this);
    initialize_support_containers();
    initialize_spinup();

    if (rank == 0) {
        cout << *seq << endl;
        cout << *resolution << endl;
    }
}

template <class eos_t>
std::string NS_XCTS_UNIFORM_ROT<eos_t>::converged_filename(
    const std::string stage) const {
    // FIXME assumes a fix resolution for all domains
    auto res = space->get_domain(0)->get_nbr_points()(0);
    const std::string eosname{extract_eos_name(*bconfig)};
    std::stringstream ss;
    ss << "NS";
    if (stage != "")
        ss << "_" << stage << ".";
    else
        ss << ".";
    ss << eosname << ".";

    // Add mass fixing parameter to filename
    auto default_idx = BCO_PARAMS::MADM;
    if (seq) {
        update_filename_from_mass_fixing(*bconfig, seq, ss);
    } else {
        auto [seq_key, tidx] =
            get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
        ss << seq_key << "." << (*bconfig)(default_idx) << ".";
    }

    default_idx = BCO_PARAMS::CHI;
    if (seq) {
        update_filename_from_spin_fixing(*bconfig, seq, ss);
    } else {
        auto [seq_key, tidx] =
            get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
        ss << seq_key << "." << (*bconfig)(default_idx) << ".";
    }
    ss << (*bconfig)(BCO_PARAMS::NSHELLS) << "." << std::setfill('0')
       << std::setw(2) << res;
    return ss.str();
}

template <class eos_t>
void NS_XCTS_UNIFORM_ROT<eos_t>::setup_syst() {
    int exit_status = EXIT_SUCCESS;
    double loghc = std::log((*bconfig)(BCO_PARAMS::HC));

    // We use `config_filename()` vs `config_filename_abs()` since
    // `solution_exists` will probe the HOME_KADATH/COs directory
    // auto const current = bconfig.config_filename();
    // if(!bconfig.control(RESOLVE) && solution_exists(stagename)) {
    //   if(rank == 0)
    //     std::cout << "Solved previously: " \
    //               << bconfig.config_filename_abs() << std::endl;
    //   return (current == bconfig.config_filename()) ? \
    //     EXIT_SUCCESS : RELOAD_FILE;
    // }

    update_fields_co(*cfields, *coord_vectors, {}, 0.);
    syst.reset(new System_of_eqs(*space));
    syst_init();

    std::string central_fixing_definition{"h - hc"};
    std::string spin_fixing_definition{"integ(intJ) - chi * Madm * Madm = 0"};
    std::string output_str_M{};
    std::string output_str_Spin{};

    if (seq) {
        central_fixing_definition =
            ::Kadath::FUKA_Syst_tools::set_ns_mass_fixing(*syst, *bconfig, seq);
        spin_fixing_definition =
            ::Kadath::FUKA_Syst_tools::set_ns_spin_fixing(*syst, *bconfig, seq);

        output_str_M =
            ::Kadath::FUKA_Syst_tools::get_ns_mass_fixing_output(*bconfig, seq);
        output_str_Spin =
            ::Kadath::FUKA_Syst_tools::get_ns_spin_fixing_output(*bconfig, seq);
    } else {
        syst->add_var("Hc", (*bconfig)(BCO_PARAMS::HC));
        syst->add_cst("chi", (*bconfig)(BCO_PARAMS::CHI));
        syst->add_var("ome", (*bconfig)(BCO_PARAMS::OMEGA));
        syst->add_cst("Madm", (*bconfig)(BCO_PARAMS::MADM));
        output_str_M = "Mass fixed using ADM mass (madm)";
        output_str_Spin = "Spin fixed using dimensionless spin parameter (chi)";
    }
    syst->add_def("omega^i = bet^i + ome * mg^i");

    for (int d = 0; d < ndom; d++) {
        if (d <= 1) {
            syst->add_def(d, "U^i = omega^i / N");
            syst->add_def(d, "Usquare = P^4 * U_i * U^i");
            syst->add_def(d, "Wsquare = 1. / (1. - Usquare)");
            syst->add_def(d, "W = sqrt(Wsquare)");

            syst->add_def(d, "Etilde = press * h * Wsquare - press * delta");
            syst->add_def(d,
                          "Stilde = 3 * press * delta + (Etilde + press * "
                          "delta) * Usquare");
            syst->add_def(d, "ptilde^i = press * h * Wsquare * U^i");

            syst->add_def(
                d,
                "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * "
                "delta + 4piG / "
                "2. * P^5 * Etilde");
            syst->add_def(d,
                          "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * "
                          "delta * A_ij *A^ij "
                          "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
            syst->add_def(
                d,
                "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
                "ptilde^i");

            syst->add_def(d, "intMb = P^6 * rho(h) * W");
            syst->add_def(d, "firstint = H + log(N) - log(W)");

        } else {
            syst->add_eq_full(d, "H = 0");

            syst->add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
            syst->add_def(
                d,
                "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
            syst->add_def(d,
                          "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                          "A^ij * D_j Ntilde");
            break;
        }
    }

    if (rank == 0) {
        std::cout << "############################" << std::endl
                  << "Uniformly Rotating NS Solver" << std::endl
                  << output_str_M << std::endl
                  << output_str_Spin << std::endl
                  << "############################" << std::endl;
    }

    // add the constraint equations and demand continuity their normal derivative
    // across domain boundaries
    space->add_eq(*syst, "eqNP= 0", "N", "dn(N)");
    space->add_eq(*syst, "eqP = 0", "P", "dn(P)");
    space->add_eq(*syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    syst->add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst->add_eq_bc(ndom - 1, OUTER_BC, "P=1");
    syst->add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

    // if the surface is resolved, define it to be where the matter vanishes
    syst->add_eq_bc(1, OUTER_BC, "H = 0");

    // first integral in the innermost domains with non-zero matter content
    // and condition on the central value, either fixed directly or by the
    // integral below
    syst->add_eq_first_integral(0,
                                1,
                                "firstint",
                                central_fixing_definition.c_str());

    // constrain stellar mass by these integrals and the central log enthalpy
    if (seq) {
        auto idx{seq->mass_idx()};
        bool add_Madm_int = true;
        switch (idx) {
            case BCO_PARAMS::MADM:
                space->add_eq_int_inf(*syst, "integ(intMadm) = Madm");
                add_Madm_int = false;
                break;
            case BCO_PARAMS::MB:
                space->add_eq_int_volume(*syst, 2, "integvolume(intMb) = Mb");
                break;
            default:
                break;
        }
        idx = seq->spin_idx();
        switch (idx) {
            case BCO_PARAMS::JADM:
                space->add_eq_int_inf(*syst, spin_fixing_definition.c_str());
                break;
            case BCO_PARAMS::CHI:
                // Since we need MADM to compute CHI, we need to ensure
                // that if it isn't a fixed quantity that it becomes a
                // variable in our system of equations and the appropriate
                // constraint equation is added
                if (add_Madm_int) {
                    syst->add_var("Madm", (*bconfig)(BCO_PARAMS::MADM));
                    space->add_eq_int_inf(*syst, "integ(intMadm) = Madm");
                }

                space->add_eq_int_inf(*syst, spin_fixing_definition.c_str());
                break;
            default:
                break;
        }
    } else {
        space->add_eq_int_inf(*syst, spin_fixing_definition.c_str());
        space->add_eq_int_inf(*syst, "integ(intMadm) = Madm");
    }
}

template <class eos_t>
void NS_XCTS_UNIFORM_ROT<eos_t>::syst_init() {
    using namespace ::Kadath::Margherita;

    // call the (flat) conformal metric "f"
    fmet->set_system(*syst, "f");

    // define numerical constants
    syst->add_cst("4piG", (*bconfig)(BCO_QPIG));

    // include the coordinate fields
    syst->add_cst("mg", *(*coord_vectors)[GLOBAL_ROT]);
    syst->add_cst("sm", *(*coord_vectors)[S_BCO1]);
    syst->add_cst("einf", *(*coord_vectors)[S_INF]);

    // the basic fields, conformal factor, lapse and (log) enthalpy
    syst->add_var("P", *conformal_factor);
    syst->add_var("N", *lapse);
    syst->add_var("H", *logh);
    syst->add_var("bet", *shift);

    // define common combinations of conformal factor and lapse
    syst->add_def("NP = P*N");
    syst->add_def("Ntilde = N / P^6");
    syst->add_def(
        "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / "
        "2. / Ntilde");

    // define quantity to be integrated at infinity
    // two (in this case) equivalent definitions of ADM mass
    // as well as the Komar mass
    syst->add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
    syst->add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
    syst->add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");

    // ADM Angular momentum
    syst->add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
    // Quasi-local spin angular momentum
    syst->add_def(2, "intS = A_ij * mg^i * sm^j / 2. / 4piG");

    // enthalpy from the logarithmic enthalpy, the latter is the actual variable
    // in this system
    syst->add_def("h = exp(H)");

    // define the EOS operators
    Param p;
    set_eos_ope<eos_t>(*syst, p);

    // define rest-mass density, internal energy and pressure through the
    // enthalpy
    syst->add_def("rho = rho(h)");
    syst->add_def("eps = eps(h)");
    syst->add_def("press = press(h)");
    syst->add_def("dHdlnrho = dHdlnrho(h)");

    // definition to rescale the equations
    // delta = p / rho
    syst->add_def("delta = h - eps - 1.");
}

template <class eos_t>
void NS_XCTS_UNIFORM_ROT<eos_t>::print_diagnostics(const int ite,
                                                   const double conv) const {
    // compute the baryonic mass at volume integral from the given integrant
    double baryonic_mass = syst->give_val_def("intMb")()(0).integ_volume() +
                           syst->give_val_def("intMb")()(1).integ_volume();

    // compute the ADM mass as surface integral at infinity
    Val_domain integMadm(syst->give_val_def("intMadm")()(ndom - 1));
    double Madm = space->get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

    // compute the Komar mass as surface integral at infinity
    Val_domain integMk(syst->give_val_def("intMk")()(ndom - 1));
    double Mk = space->get_domain(ndom - 1)->integ(integMk, OUTER_BC);

    // get the maximum and minimum coordinate radius along the surface,
    // i.e. the adapted domain boundary
    auto rs = Kadath::bco_utils::get_rmin_rmax(*space, 1);

    // alternative, equivalent ADM mass integral
    Val_domain integMadmalt(syst->give_val_def("intMadmalt")()(ndom - 1));
    double Madmalt = space->get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

    // compute the ADM angular momentum as surface integral at infinity
    Val_domain integJ(syst->give_val_def("intJ")()(ndom - 1));
    double J = space->get_domain(ndom - 1)->integ(integJ, OUTER_BC);

    std::ios_base::fmtflags f(std::cout.flags());
    std::cout << "=======================================" << std::endl
              << FORMAT << "Iter: " << ite << std::endl
              << FORMAT << "Error: " << conv << std::endl
              << FORMAT << "Mb: " << baryonic_mass << std::endl
              << FORMAT << "Madm: " << Madm << std::endl
              << FORMAT << "Madm_ql: " << Madmalt << " ["
              << std::abs(Madm - Madmalt) / Madm << "]" << std::endl
              << FORMAT << "Mk: " << Mk << " [" << std::abs(Madm - Mk) / Madm
              << "]" << std::endl
              << FORMAT << "R: " << rs[0] << " " << rs[1] << std::endl;
    std::cout << FORMAT << "Jadm: " << J << std::endl
              << FORMAT << "Chi: " << J / Madm / Madm << " [" << (*bconfig)(CHI)
              << "]\n"
              << FORMAT << "Omega: " << (*bconfig)(OMEGA) << std::endl;
    std::cout.flags(f);
#undef FORMAT
    std::cout << "=======================================" << "\n\n";
}  // end print diagnostics norot

template <class eos_t>
void NS_XCTS_UNIFORM_ROT<eos_t>::update_config_quantities() {
    auto rs = bco_utils::get_rmin_rmax(*space, 1);
    bconfig->set(BCO_PARAMS::RMID) = rs[0];

    // compute the ADM mass as surface integral at infinity
    Val_domain integMadm(syst->give_val_def("intMadm")()(ndom - 1));
    double Madm = space->get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

    // compute the baryonic mass at volume integral from the given integrant
    double baryonic_mass = syst->give_val_def("intMb")()(0).integ_volume() +
                           syst->give_val_def("intMb")()(1).integ_volume();

    auto loghc = bco_utils::get_boundary_val(0, *logh, INNER_BC);

    // compute the ADM Angular Momentum as surface integral at infinity
    Val_domain integJ(syst->give_val_def("intJ")()(ndom - 1));
    double const Jadm = space->get_domain(ndom - 1)->integ(integJ, OUTER_BC);
    double const chi = Jadm / Madm / Madm;

    if (seq) {
        auto idx{seq->mass_idx()};
        switch (idx) {
            case BCO_PARAMS::HC:
                bconfig->set(BCO_PARAMS::NC) =
                    EOS<eos_t, DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
                bconfig->set(BCO_PARAMS::MADM) = Madm;
                bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
                break;
            case BCO_PARAMS::NC:
                bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
                bconfig->set(BCO_PARAMS::MADM) = Madm;
                bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
                break;
            case BCO_PARAMS::MB:
                bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
                bconfig->set(BCO_PARAMS::NC) =
                    EOS<eos_t, DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
                bconfig->set(BCO_PARAMS::MADM) = Madm;
                break;
            default:
                bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
                bconfig->set(BCO_PARAMS::NC) =
                    EOS<eos_t, DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
                bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
                break;
        }

        idx = seq->spin_idx();
        switch (idx) {
            case BCO_PARAMS::CHI:
                break;
            default:
                bconfig->set(BCO_PARAMS::CHI) = chi;
                break;
        }
    } else {
        bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
        bconfig->set(BCO_PARAMS::NC) =
            EOS<eos_t, DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
        bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
        bconfig->set(BCO_PARAMS::CHI) = chi;
    }
    bconfig->set(BCO_PARAMS::QLMADM) = bconfig->set(BCO_PARAMS::MADM);
}

template <class eos_t>
void NS_XCTS_UNIFORM_ROT<eos_t>::initialize_spinup() {
    auto const spinidx = seq->spin_idx();
    auto finalspin = (*bconfig)(spinidx);

    // No need to spinup if we're computing a sequence of
    // spinning NS starting from ~zero
    if ((ns_seq_is_spin_fixing(*seq) && std::fabs(seq->init()) < 1e-2) ||
        ((seq->spin_idx() == BCO_PARAMS::CHI) &&
         (std::fabs(seq->spin_val()) < 0.2))) {
        return;
    } else if (ns_seq_is_spin_fixing(*seq)) {
        finalspin = seq->init();
    }

    auto npts = space->get_domain(1)->get_nbr_points();
    Index pos_eq(npts);
    pos_eq.set(0) = npts(0) - 1;  /// Set to outer radius
    pos_eq.set(1) = npts(1) - 1;  /// Set theta to be on the xy plane.

    Index pos_pole(npts);
    pos_pole.set(0) = npts(0) - 1;  /// Set to outer radius

    auto adpt_dom = space->get_domain(1);
    double const R0 = adpt_dom->get_radius()(pos_eq);
    double const Rp = adpt_dom->get_radius()(pos_pole);
    double const axis_ratio = Rp / R0;

    // no reason to spinup if the solution is already sufficiently rotating
    if (1. - axis_ratio > 1e-3) {
        return;
    }
    spinup.reset(
        new Parameter_sequence<BCO_PARAMS>(seq->spin_str(), seq->spin_idx()));
    spinup->set(0., 0., finalspin);
    spinup->set_N(3);
    bconfig->set(spinidx) = 0.;
}

template <class eos_t>
bool NS_XCTS_UNIFORM_ROT<eos_t>::increment_spin() {
    if (!spinup || !spinup->is_set() || !spinup->is_varying()) {
        return false;
    }
    auto sequence_var_indices = spinup->get_indices();
    auto const& dx = spinup->step_size();
    auto x = bconfig->set(sequence_var_indices);
    if (spinup->loop_condition(x)) {
        x += dx;
        x = (std::fabs(1. - x / spinup->final()) < 1e-4) ? spinup->final() : x;
        bconfig->set(sequence_var_indices) = x;
        return true;
    }
    return false;
}
}  // namespace Kadath::FUKA_Solvers
