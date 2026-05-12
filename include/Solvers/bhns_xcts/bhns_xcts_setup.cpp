/**
 * \addtogroup BHNS_XCTS
 * \ingroup FUKA
 * @{*/
#include "EOS/FUKA_EOS_Utilities.hh"
using namespace ::Kadath::FUKA_Config;

namespace Kadath {
namespace FUKA_Solvers {
template <class config_t>
inline void bhns_xcts_setup_headon_config(config_t& bconfig) {
    using namespace ::Kadath::bco_utils;
    check_dist(bconfig(BIN_PARAMS::DIST),
               bconfig(BCO_PARAMS::MADM, NODES::BCO1),
               bconfig(BCO_PARAMS::MCH, NODES::BCO2));

    // Binary Parameters
    bconfig.set(BIN_PARAMS::REXT) = 2 * bconfig(BIN_PARAMS::DIST);

    bconfig.set(BIN_PARAMS::Q) = bconfig(BCO_PARAMS::MADM, NODES::BCO1) /
                                 bconfig(BCO_PARAMS::MCH, NODES::BCO2);
    bconfig.set(BIN_PARAMS::Q) = (bconfig.set(BIN_PARAMS::Q) > 1.0)
                                     ? 1.0 / bconfig.set(BIN_PARAMS::Q)
                                     : bconfig.set(BIN_PARAMS::Q);
    // classical Newtonian estimate
    bconfig.set(BIN_PARAMS::COM) =
        com_estimate(bconfig(BIN_PARAMS::DIST),
                     bconfig(BCO_PARAMS::MADM, NODES::BCO1),
                     bconfig(BCO_PARAMS::MCH, NODES::BCO2));
}

template <class config_t>
inline void bhns_xcts_setup_bin_config(config_t& bconfig) {
    using namespace ::Kadath::bco_utils;
    bhns_xcts_setup_headon_config(bconfig);
    auto& stages = bconfig.return_stages();

    // In the case of a head-on collision, we need do not
    // need 3.5PN estimates for the orbital velocity nor
    // the radial infall velocity.
    if (!bconfig.set_stage(STAGES::HEADON)) {
        // obtain 3PN estimate for the global, orbital omega
        KadathPNOrbitalParams(bconfig,
                              bconfig(BCO_PARAMS::MADM, NODES::BCO1),
                              bconfig(BCO_PARAMS::MCH, NODES::BCO2));
    }

    // delete ADOT, this can always be recalculated during
    // the eccentricity reduction stage
    bconfig.reset(BIN_PARAMS::ADOT);
}

template <class config_t>
void bhns_xcts_setup_space(config_t& bconfig) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::array<std::string, 2> filenames;

    filenames[0] = solve_NS_from_binary(bconfig, NODES::BCO1);
    filenames[1] = solve_BH_from_binary(bconfig, NODES::BCO2);

    // debugging only
    for (auto& f : filenames)
        if (rank == 0)
            std::cout << f << std::endl;

    if (rank == 0)
        bhns_xcts_superimposed_import(bconfig, filenames);
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.control(CONTROLS::SEQUENCES) = false;
    bconfig.open_config();
}

template <typename eos_t>
struct bhns_setup_boosted_3d {
    inline void operator()(kadath_config_boost<BCO_NS_INFO>& NSconfig,
                           kadath_config_boost<BCO_BH_INFO>& BHconfig,
                           kadath_config_boost<BIN_INFO>& bconfig) {
        using namespace ::Kadath::bco_utils;
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // open previous ns solution
        std::string nsspaceinf = NSconfig.space_filename();

        FILE* ff1 = fopen(nsspaceinf.c_str(), "r");
        Space_spheric_adapted nsspacein(ff1);
        Scalar nsconf(nsspacein, ff1);
        Scalar nslapse(nsspacein, ff1);
        Vector nsshift(nsspacein, ff1);
        Scalar nslogh(nsspacein, ff1);
        Scalar nsphi(nsspacein, ff1);
        fclose(ff1);
        // end opening bns solution

        // update NSconfig quantities before updating binary configuration file
        NSconfig.set(BCO_PARAMS::HC) =
            std::exp(get_boundary_val(0, nslogh, INNER_BC));
        NSconfig.set(BCO_PARAMS::NC) =
            EOS<eos_t, eos_var_t::DENSITY>::get(NSconfig(BCO_PARAMS::HC));
        update_config_NS_radii(nsspacein, NSconfig, 1);

        // obtain adapted NS shells for radius information and copying adapted mapping
        // later
        const Domain_shell_outer_adapted* old_outer_adapted1 =
            dynamic_cast<const Domain_shell_outer_adapted*>(
                nsspacein.get_domain(1));
        const Domain_shell_inner_adapted* old_inner_adapted1 =
            dynamic_cast<const Domain_shell_inner_adapted*>(
                nsspacein.get_domain(2));

        // setup radius field - needed for copying the adapted domain mappings.
        Scalar old_space_radius(nsspacein);
        old_space_radius.annule_hard();

        int ndominns = nsspacein.get_nbr_domains();

        for (int d = 0; d < ndominns; ++d)
            old_space_radius.set_domain(d) =
                nsspacein.get_domain(d)->get_radius();

        old_space_radius.set_domain(1) = old_outer_adapted1->get_outer_radius();
        old_space_radius.std_base();
        // end setup radius field

        // open old BH solution
        std::string bhspaceinf = BHconfig.space_filename();
        FILE* ff2 = fopen(bhspaceinf.c_str(), "r");
        Space_adapted_bh bhspacein(ff2);
        Scalar bhconf(bhspacein, ff2);
        Scalar bhlapse(bhspacein, ff2);
        Vector bhshift(bhspacein, ff2);
        fclose(ff2);
        // end open BH solution
        update_config_BH_radii(bhspacein, BHconfig, 1, bhconf);

        auto interp_field = [&](auto& space, int outer_dom, auto& old_phi) {
            const int d = outer_dom;
            const Domain_shell_inner_adapted* old_inner =
                dynamic_cast<const Domain_shell_inner_adapted*>(
                    space.get_domain(d + 1));
            // interpolate old_phi field outside of the star for import
            update_adapted_field(old_phi, d, d + 1, old_inner, INNER_BC);
        };

        // in case we used boosted TOVs, we need to import PHI
        if (NSconfig.set_field(PHI) == true)
            interp_field(nsspacein, 1, nsphi);

        // start Update config vars
        double r_max_tot = std::max(bconfig(BCO_PARAMS::RMID, NODES::BCO1),
                                    bconfig(BCO_PARAMS::RMID, NODES::BCO2));
        const double rout_sep_est =
            (bconfig(BIN_PARAMS::DIST) / 2. - r_max_tot) / 3. + r_max_tot;
        bconfig.set(BCO_PARAMS::ROUT, NODES::BCO1) = rout_sep_est;
        bconfig.set(BCO_PARAMS::ROUT, NODES::BCO2) =
            bconfig(BCO_PARAMS::ROUT, NODES::BCO1);

        const double q = bconfig.set(BIN_PARAMS::Q);
        bconfig.set(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO1) =
            1.7 * std::pow(2.0, shell_factor / q);
        bconfig.set(BCO_PARAMS::MIN_SHELL_DR, NODES::BCO2) =
            1.7 * std::pow(2.0, shell_factor / q);
        // end updating config vars

        // setup domain boundaries
        std::vector<double> out_bounds(1 + bconfig(BIN_PARAMS::OUTER_SHELLS));

        // scale outer shells by constant steps of 1/4 for the time being
        for (int e = 0; e < out_bounds.size(); ++e)
            out_bounds[e] = bconfig(BIN_PARAMS::REXT) * (1. + e * 0.25);
        // Determine optimal grid struction based on isolated solutions
        std::vector<double> NS_bounds;
        {
            auto drPsi(
                compute_drPsi(nsspacein,
                              nsconf,
                              Metric_flat(nsspacein, nsshift.get_basis()),
                              {0, 1}));
            NS_bounds = set_arb_boundsv3(bconfig, drPsi, 2, NODES::BCO1);
        }
        std::vector<double> BH_bounds;
        {
            auto drPsi(
                compute_drPsi(bhspacein,
                              bhconf,
                              Metric_flat(bhspacein, bhshift.get_basis()),
                              {0, 1}));
            BH_bounds = set_arb_boundsv3(bconfig, drPsi, 2, NODES::BCO2);
        }
        // end setup domain boundaries

        if (rank == 0) {
            std::cout << "Local bounds:" << std::endl;
            print_bounds("NS-bounds", NS_bounds);
            print_bounds("BH-bounds", BH_bounds);
            print_bounds("outer-bounds", out_bounds);
        }

        // Setup actual space
        int typer = CHEB_TYPE;
        Space_bhns space(typer,
                         bconfig(BIN_PARAMS::DIST),
                         NS_bounds,
                         BH_bounds,
                         out_bounds,
                         bconfig(BIN_PARAMS::BIN_RES),
                         bconfig(BCO_PARAMS::NINSHELLS, NODES::BCO1));
        Base_tensor basis(space, CARTESIAN_BASIS);

        const Domain_shell_inner_adapted* new_ns_inner =
            dynamic_cast<const Domain_shell_inner_adapted*>(
                space.get_domain(space.ADAPTEDNS + 1));
        const Domain_shell_outer_adapted* new_ns_outer =
            dynamic_cast<const Domain_shell_outer_adapted*>(
                space.get_domain(space.ADAPTEDNS));

        const Domain_shell_outer_homothetic* old_bh_outer =
            dynamic_cast<const Domain_shell_outer_homothetic*>(
                bhspacein.get_domain(1));

        // Update BH fields based to help with interpolation later
        update_adapted_field(bhconf, 2, 1, old_bh_outer, OUTER_BC);
        update_adapted_field(bhlapse, 2, 1, old_bh_outer, OUTER_BC);

        // Updated mapping for NS
        interp_adapted_mapping(new_ns_inner, 1, old_space_radius);
        interp_adapted_mapping(new_ns_outer, 1, old_space_radius);

        double xc1 = get_center(space, space.NS);
        double xc2 = get_center(space, space.BH);

        if (rank == 0) {
            std::cout << "xc1: " << xc1 << std::endl;
            std::cout << "xc2: " << xc2 << std::endl;
        }

        if (NSconfig.set_field(PHI) == true)
            update_adapted_field(nsphi, 1, 2, old_inner_adapted1, INNER_BC);

        // start to create the new fields
        // initialized to zero globally
        Scalar logh(space);
        logh.annule_hard();
        logh.std_base();

        Scalar conf(space);
        conf.annule_hard();

        Scalar lapse(space);
        lapse.annule_hard();

        Vector shift(space, CON, basis);
        for (int i = 1; i <= 3; i++)
            shift.set(i).annule_hard();

        Scalar phi(space);
        phi.annule_hard();
        // end create new fields

        const double ns_invw4 = set_decay(bconfig, NODES::BCO1);
        const double bh_invw4 = set_decay(bconfig, NODES::BCO2);

        if (rank == 0)
            std::cout << "WeightNS: " << bconfig(BCO_PARAMS::DECAY, NODES::BCO1)
                      << ", "
                      << "WeightBH: " << bconfig(BCO_PARAMS::DECAY, NODES::BCO2)
                      << std::endl;

        int ndom = space.get_nbr_domains();
        for (int dom = 0; dom < ndom; dom++) {
            // get an index in each domain to iterate over all colocation points
            Index new_pos(space.get_domain(dom)->get_nbr_points());

            do {
                // get cartesian coordinates of the current colocation point
                double x = space.get_domain(dom)->get_cart(1)(new_pos);
                double y = space.get_domain(dom)->get_cart(2)(new_pos);
                double z = space.get_domain(dom)->get_cart(3)(new_pos);

                // define a point shifted suitably to the stellar centers in the binary
                Point absol1(3);
                absol1.set(1) = (x - xc1);
                absol1.set(2) = y;
                absol1.set(3) = z;
                double r2 = y * y + z * z;
                double r2_1 = (x - xc1) * (x - xc1) + r2;
                double r4_1 = r2_1 * r2_1;
                double r4_invw4_1 = r4_1 * ns_invw4;
                double decay_1 = std::exp(-r4_invw4_1);

                Point absol2(3);
                absol2.set(1) = (x - xc2);
                absol2.set(2) = y;
                absol2.set(3) = z;
                double r2_2 = (x - xc2) * (x - xc2) + r2;
                double r4_2 = r2_2 * r2_2;
                double r4_invw4_2 = r4_2 * bh_invw4;
                double decay_2 = std::exp(-r4_invw4_2);

                if (dom < ndom - 1) {
                    conf.set_domain(dom).set(new_pos) =
                        1. + decay_1 * (nsconf.val_point(absol1) - 1.) +
                        decay_2 * (bhconf.val_point(absol2) - 1.);

                    lapse.set_domain(dom).set(new_pos) =
                        1. + decay_1 * (nslapse.val_point(absol1) - 1.) +
                        decay_2 * (bhlapse.val_point(absol2) - 1.);

                    logh.set_domain(dom).set(new_pos) =
                        0. + decay_1 * nslogh.val_point(absol1);

                    phi.set_domain(dom).set(new_pos) = 0;
                    if (NSconfig.set_field(PHI) == true)
                        phi.set_domain(dom).set(new_pos) +=
                            decay_1 * nsphi.val_point(absol1);

                    for (int i = 1; i <= 3; i++)
                        shift.set(i).set_domain(dom).set(new_pos) =
                            0. + decay_1 * nsshift(i).val_point(absol1) +
                            decay_2 * bhshift(i).val_point(absol2);

                } else {
                    // We have to set the compactified domain manually
                    // since the outer collocation point is always
                    // at inf which is undefined numerically
                    conf.set_domain(dom).set(new_pos) = 1.;
                    lapse.set_domain(dom).set(new_pos) = 1.;
                }
                // loop over all colocation points
            } while (new_pos.inc());
        }  // end importing fields

        // Want to make sure there is no matter around the BH
        for (int d = space.BH; d < space.OUTER; ++d) {
            logh.set_domain(d).annule_hard();
            phi.set_domain(d).annule_hard();
        }
        for (int d = space.BH; d < space.BH + 2; ++d) {
            conf.set_domain(d).annule_hard();
            lapse.set_domain(d).annule_hard();
        }

        // employ standard spectral expansion, compatible with the given paraties
        conf.std_base();
        lapse.std_base();
        logh.std_base();
        shift.std_base();
        phi.std_base();

        // save everything to a binary file
        save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
};

template <class config_t>
void bhns_xcts_superimposed_import(config_t& bconfig,
                                   std::array<std::string, 2> co_filenames) {
    using namespace Kadath::FUKA_EOS;

    // load single NS configuration
    std::string nsfilename{co_filenames[0]};
    kadath_config_boost<BCO_NS_INFO> NSconfig(nsfilename);

    // load single BH configuration
    std::string bhfilename{co_filenames[1]};
    kadath_config_boost<BCO_BH_INFO> BHconfig(bhfilename);

    // setup eos and update central density
    const std::string eos_type = NSconfig.eos<std::string>(EOS_PARAMS::EOSTYPE);
    EOS_Function_Dispatcher::dispatch<bhns_setup_boosted_3d>(bconfig,
                                                             eos_type,
                                                             NSconfig,
                                                             BHconfig,
                                                             bconfig);
};

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
