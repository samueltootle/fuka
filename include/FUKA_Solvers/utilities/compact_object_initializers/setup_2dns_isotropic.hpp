
#pragma once
#include "Configurator/config_bco.hpp"
#include "kadath_spheric.hpp"
#include "setup_interpolator_from_TOV.hpp"
#include "setup_ns_config_from_TOV.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/
namespace Kadath::FUKA_Solvers {
using namespace Kadath::FUKA_Config;

/**
 * @brief Write the initial guess to file
 * for a 2D Neutron star in
 * isotropic coordinates based on
 * arxiv.org:1003.5015
 *
 * @tparam eos_t EOS wrapper type
 * @tparam config_t configuration file type
 * @param space numerical space
 * @param bconfig the configuration file
 * @param tov the TOV solution to use for interpolation
 */
template <typename tov_t, typename config_t>
void write_2dns_isotropic_init_setup_tofile(Space_polar_adapted& space,
                                            config_t& bconfig,
                                            tov_t& tov) {
    using eos_t = typename tov_t::eos_t;

    enum ltpQ { LAPSE = 0, RHO, CONF };

    const int ndom = space.get_nbr_domains();

    // setup fields
    Scalar lapse(space);
    lapse = 1.;

    Scalar conf(lapse);

    // H = log(h), the logarithm of the specific enthalpy
    Scalar logh(space);
    logh.annule_hard();

    // interpolate TOV solution
    auto lintp = setup_interpolator_from_TOV(tov);

    Scalar r_field(space);
    r_field.annule_hard();
    for (auto i = 0; i < ndom; ++i)
        r_field.set_domain(i) = space.get_domain(i)->get_radius();
    auto npts = space.get_domain(ndom - 1)->get_nbr_points();
    Index pos_c(npts);
    pos_c.set(0) = npts(0) - 1;  // outer_bc

    for (auto i = 0; i < npts(1); ++i) {
        pos_c.set(1) = i;
        r_field.set_domain(ndom - 1).set(pos_c) = 1e10;
    }
    pos_c.set_start();

    // update the fields based on TOV solution for a given domain
    auto update_fields = [&](const size_t dom) {
        Index pos(space.get_domain(dom)->get_nbr_points());
        do {
            double rval = r_field(dom)(pos);
            auto all_ltp = lintp.interpolate_all(rval);
            auto rho = (all_ltp[ltpQ::RHO] <= 0) ? 1e-15 : all_ltp[ltpQ::RHO];
            auto h = EOS<eos_t, eos_var_t::DENSITY>::h_cold__rho(rho);

            if (dom == 0 && pos(0) == 0 && pos(1) == 0)
                bconfig.set(BCO_PARAMS::HC) = h;

            logh.set_domain(dom).set(pos) = (h <= 1) ? 0. : std::log(h);
            lapse.set_domain(dom).set(pos) = all_ltp[ltpQ::LAPSE];
            conf.set_domain(dom).set(pos) = all_ltp[ltpQ::CONF];
        } while (pos.inc());
    };
    for (int i = 0; i < ndom; ++i)
        update_fields(i);

    for (int d = space.ADAPTED_INNER; d < ndom; ++d)
        logh.set_domain(d).annule_hard();

    // Fix compactified domain and additional shells
    for (auto d = 3; d < ndom; ++d) {
        pos_c.set_start();
        auto decay_factor = r_field(d)(pos_c) / r_field(d);
        if (d > 3) {
            pos_c.set(0) = npts(0) - 1;
            conf.set_domain(d) = 1 + (conf(d - 1)(pos_c) - 1) * decay_factor;
            lapse.set_domain(d) = 1 + (lapse(d - 1)(pos_c) - 1) * decay_factor;
        } else {
            conf.set_domain(d) = 1 + (conf(d)(pos_c) - 1) * decay_factor;
            lapse.set_domain(d) = 1 + (lapse(d)(pos_c) - 1) * decay_factor;
        }
    }

    Scalar A(conf * conf);
    A.std_base();
    lapse.std_base();

    // We store the quantities that appear in the lapace terms
    // for solver stability.  These will need to be transformed
    // to obtain the full metric.
    Scalar nu(log(lapse));
    Scalar lap_aterm(nu + log(A));

    logh.std_base();
    nu.std_base();
    lap_aterm.std_base();

    Scalar tmp(lapse * A - 1);
    Scalar lap_Bterm = Scalar(tmp.mult_r().mult_sin_theta());
    // end setup fields

    bco_utils::save_to_file(space, bconfig, lap_aterm, nu, logh, lap_Bterm);
}

/**
 * @brief Functor for setting up 2D Neutron star in isotropic coordinates
 * @tparam eos_t EOS wrapper type
 * @param bconfig Configuration object
 * @param mass_fixing_idx Index of the mass
*/
template <class eos_t>
struct setup_2dns_isotropic_functor {

    template <typename config_t>
    void operator()(config_t& bconfig, size_t mass_fixing_idx) {
        using namespace Kadath::FUKA_EOS;

        int type_coloc = CHEB_TYPE;
        auto const& dim = bconfig(BCO_PARAMS::DIM);
        Dim_array res(dim);
        res.set(0) = bconfig(BCO_PARAMS::BCO_RES);
        res.set(1) = bconfig(BCO_PARAMS::BCO_RES);

        Point center(dim);
        for (int i = 1; i <= dim; i++)
            center.set(i) = 0;

        const int shells = (int)bconfig(BCO_PARAMS::NSHELLS);
        int ndom = 4 + shells;

        std::unique_ptr<Kadath::Margherita::MargheritaTOV<eos_t>> tov =
            setup_ns_config_from_TOV<eos_t>(bconfig, mass_fixing_idx);

        Array<double> bounds(ndom - 1);

        bounds.set(0) = bconfig(RIN);
        bounds.set(1) = bconfig(RMID);
        bounds.set(2) = bconfig(ROUT);

        for (int shell = 1, b = 3; shell <= shells; ++shell, ++b) {
            bounds.set(b) = bounds(b - 1) * 2;
        }

        // generate a full single star space including compactification to infinity
        Space_polar_adapted space(type_coloc, center, res, bounds);

        write_2dns_isotropic_init_setup_tofile(space, bconfig, *tov);
    }
};

/** @}*/
}  // namespace Kadath::FUKA_Solvers