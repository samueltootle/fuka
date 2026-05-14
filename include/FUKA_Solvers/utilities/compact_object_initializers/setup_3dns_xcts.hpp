#pragma once
#include "EOS/EOS.hh"
#include "EOS/FUKA_EOS_Utilities.hh"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "setup_interpolator_from_TOV.hpp"
#include "setup_ns_config_from_TOV.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/
using namespace Kadath::FUKA_Config;

namespace Kadath::FUKA_Solvers {

template <typename tov_t, typename config_t>
static void write_ns_init_setup_tofile__xcts(Space_spheric_adapted& space,
                                     config_t& bconfig,
                                     tov_t& tov) {
    using eos_t = typename tov_t::eos_t;

    enum ltpQ { LAPSE = 0, RHO, CONF };

    const int ndom = space.get_nbr_domains();
    Base_tensor basis(space, CARTESIAN_BASIS);
    // setup fields
    Scalar lapse(space);
    lapse = 1.;

    Scalar conf(lapse);

    // H = log(h), the logarithm of the specific enthalpy
    Scalar logh(space);
    logh.annule_hard();

    // interpolate TOV solution
    auto lintp = setup_interpolator_from_TOV(tov);

    // generate a radius field for use with the 1D interpolators
    CoordFields<Space_spheric_adapted> cfg(space);
    Scalar r_field(cfg.radius());

    // update the fields based on TOV solution for a given domain
    auto update_fields = [&](const size_t dom) {
        Index pos(space.get_domain(dom)->get_nbr_points());
        do {
            double rval = r_field(dom)(pos);
            auto all_ltp = lintp.interpolate_all(rval);
            auto rho = (all_ltp[ltpQ::RHO] <= 0) ? 1e-15 : all_ltp[ltpQ::RHO];

            auto h = EOS<eos_t, eos_var_t::DENSITY>::h_cold__rho(rho);
            logh.set_domain(dom).set(pos) = (h <= 1) ? 0. : std::log(h);
            lapse.set_domain(dom).set(pos) = all_ltp[ltpQ::LAPSE];
            conf.set_domain(dom).set(pos) = all_ltp[ltpQ::CONF];
        } while (pos.inc());
    };
    for (int i = 0; i < 3; ++i)
        update_fields(i);

    for (int d = 2; d < ndom; ++d)
        logh.set_domain(d).annule_hard();

    Vector shift(space, CON, basis);
    for (int i = 1; i <= 3; i++)
        shift.set(i).annule_hard();
    shift.std_base();
    logh.std_base();
    conf.std_base();
    lapse.std_base();
    // end setup fields

    Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
}

/**
 * @brief Set the initial guess for 3D Neutron star in
 * XCTS coordinates
 *
 * @tparam eos_t EOS wrapper type
 */
template <class eos_t>
struct setup_3dns_xcts_functor {

    template <typename config_t>
    void operator()(config_t& bconfig, size_t mass_fixing_idx) {
        auto& fields = bconfig.return_fields();

        int type_coloc = CHEB_TYPE;
        auto const dim = bconfig(BCO_PARAMS::DIM);
        Dim_array res(dim);
        res.set(0) = bconfig(BCO_PARAMS::BCO_RES);
        res.set(1) = bconfig(BCO_PARAMS::BCO_RES);
        res.set(2) = bconfig(BCO_PARAMS::BCO_RES) - 1;

        Point center(dim);
        for (int i = 1; i <= dim; i++)
            center.set(i) = 0;

        const int shells = (int)bconfig(BCO_PARAMS::NSHELLS);
        int ndom = 4 + shells;

        std::unique_ptr<Kadath::Margherita::MargheritaTOV<eos_t>> tov =
            setup_ns_config_from_TOV<eos_t>(bconfig, mass_fixing_idx);
        std::vector<double> bounds(ndom - 1);

        Kadath::bco_utils::set_NS_bounds(bounds, bconfig);

        // generate a full single star space including compactification to infinity
        Space_spheric_adapted space(type_coloc, center, res, bounds);

        write_ns_init_setup_tofile__xcts(space, bconfig, *tov);
    }
};
/** @}*/
}  // namespace Kadath::FUKA_Solvers
