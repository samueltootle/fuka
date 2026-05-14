#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath_adapted_bh.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/
using namespace Kadath::FUKA_Config;

namespace Kadath::FUKA_Solvers {

template <typename config_t>
static void write_bh_init_setup_tofile_XCTS(Space_adapted_bh& space,
                                     config_t& bconfig) {
    Base_tensor basis(space, CARTESIAN_BASIS);

    // setup fields
    Scalar lapse(space);
    lapse = 1.;
    lapse.set_domain(0).annule_hard();
    lapse.set_domain(1).annule_hard();

    Scalar conf(lapse);

    // set a better estimate for PSI on the horizon
    conf.set_domain(2) = Kadath::bco_utils::psi;
    conf.std_base();
    lapse.std_base();

    Vector shift(space, CON, basis);
    for (int i = 1; i <= 3; i++)
        shift.set(i).annule_hard();
    shift.std_base();
    // end setup fields

    Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift);
}

/**
 * @brief Set the initial guess for a BH
 *
 * @tparam config_t Config object type
 * @param bconfig Config object
 */
template <typename config_t>
void setup_3d_BH_xcts(config_t& bconfig) {
    auto& fields = bconfig.return_fields();

    int type_coloc = CHEB_TYPE;
    Dim_array res(bconfig(BCO_PARAMS::DIM));
    res.set(0) = bconfig(BCO_PARAMS::BCO_RES);
    res.set(1) = bconfig(BCO_PARAMS::BCO_RES);
    res.set(2) = bconfig(BCO_PARAMS::BCO_RES) - 1;

    Point center(bconfig(BCO_PARAMS::DIM));
    for (int i = 1; i <= bconfig(BCO_PARAMS::DIM); i++)
        center.set(i) = 0;

    const int shells = (int)bconfig(BCO_PARAMS::NSHELLS);
    int ndom = 4 + shells;
    std::vector<double> bounds(ndom - 1);

    // Estimate Radius of BH based on Schwarzschild radius and an
    // estimate for Psi on the horizon based on prev. BH solutions
    if (!bconfig.control(CONTROLS::USE_CONFIG_VARS)) {
        bconfig.set(BCO_PARAMS::RIN) =
            bconfig(BCO_PARAMS::MCH) * Kadath::bco_utils::invpsisq;
        bconfig.set(BCO_PARAMS::RMID) =
            2 * bconfig(BCO_PARAMS::MCH) * Kadath::bco_utils::invpsisq;
        bconfig.set(BCO_PARAMS::ROUT) = 4 * bconfig(BCO_PARAMS::RMID);
    }
    Kadath::bco_utils::set_isolated_BH_bounds(bounds, bconfig);
#ifdef DEBUG
    Kadath::bco_utils::print_bounds("BH", bounds);
#endif
    Space_adapted_bh space(type_coloc, center, res, bounds);

    write_bh_init_setup_tofile_XCTS(space, bconfig);
}

/** @}*/
}  // namespace Kadath::FUKA_Solvers