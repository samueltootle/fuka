#include "Configurator/config_bco.hpp"
#include "EOS/EOS.hh"
#include "EOS/FUKA_EOS_Utilities.hh"
#include "Solvers/ns_isotropic/ns_isotropic_exporter.hpp"
#include "bco_utilities.hpp"
// #include "EOS/standalone/tov.hh"

// Kadath includes
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_polar.hpp"

#include <mpi.h>

using namespace Kadath;
using namespace Kadath::FUKA_Config;
using namespace Kadath::FUKA_EOS;

using config_t = kadath_config_boost<BCO_NS_INFO>;
using input_reader_t = Kadath::FUKA_Solvers::CFMS_NS_ISO_Exporter;

template <class eos_t>
struct ISO_XCTS_convert {
    void operator()(config_t& in_bconfig) {
        // Convert from QI solution
        input_reader_t input_reader(in_bconfig.config_filename_abs());

        auto& old_space = input_reader.get_space();
        Scalar old_space_radius(*old_space);
        old_space_radius = 0.;
        // get the radius from each domain
        for (int i = 0; i < old_space->get_nbr_domains(); ++i)
            old_space_radius.set_domain(i) =
                old_space->get_domain(i)->get_radius();
        // get the adapted radius of the adapted domain
        const Domain_polar_shell_outer_adapted* old_outer_adapted =
            dynamic_cast<const Domain_polar_shell_outer_adapted*>(
                old_space->get_domain(old_space->ADAPTED_OUTER));
        old_space_radius.set_domain(old_space->ADAPTED_OUTER) =
            old_outer_adapted->get_outer_radius();

        // define a standard decomposition, compatible with the parity of this field
        old_space_radius.std_base();
        //end setup old radius field

        // get the minimal and maximal radius from the adapted domain
        auto [r_min, r_max] = Kadath::bco_utils::get_rmin_rmax(*old_space, 1);

        config_t bconfig = in_bconfig;
        bconfig.set(BCO_PARAMS::DIM) = 3;
        bconfig.set_filename("./initns_xcts.info");

        // number of dimensions
        const int dim = bconfig(DIM);

        // collocation point type and number of collocation points per domain
        int type_coloc = CHEB_TYPE;
        Dim_array res(dim);
        res.set(0) = bconfig(BCO_RES);
        res.set(1) = bconfig(BCO_RES);
        res.set(2) = bconfig(BCO_RES) - 1;

        // physical center of the system, arbitrary
        Point center(dim);
        for (int i = 1; i <= dim; i++)
            center.set(i) = 0;

        // domain boundaries, i.e. radii for each domain of the single star space
        int ndom = 4;
        Array<double> bounds(ndom - 1);
        bounds.set(0) = bconfig(RIN);
        bounds.set(1) = bconfig(RMID);
        bounds.set(2) = bconfig(ROUT);

        // generate a full single star space including compactification to infinity
        Space_spheric_adapted space(type_coloc, center, res, bounds);

        // use a basis of cartesian type
        Base_tensor basis(space, CARTESIAN_BASIS);

        // get adapted domains to update the radius
        const Domain_shell_outer_adapted* new_outer_adapted =
            dynamic_cast<const Domain_shell_outer_adapted*>(
                space.get_domain(1));
        const Domain_shell_inner_adapted* new_inner_adapted =
            dynamic_cast<const Domain_shell_inner_adapted*>(
                space.get_domain(2));

        // update adapted domain mapping
        Kadath::bco_utils::interp_adapted_mapping(new_outer_adapted,
                                                  1,
                                                  old_space_radius);
        Kadath::bco_utils::interp_adapted_mapping(new_inner_adapted,
                                                  1,
                                                  old_space_radius);

        Scalar logh(space);
        logh.annule_hard();

        Scalar omega(logh);

        // conformal factor is initialized to one, representing flat space as an initial guess
        Scalar conf(space);
        conf = 1.;

        // same for the lapse
        Scalar lapse(conf);

        Vector shift(space, CON, basis);
        shift.annule_hard();

        auto convert_fields = [&](const size_t dom) {
            auto npts = space.get_domain(dom)->get_nbr_points();
            Index pos(npts);
            Val_domain xx = space.get_domain(dom)->get_cart(1);
            Val_domain yy = space.get_domain(dom)->get_cart(2);
            Val_domain zz = space.get_domain(dom)->get_cart(3);
            do {
                if (!(dom == ndom - 1 && pos(0) == npts(0) - 1)) {
                    auto x = xx(pos);
                    auto y = yy(pos);
                    auto z = zz(pos);
                    auto all_data = input_reader.export_pointwise(x, y, z);
                    auto const rho = all_data[input_reader_t::OUTPUT_VARS::RHO];
                    auto const h = EOS<eos_t, DENSITY>::h_cold__rho(rho);
                    logh.set_domain(dom).set(pos) =
                        (std::log(h) < 0) ? 0. : std::log(h);
                    lapse.set_domain(dom).set(pos) =
                        all_data[input_reader_t::OUTPUT_VARS::ALPHA];

                    // approximation
                    auto const psi4 =
                        all_data[input_reader_t::OUTPUT_VARS::G11];
                    auto const psi2 = std::sqrt(psi4);
                    auto const psi = std::sqrt(psi2);
                    conf.set_domain(dom).set(pos) = psi;

                    auto iso_omega(input_reader.get_omega());
                    double r2_xy = x * x + y * y;
                    double r_xy = std::sqrt(r2_xy);
                    Point abs_coords(2);
                    abs_coords.set(1) = r_xy;
                    abs_coords.set(2) = z;
                    omega.set_domain(dom).set(pos) =
                        iso_omega->val_point(abs_coords);

                    shift.set(1).set_domain(dom).set(pos) =
                        all_data[input_reader_t::OUTPUT_VARS::BETA1];
                    shift.set(2).set_domain(dom).set(pos) =
                        all_data[input_reader_t::OUTPUT_VARS::BETA2];
                    shift.set(3).set_domain(dom).set(pos) =
                        all_data[input_reader_t::OUTPUT_VARS::BETA3];
                } else {
                    lapse.set_domain(dom).set(pos) = 1.;
                    conf.set_domain(dom).set(pos) = 1.;
                    shift.set(1).set_domain(dom).set(pos) = 0;
                    shift.set(2).set_domain(dom).set(pos) = 0;
                    shift.set(3).set_domain(dom).set(pos) = 0;
                    logh.set_domain(dom).set(pos) = 0.;
                }
            } while (pos.inc());
        };
        for (int i = 0; i < ndom; ++i)
            convert_fields(i);

        for (int d = 2; d < ndom; ++d)
            logh.set_domain(d).annule_hard();

        // shift is zero in the beginning and in general for TOV solutions

        shift.std_base();
        logh.std_base();
        conf.std_base();
        lapse.std_base();
        omega.std_base();

        for (int i = 1; i <= 3; ++i)
            shift.set(i).coef();
        logh.coef();
        conf.coef();
        lapse.coef();
        omega.coef();
        // finished create the fields
        for (int i = 0; i < BCO_FIELDS::NUM_BCO_FIELDS; ++i)
            bconfig.set_field(i) = false;
        bconfig.set_field(BCO_FIELDS::LAPSE) = true;
        bconfig.set_field(BCO_FIELDS::CONF) = true;
        bconfig.set_field(BCO_FIELDS::SHIFT) = true;
        bconfig.set_field(BCO_FIELDS::LOGH) = true;
        bconfig.set_field(BCO_FIELDS::DIFF_OMEGA) = true;
        // write space and config to files on one processor
        bconfig.set_stage(PRE) = false;
        bco_utils::save_to_file(space,
                                bconfig,
                                conf,
                                lapse,
                                shift,
                                logh,
                                omega);
    }
};

int main(int argc, char** argv) {
    // initialize MPI
    int rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) {
        cerr << "Error starting MPI" << endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // expecting a configuration file on execution
    if (argc < 2) {
        std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
        std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
        std::_Exit(EXIT_FAILURE);
    }

    // load the INPUT solution in QI coordinates
    std::string ifilename{argv[1]};
    config_t bconfig(ifilename);

    // setup the EOS
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);
    EOS_Function_Dispatcher::dispatch<ISO_XCTS_convert>(bconfig,
                                                        eos_type,
                                                        bconfig);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
