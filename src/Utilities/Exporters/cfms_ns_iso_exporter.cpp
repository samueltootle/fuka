#include "Solvers/ns_isotropic/ns_isotropic_exporter.hpp"

namespace Kadath::FUKA_Solvers {
#ifdef DEFAULT_KAD_MEM
CFMS_NS_ISO_Exporter::CFMS_NS_ISO_Exporter(CFMS_NS_ISO_Exporter const& r) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;
    h_cut = r.h_cut;
    eos_file = r.eos_file;
    eos_type = r.eos_type;

    space.reset(new space_t((*r.get_space())));
    // Copy ID fields
    lap_Aterm.reset(new Scalar(*space, *r.lap_Aterm.get()));
    Nu.reset(new Scalar(*space, *r.Nu.get()));
    logh.reset(new Scalar(*space, *r.logh.get()));
    lap_Bterm.reset(new Scalar(*space, *r.lap_Bterm.get()));

    // For uniform rotation solutions
    if (r.lap_omega_term)
        lap_omega_term.reset(new Scalar(*space, *r.lap_omega_term.get()));
    else
        lap_omega_term = nullptr;
    // For differential rotation solutions
    if (r.omega)
        omega.reset(new Scalar(*space, *r.omega.get()));
    else
        omega = nullptr;

    // Copy Computed Terms
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    metric_A.reset(new Scalar(*space, *r.metric_A.get()));
    metric_B.reset(new Scalar(*space, *r.metric_B.get()));
    metric_omega.reset(new Scalar(*space, *r.metric_omega.get()));
    domega_dr.reset(new Scalar(*space, *r.domega_dr.get()));
    domega_dt.reset(new Scalar(*space, *r.domega_dt.get()));
    fluidvel.reset(new Scalar(*space, *r.fluidvel.get()));

    bconfig.reset(new config_t(*r.bconfig));

    export_ready = false;

    populate_quants();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor,
    // *lapse, *shift); std::cout << "copy\n";
}

CFMS_NS_ISO_Exporter& CFMS_NS_ISO_Exporter::operator=(
    const CFMS_NS_ISO_Exporter& b) {
    if (this == &b)
        return *this;

    CFMS_NS_ISO_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
}
#else
CFMS_NS_ISO_Exporter::CFMS_NS_ISO_Exporter(CFMS_NS_ISO_Exporter const& r) {
    std::string error_msg = export_utils::throw_no_multithreaded_support_error(
        "CFMS_NS_ISO_Exporter - Copy Constructor");
    throw std::runtime_error(error_msg);
}

CFMS_NS_ISO_Exporter& CFMS_NS_ISO_Exporter::operator=(
    const CFMS_NS_ISO_Exporter& b) {
    std::string error_msg = export_utils::throw_no_multithreaded_support_error(
        "CFMS_NS_ISO_Exporter - Assignment operator");
    throw std::runtime_error(error_msg);
}
#endif

void CFMS_NS_ISO_Exporter::initialize_eos() {
    using namespace Kadath::FUKA_EOS;
    EOS_initialize::init(*bconfig);
    eos_file = (*bconfig).template eos<std::string>(
        Kadath::FUKA_Config::EOS_PARAMS::EOSFILE);
    eos_type = (*bconfig).template eos<std::string>(
        Kadath::FUKA_Config::EOS_PARAMS::EOSTYPE);
}

void CFMS_NS_ISO_Exporter::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen(spacein.c_str(), "r");

    space.reset(new space_t{ff1});
    lap_Aterm.reset(new Scalar(*space.get(), ff1));
    Nu.reset(new Scalar(*space.get(), ff1));
    logh.reset(new Scalar(*space.get(), ff1));
    lap_Bterm.reset(new Scalar(*space.get(), ff1));

    lap_Aterm->coef();
    Nu->coef();
    logh->coef();
    lap_Bterm->coef();

    if (bconfig->set_field(Kadath::FUKA_Config::BCO_FIELDS::LAP_WTERM)) {
        lap_omega_term.reset(new Scalar(*space.get(), ff1));
        lap_omega_term->coef();

        if (bconfig->field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)) {
            omega.reset(new Scalar(*space.get(), ff1));
            omega->coef();
#ifdef DEBUG
            std::cout << "**** Reading Differentially rotating solution ****\n";
#endif
        } else {
#ifdef DEBUG
            std::cout << "**** Reading Uniformly rotating solution ****\n";
#endif
        }
    } else {
#ifdef DEBUG
        std::cout << "**** Reading non-rotating, spherical solution ****\n";
#endif
    }

    fclose(ff1);
    ndom = space->get_nbr_domains();
}

void CFMS_NS_ISO_Exporter::extract_computed_grid_functions() {
    using namespace Kadath::FUKA_EOS;

    System_of_eqs syst(*space);

    Param p;
    EOS_Function_Dispatcher::dispatch<set_eos_ope_struct>(*bconfig,
                                                          eos_type,
                                                          syst,
                                                          p);

    // Fields - must be initialized before common setup
    syst.add_cst("nu", *Nu);
    syst.add_cst("lapAterm", *lap_Aterm);
    syst.add_cst("lapBterm", *lap_Bterm);

    if (bconfig->field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)) {
#ifdef DEBUG
        std::cout << "**** Importing differential rotation profile ****\n";
#endif
        syst.add_cst("Omega", *omega);
    } else {
#ifdef DEBUG
        std::cout << "**** Importing uniform rotation profile ****\n";
#endif
        syst.add_cst("Omega",
                     (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA));
    }
    syst.add_def("N = exp(nu)");
    syst.add_def("A = exp(lapAterm - nu)");
    syst.add_def("B = (divrsint(lapBterm) + 1) / N");

    if (bconfig->set_field(Kadath::FUKA_Config::BCO_FIELDS::LAP_WTERM)) {
        syst.add_cst("wrsint", *lap_omega_term);
        syst.add_def("w = divrsint(wrsint)");
        syst.add_def("drw = dr(w)");
        syst.add_def("dtw = dt(w)");

        syst.add_def("UphiU = 1 / N * (Omega - w)");

        fluidvel.reset(new Scalar(syst.give_val_def("UphiU")));
        fluidvel->coef();

        metric_omega.reset(new Scalar(syst.give_val_def("w")));
        metric_omega->coef();
        domega_dr.reset(new Scalar(syst.give_val_def("drw")));
        domega_dr->coef();
        domega_dt.reset(new Scalar(syst.give_val_def("dtw")));
        domega_dt->coef();
    } else {
#ifdef DEBUG
        std::cout << "**** Importing non-rotation profile ****\n";
#endif
        fluidvel.reset(new Scalar(*space));
        fluidvel->annule_hard();
        fluidvel->std_base();

        metric_omega.reset(new Scalar(*space));
        metric_omega->annule_hard();
        metric_omega->std_base();
        domega_dr.reset(new Scalar(*space));
        domega_dr->annule_hard();
        domega_dr->std_base();
        domega_dt.reset(new Scalar(*space));
        domega_dt->annule_hard();
        domega_dt->std_base();
    }
    metric_A.reset(new Scalar(syst.give_val_def("A")));
    metric_A->coef();
    metric_B.reset(new Scalar(syst.give_val_def("B")));
    metric_B->coef();
    lapse.reset(new Scalar(syst.give_val_def("N")));
    lapse->coef();
}

void CFMS_NS_ISO_Exporter::populate_quants() {
    if (quants.capacity() != ISO_VARS::NUM_ISO_VARS) {
        for (size_t i = 0; i < ISO_VARS::NUM_ISO_VARS; ++i)
            quants.push_back(std::cref(*Nu));
    }
    quants[ISO_VARS::ISO_METRIC_A] = std::cref(*metric_A);
    quants[ISO_VARS::ISO_METRIC_B] = std::cref(*metric_B);
    quants[ISO_VARS::ISO_ALPHA] = std::cref(*lapse);
    quants[ISO_VARS::ISO_METRIC_OMEGA] = std::cref(*metric_omega);
    quants[ISO_VARS::ISO_DOMEGA_DR] = std::cref(*domega_dr);
    quants[ISO_VARS::ISO_DOMEGA_DTHETA] = std::cref(*domega_dt);

    // Fluid related quantities
    quants[ISO_VARS::ISO_H] = std::cref(*logh);
    quants[ISO_VARS::ISO_U] = std::cref(*fluidvel);
    export_ready = true;
}

CFMS_NS_ISO_Exporter::interp_ary_t CFMS_NS_ISO_Exporter::interpolate_pointwise(
    double const& x,
    double const& y,
    double const& z) {
    double r2_xy = x * x + y * y;
    double r_xy = std::sqrt(r2_xy);
    Point abs_coords(ndim);
    abs_coords.set(1) = r_xy;
    abs_coords.set(2) = z;

    for (size_t k = 0; k < ISO_VARS::NUM_ISO_VARS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }

    return quant_vals;
}

CFMS_NS_ISO_Exporter::interp_ary_t
CFMS_NS_ISO_Exporter::interpolate_pointwise_subset(
    double const& x,
    double const& y,
    double const& z,
    std::vector<CFMS_NS_ISO_Exporter::ISO_VARS> slice) {
    double r2_xy = x * x + y * y;
    double r_xy = std::sqrt(r2_xy);
    Point abs_coords(ndim);
    abs_coords.set(1) = r_xy;
    abs_coords.set(2) = z;

    for (const auto k : slice) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }

    return quant_vals;
}

CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise(
    double const& x,
    double const& y,
    double const& z) {
    using namespace Kadath::FUKA_EOS;
    return EOS_Function_Dispatcher::dispatch<
        CFMS_NS_ISO_Exporter::export_pointwise_imp>(*bconfig,
                                                    eos_type,
                                                    *this,
                                                    x,
                                                    y,
                                                    z);
}

CFMS_NS_ISO_Exporter::output_ary_t
CFMS_NS_ISO_Exporter::export_pointwise__spherical(double const& x,
                                                  double const& y,
                                                  double const& z) {
    using namespace Kadath::FUKA_EOS;
    return EOS_Function_Dispatcher::dispatch<
        CFMS_NS_ISO_Exporter::export_pointwise__spherical_imp>(*bconfig,
                                                               eos_type,
                                                               *this,
                                                               x,
                                                               y,
                                                               z);
}
}  // namespace Kadath::FUKA_Solvers