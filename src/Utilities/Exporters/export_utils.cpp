#include <cmath>
#include <exporter_utilities.hpp>
#include <point.hpp>

Kadath::Point export_utils::point_spherical(double r,
                                            double theta,
                                            double phi,
                                            double shift_x) {
    Kadath::Point abs_coords(3);

    abs_coords.set(1) = r * std::sin(theta) * std::cos(phi) + shift_x;
    abs_coords.set(2) = r * std::sin(theta) * std::sin(phi);
    abs_coords.set(3) = r * std::cos(theta);

    return abs_coords;
}

double export_utils::lagrange_gen_k(int interp_order,
                                    double x,
                                    const double* xp,
                                    const double* yp) {

    double temp = 0;
    for (int i = 0; i < interp_order; ++i) {
        double temp2 = 1.;
        for (int j = 0; j < interp_order; ++j) {
            if (i == j)
                continue;
            temp2 *= (x - xp[j]) / (xp[i] - xp[j]);
        }
        temp += yp[i] * temp2;
    }

    return temp;
}

std::string export_utils::throw_no_multithreaded_support_error(
    const std::string not_implemented) {
    std::string error_msg = "\n" + not_implemented +
                            " is not available when the compile time "
                            "definition DEFAULT_KAD_MEM is not defined.\n";
    error_msg +=
        "DEFAULT_KAD_MEM is documented in $HOME_KADATH/include/memory.hpp and "
        "should only be defined when multi-threaded\n";
    error_msg +=
        "support is needed for exporters.  It is highly recommended to have a "
        "separate build of FUKA with DEFAULT_KAD_MEM\n";
    error_msg +=
        "for export only as enabling DEFAULT_KAD_MEM and attempting to using "
        "the initial data codes will be drastically slower.\n";
    return error_msg;
}