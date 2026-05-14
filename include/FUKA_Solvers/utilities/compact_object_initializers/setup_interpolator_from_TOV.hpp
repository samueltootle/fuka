#pragma once
#include<memory>
#include<cmath>
template <typename tov_t>
auto setup_interpolator_from_TOV(tov_t& tov) {
    const size_t max_iter = tov.state.size();

    std::unique_ptr<double[]> radius_lin_ptr{new double[max_iter]};
    std::unique_ptr<double[]> conf_lin_ptr{new double[max_iter]};
    std::unique_ptr<double[]> lapse_lin_ptr{new double[max_iter]};
    std::unique_ptr<double[]> rho_lin_ptr{new double[max_iter]};

    for (size_t j = 0; j < max_iter; ++j) {
        // lapse est.
        auto phi = tov.state[j][tov.PHI];

        lapse_lin_ptr[j] = std::exp(phi);

        conf_lin_ptr[j] = tov.state[j][tov.CONF];

        radius_lin_ptr[j] = tov.state[j][tov.RISO];

        rho_lin_ptr[j] = tov.state[j][tov.RHOB];
    }

    // linear interpolation of conf(r_isotropic), lapse(r_isotropic), and
    // rho(r_isotropic) The order here dictates the enum ltpQ enum
    linear_interp_t<double, 3> ltp(max_iter,
                                   std::move(radius_lin_ptr),
                                   std::move(lapse_lin_ptr),
                                   std::move(rho_lin_ptr),
                                   std::move(conf_lin_ptr));
    return ltp;
}