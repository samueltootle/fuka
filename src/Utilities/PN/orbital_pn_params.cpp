#include "bco_utilities.hpp"

namespace Kadath {
namespace bco_utils {

void KadathPNOrbitalParams(kadath_config_boost<BIN_INFO>& bconfig,
                           double Madm1,
                           double Madm2) {

    // FIXME rescale if non-geometric units
    //double units = bconfig(QPIG);

    double r = bconfig(DIST);
    double axis = bconfig(COM);

    double r1 = std::abs(-r / 2.0 + axis);
    double r2 = std::abs(r / 2.0 + axis);

    double M = Madm1 + Madm2;
    double mu = Madm1 * Madm2 / M;
    double nu = mu / M;
    double gam = M / r;
    double X1 = Madm1 / M;
    double X2 = Madm2 / M;
    double nab = X1 - X2;

    double r1p = r1 + 1;
    double r2p = r2 + 1;
    double lnr0 = X1 * std::log(r1p) + X2 * std::log(r2p);

    double Omegasq =
        M / std::pow(r, 3) *
        (1.0 + (-3.0 + nu) * gam +
         (6.0 + 41.0 / 4.0 * nu + std::pow(nu, 2)) * std::pow(gam, 2) +
         (-10.0 +
          (-75707.0 / 840.0 + 41.0 / 64.0 * std::pow(M_PI, 2) +
           22.0 * (std::log(r) - lnr0)) *
              nu +
          19.0 / 2.0 * std::pow(nu, 2) + std::pow(nu, 3)) *
             std::pow(gam, 3));

    double omega_pn = std::sqrt(Omegasq);
    double rdot = -64.0 / 5.0 * std::pow(M, 3) * nu / std::pow(r, 3) *
                  (1.0 + gam * (-1751.0 / 336.0 - 7.0 / 4.0 * nu));

    double adot_pn = rdot / r;

    bconfig.set(ADOT) = adot_pn;
    bconfig.set(GOMEGA) = omega_pn;
}
}  // namespace bco_utils
}  // namespace Kadath