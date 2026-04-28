#pragma once

namespace Kadath::FUKA_Solvers {
inline void check_dist(double dist, double M1, double M2, double garbage_factor) {
    auto M = M1 + M2;
    auto const garbage_dist = garbage_factor * M;
    auto const recommended_dist = 8. * M;
    if (dist <= garbage_dist) {
        std::cerr << "Distance is set to (" << dist << ") which may not give results. \nSet to ("
                  << recommended_dist << ") for something reasonable.\n";
        std::_Exit(EXIT_FAILURE);
    }
}
}    // namespace Kadath::FUKA_Solvers