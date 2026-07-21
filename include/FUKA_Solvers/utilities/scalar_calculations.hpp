#pragma once
#include "kadath.hpp"

namespace Kadath::FUKA_Solvers {
template <class space_t, class metric_t>
Scalar compute_drPsi(space_t& space,
                     Scalar& conf,
                     metric_t metric,
                     std::vector<int> excluded_doms,
                     int bound_dom = -1);

template <class space_t, class metric_t>
Scalar compute_ddrPsi(space_t& space,
                      Scalar& conf,
                      metric_t metric,
                      std::vector<int> excluded_doms,
                      int bound_dom = -1);
}  // namespace Kadath::FUKA_Solvers

#include "scalar_calculations.cpp"