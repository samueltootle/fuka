/*
 * Copyright 2023
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author:
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once
#include "FUKA_Solvers/utilities/compact_object_initializers/setup_2dns_isotropic.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/sequences/sequence_utilities.hpp"
#include "ns_isotropic_diff_rot_regrid.hpp"
#include "ns_isotropic_diff_rot_solver.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief Driver to compute a stationary solution for a given resolution
 *
 * @tparam config_t Configurator type
 * @param bconfig BH config file
 * @param outputdir output location
 * @return int error code
 */
template <typename config_t>
int ns_isotropic_diff_rot_stationary_driver(config_t& bconfig,
                                            std::string outputdir,
                                            ns_sequence const* seq = nullptr);

/**
 * @brief Driver for computing a NS solution including increasing resolution
 *
 * @tparam config_t Configurator type
 * @tparam Res_t Parameter_sequence for resolution
 * @param bconfig NS config file
 * @param resolution Resolution sequence
 * @param outputdir output location
 * @return int error code
 */
template <class config_t, class Res_t>
inline int ns_isotropic_diff_rot_driver(config_t& bconfig,
                                        Res_t& resolution,
                                        std::string outputdir,
                                        ns_sequence const* seq = nullptr);
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath

#include "ns_isotropic_diff_rot_driver_imp.cpp"