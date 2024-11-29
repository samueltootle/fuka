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
#include "Solvers/co_solver_utils.hpp"
#include "Solvers/sequences/ns_sequence.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/sequences/sequence_utilities.hpp"
#include "ns_3d_xcts_solver.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief Setup complete Config file for an NS sequence
 *
 * @tparam config_t Configurator type
 * @param seqconfig Sequence Config object
 * @param outputdir output location
 * @return config_t
 */
template <class config_t>
config_t ns_3d_xcts_sequence_setup(config_t& seqconfig, std::string outputdir);

/**
 * @brief Sequence driver for an NS solution
 *
 * @tparam ns_sequence Parameter_sequence for an NS sequence
 * @tparam Res_t Parameter_sequence for resolution
 * @tparam config_t Configurator type
 * @param seqconfig Sequence Config object
 * @param seq Parameter sequence
 * @param resolution Resolution sequence
 * @param outputdir output location
 */
template <class Res_t, class config_t>
config_t ns_3d_xcts_sequence(config_t& seqconfig,
                             ns_sequence const& seq,
                             Res_t& resolution,
                             std::string outputdir);

/**
 * @brief Driver to compute a stationary solution for a given resolution
 *
 * @tparam config_t Configurator type
 * @param bconfig BH config file
 * @param outputdir output location
 * @param seq optional sequence to enable sequenced NS solutions
 * @return int error code
 */
template <typename config_t>
int ns_3d_xcts_stationary_driver(config_t& bconfig,
                                 std::string outputdir,
                                 ns_sequence const* seq = nullptr);

/**
 * @brief Driver to compute a base solution for a given resolution to build on
 *
 * @tparam config_t Configurator type
 * @param bconfig BH config file
 * @param outputdir output location
 * @return int error code
 */
template <typename config_t>
int ns_3d_xcts_base_solution_driver(config_t& bconfig,
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
 * @param seq optional sequence to enable sequenced NS solutions
 * @return int error code
 */
template <class config_t, class Res_t>
inline int ns_3d_xcts_driver(config_t& bconfig,
                             Res_t& resolution,
                             std::string outputdir,
                             ns_sequence const* seq = nullptr);

/**
 * @brief Driver for computing a boosted NS solution based on binary parameters
 *
 * @tparam config_t Configurator type
 * @tparam Res_t Parameter_sequence for resolution
 * @param bconfig NS config object
 * @param resolution Resolution sequence
 * @param outputdir output location
 * @param binconfig Binary config object
 * @param bco Index of object in binary config, e.g. NODES::BCO1, NODES::BCO2
 * @return int
 */
template <typename config_t, class Res_t>
inline int ns_3d_xcts_binary_boost_driver(
    config_t& bconfig,
    Res_t& resolution,
    std::string outputdir,
    kadath_config_boost<BIN_INFO> binconfig,
    const size_t bco);
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
#include "ns_3d_xcts_driver.cpp"
