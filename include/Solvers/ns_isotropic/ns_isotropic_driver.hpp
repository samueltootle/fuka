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
#include "norot/ns_isotropic_norot_solver.hpp"
#include "norot/ns_isotropic_norot_driver.hpp"
#include "uniform_rot/ns_isotropic_uniform_rot_driver.hpp"
#include "Solvers/co_solver_utils.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/sequences/sequence_utilities.hpp"

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
template<class config_t>
config_t ns_isotropic_sequence_setup (config_t & seqconfig, std::string outputdir);


/**
 * @brief Sequence driver for an NS solution
 * 
 * @tparam Seq_t Parameter_sequence for an NS sequence
 * @tparam Res_t Parameter_sequence for resolution
 * @tparam config_t Configurator type
 * @param seqconfig Sequence Config object
 * @param seq Parameter sequence
 * @param resolution Resolution sequence
 * @param outputdir output location
 */
template<class Seq_t, class Res_t, class config_t>
config_t ns_isotropic_sequence (config_t & seqconfig, 
                          Seq_t const & seq,
                          Res_t const & resolution,
                          std::string outputdir);

/**
 * @brief Driver to compute a stationary solution for a given resolution
 * 
 * @tparam config_t Configurator type
 * @param bconfig BH config file
 * @param outputdir output location
 * @return int error code
 */
template<typename config_t>
int ns_isotropic_stationary_driver (config_t& bconfig, std::string outputdir);

/**
 * @brief Driver to compute a stationary solution for a given resolution
 * 
 * @tparam config_t Config file type
 * @param bconfig NS config file
 * @param outputdir directory to store solutions in
 * @return int error code
 */
template<typename config_t>
int ns_isotropic_base_solution_driver (config_t& bconfig, std::string outputdir);

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
template<class config_t, class Res_t>
inline int ns_isotropic_driver (config_t& bconfig, Res_t& resolution, std::string outputdir, Parameter_sequence<BCO_PARAMS> const * seq=nullptr);
/** @}*/
}}
#include "ns_isotropic_driver_imp.cpp"