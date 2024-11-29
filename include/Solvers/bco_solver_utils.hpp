/*
 * Copyright 2022
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
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "EOS/standalone/tov.hh"
#include "bco_utilities.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

/**
 * solve_NS_from_binary
 *
 * -Solve TOV solution based on binary Configurator input
 * -Update configurator based on TOV solution
 *
 * @tparam config_t: binary config type
 * @param[input] bconfig: binary Configurator file
 * @param[return] TOV solution filename
 */

template <typename config_t>
std::string solve_NS_from_binary(config_t& bconfig, const size_t bco);

/**
 * solve_BH_from_binary
 *
 * -Solve BH solution based on binary Configurator input
 * -Update configurator based on BH solution
 *
 * @tparam config_t: binary config type
 * @param[input] bconfig: binary Configurator file
 * @param[return] TOV solution filename
 */

template <typename config_t>
std::string solve_BH_from_binary(config_t& bconfig, const size_t bco);

/**
 * check_dist
 *
 * Check to ensure separation distance is reasonable
 *
 * @param[input] dist: separation distance
 * @param[input] M1: Gravitational Mass of object 1
 * @param[input] M2: Gravitational Mass of object 2
 * @param[input] garbage_factor: factor used to set garbage distance
 */
inline void check_dist(double dist,
                       double M1,
                       double M2,
                       double garbage_factor = 2.5);
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
#include "bco_solver_utils_imp.cpp"