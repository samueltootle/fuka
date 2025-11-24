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


/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
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

#include "simple_calculations.cpp"