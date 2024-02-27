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
#include "mpi.h"
#include "bhns_xcts_solver.hpp"
#include <array>
#include <string>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

/**
 * \addtogroup BHNS_XCTS
 * \ingroup FUKA
 * @{*/

// using namespace Kadath::FUKA_Config;
namespace Kadath {
namespace FUKA_Solvers {

/**
 * bhns_xcts_setup_bin
 *
 * Ensure proper binary settings
 *
 * @param[input] bconfig: BHNS Configurator file
 */
template<class config_t>
inline void bhns_xcts_setup_bin_config(config_t& bconfig);

/**
 * bhns_xcts_setup_space
 *
 * Create the numerical space and initial guess for the BHNS by
 * obtaining the isolated BH solutions, updating the BHNS
 * config file, and importing the BH solutions into the BHNS space.
 *
 * @param[input] bconfig: BHNS Configurator file
 */
template<class config_t>
void bhns_xcts_setup_space (config_t& bconfig);

/**
 * bhns_xcts_superimposed_import
 *
 * Flow control to initialize the binary space based on
 * superimposing two isolated compact object solutions.
 * These are simply read from file here before being passed
 * to bhns_xcts_setup_boosted_3d for computing the initial
 * guess
 *
 * @param[input] bconfig: binary configurator
 * @param[input] co_filenames: array of filenames for isolated solutions
 */
template<class config_t>
void bhns_xcts_superimposed_import(config_t& bconfig,std::array<std::string, 2> co_filenames);

/**
 * @brief Generate superimposed guess from 3D isolated solutions
 * 
 * @param NSconfig NS solution
 * @param BHconfig BH solution
 * @param bconfig Binary Config
 */
template<typename eos_t>
inline void bhns_setup_boosted_3d(
  kadath_config_boost<BCO_NS_INFO>& NSconfig, 
  kadath_config_boost<BCO_BH_INFO>& BHconfig,
  kadath_config_boost<BIN_INFO>& bconfig);
/** @}*/
}}
#include "bhns_xcts_setup.cpp"
