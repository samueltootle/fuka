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
#include "kadath_adapted_bh.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief Set the initial guess for a compact object
 *
 * @tparam s_type NODES::NS or NODES::BH
 * @tparam config_t Config object type
 * @param bconfig Config object
 */
template <std::size_t s_type, typename config_t>
void setup_co(config_t& bconfig);

template <typename config_t>
void setup_ns_3d_xcts(config_t& bconfig, size_t mass_fixing_idx);

/**
 * @brief Set the initial guess for 2D Neutron star in
 * isotropic coordinates based on
 * arxiv.org:1003.5015
 *
 * @tparam eos_t EOS wrapper type
 */
template <class eos_t>
struct setup_2dns_isotropic_functor;

/**
 * write_bh_init_setup_tofile_XCTS
 *
 * Takes a numerical space and configuration and constructs the initial guess
 * for the fields associated with the XCTS formalism and writes everything to
 * file at the end.
 *
 * @tparam config_t configuration file type
 * @param [input] space numerical space
 * @param [input] bconfig the configuration file
 */
template <typename config_t>
void write_bh_init_setup_tofile_XCTS(Space_adapted_bh& space,
                                     config_t& bconfig);

/**
 * write_ns_init_setup_tofile_XCTS
 *
 * Takes a numerical space and configuration and constructs the initial guess
 * for the fields associated with the XCTS formalism and writes everything to
 * file at the end.
 *
 * @tparam config_t configuration file type
 * @tparam tov_uptr unique_ptr to the given 1D tov solution type
 * @param [input] space numerical space
 * @param [input] bconfig the configuration file
 */
template <typename tov_t, typename config_t>
void write_ns_init_setup_tofile_XCTS(Space_spheric_adapted& space,
                                     config_t& bconfig,
                                     tov_t& tov);

/**
 * write_ns2d_isotropic_init_setup_tofile
 *
 * Takes a numerical space and configuration and constructs the initial guess
 * for the fields associated with a Neutron star in isotropic coordinates and
 * writes everything to file at the end.
 *
 * @tparam config_t configuration file type
 * @tparam tov_uptr unique_ptr to the given 1D tov solution type
 * @param [input] space numerical space
 * @param [input] bconfig the configuration file
 */
template <typename tov_t, typename config_t>
void write_ns2d_isotropic_init_setup_tofile(Space_polar_adapted& space,
                                            config_t& bconfig,
                                            tov_t& tov);

/**
 * setup_ns_config_from_TOV
 *
 * For a given EOS, solve the 1D tov equations in order to find the relevant TOV
 * quantities for the initial guess. From this initial guess, we extract NC, HC,
 * and Radii quantities needed for setting up the initial guess and Space.
 *
 * @tparam eos_t type of EOS used
 * @param [input] bconfig Input configuration file to update
 * @return tov the 1D tov solution
 */
template <typename eos_t, typename config_t>
auto setup_ns_config_from_TOV(config_t& bconfig,
                              size_t mass_fixing_idx = ::Kadath::FUKA_Config::BCO_PARAMS::MADM);

/**
 * setup_interpolator_from_TOV
 *
 * For a given TOV solution, setup a linear interpolator for the relevant
 * quantities (R, PSI, ALPHA, RHO)
 *
 * @tparam tov_uptr unique_ptr to the given 1D tov solution type
 * @param [input] tov 1D tov solution
 * @return ltp the linear interpolator
 */
template <typename tov_t>
auto setup_interpolator_from_TOV(tov_t& tov);

/**
 * setup_KerrSchild_BH
 *
 * Initialize a new KerrSchild space and initial guess before writing to file
 * @tparam space_t: numerical space type
 * @tparam config_t: Configurator type
 * @param[input] bconfig: Configurator file
 */
template <typename space_t, typename config_t>
void setup_KerrSchild_BH(config_t& bconfig);

/**
 * setup_KerrSchild_BH
 *
 * Write initial KS setup to file
 * @tparam space_t: numerical space type
 * @tparam config_t: Configurator type
 * @param[input] space: numerical space
 * @param[input] bconfig: Configurator file
 */
template <class space_t, class config_t>
void write_KerrSchild_bh_init_setup_tofile_XCTS(space_t& space,
                                                config_t& bconfig);
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
#include "co_solver_utils_imp.cpp"
