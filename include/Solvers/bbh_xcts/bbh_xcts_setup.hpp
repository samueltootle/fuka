#pragma once
#include <array>
#include <string>
#include "bbh_xcts_solver.hpp"
#include "mpi.h"

/**
 * \addtogroup BBH_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

/**
 * bbh_xcts_setup_bin
 *
 * Ensure proper binary settings
 *
 * @param[input] bconfig: BBH Configurator file
 */
template <class config_t>
void bbh_xcts_setup_bin_config(config_t& bconfig);

/**
 * bbh_xcts_setup_space
 *
 * Create the numerical space and initial guess for the BBH by
 * obtaining the isolated BH solutions, updating the BBH
 * config file, and importing the BH solutions into the BBH space.
 *
 * @param[input] bconfig: BBH Configurator file
 */
template <class config_t>
void bbh_xcts_setup_space(config_t& bconfig);

/**
 * bbh_xcts_superimposed_import
 *
 * Flow control to initialize the binary space based on
 * superimposing two isolated compact object solutions.
 * These are simply read from file here before being passed
 * to bbh_xcts_setup_boosted_3d for computing the initial
 * guess
 *
 * @param[input] bconfig: binary configurator
 * @param[input] BHfilenames: array of filenames for isolated solutions
 */
template <class config_t>
void bbh_xcts_superimposed_import(config_t& bconfig,
                                  std::array<std::string, 2> BHfilenames);

/**
 * @brief Generate superimposed guess from 3D isolated solutions
 *
 * @param BH1config First BH solution
 * @param BH2config Second BH solution
 * @param bconfig Binary Config
 */
inline void bbh_xcts_setup_boosted_3d(
    kadath_config_boost<BCO_BH_INFO>& BH1config,
    kadath_config_boost<BCO_BH_INFO>& BH2config,
    kadath_config_boost<BIN_INFO>& bconfig);
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath

#include "bbh_xcts_setup.cpp"
