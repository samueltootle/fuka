#pragma once
#include "Solvers/co_solver_utils.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/sequences/sequence_utilities.hpp"
#include "bh_3d_xcts_solver.hpp"

/**
 * \addtogroup BH_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief Setup complete Config file for a BH sequence
 *
 * @tparam config_t Configurator type
 * @param seqconfig Sequence Config object
 * @param outputdir output location
 * @return config_t
 */
template <class config_t>
config_t bh_3d_xcts_sequence_setup(config_t& seqconfig, std::string outputdir);

/**
 * @brief Sequence driver for a BH solution
 *
 * @tparam Seq_t Parameter_sequence for a BH sequence
 * @tparam Res_t Parameter_sequence for resolution
 * @tparam config_t Configurator type
 * @param seqconfig Sequence Config object
 * @param seq Parameter sequence
 * @param resolution Resolution sequence
 * @param outputdir output location
 */
template <class Seq_t, class Res_t, class config_t>
config_t bh_3d_xcts_sequence(config_t& seqconfig,
                             Seq_t const& seq,
                             Res_t const& resolution,
                             std::string outputdir);

/**
 * @brief Driver to compute a stationary solution for a given resolution
 *
 * @tparam config_t Configurator type
 * @param bconfig BH config file
 * @param outputdir output location
 * @return int error code
 */
template <typename config_t>
int bh_3d_xcts_stationary_driver(config_t& bconfig, std::string outputdir);

/**
 * @brief Driver to compute linear boosted BH solutions in the XY plane
 * Note: Requires a stationary solution
 *
 * @tparam config_t Configurator type
 * @param bconfig BH Config file
 * @param outputdir output location
 * @return int error code
 */
template <typename config_t>
int bh_3d_xcts_linear_boost_driver(config_t& bconfig, std::string outputdir);

/**
 * @brief Driver for computing a BH solution including increasing resolution
 *
 * @tparam config_t Configurator type
 * @tparam Res_t Parameter_sequence for resolution
 * @param bconfig BH config file
 * @param resolution Resolution sequence
 * @param outputdir output location
 * @return int error code
 */
template <class config_t, class Res_t>
inline int bh_3d_xcts_driver(config_t& bconfig,
                             Res_t& resolution,
                             std::string outputdir);

/**
 * @brief Driver for computing a boosted BH solution based on binary parameters
 *
 * @tparam config_t Configurator type
 * @tparam Res_t Parameter_sequence for resolution
 * @param bconfig BH config object
 * @param resolution Resolution sequence
 * @param outputdir output location
 * @param binconfig Binary config object
 * @param bco Index of object in binary config, e.g. NODES::BCO1, NODES::BCO2
 * @return int
 */
template <typename config_t, class Res_t>
inline int bh_3d_xcts_binary_boost_driver(
    config_t& bconfig,
    Res_t& resolution,
    std::string outputdir,
    kadath_config_boost<BIN_INFO> binconfig,
    const size_t bco);
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
#include "bh_3d_xcts_driver.cpp"