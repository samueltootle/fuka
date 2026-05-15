#pragma once
#include "Solvers/ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "FUKA_Solvers/utilities/fuka_path_tools.hpp"



 /**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/
namespace Kadath::FUKA_Solvers {
/**
 * @brief
 * -Solve TOV solution based on binary Configurator input
 * -Update configurator based on TOV solution
 *
 * @tparam config_t: binary config type
 * @param[input] bconfig: binary Configurator file
 * @param[return] TOV solution filename
 */
template <typename config_t>
std::string solve_NS_for_binary_v1(config_t& bconfig, const size_t bco) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string output_path = (bconfig.control(CONTROLS::SAVE_COS))
                                  ? get_cos_path()
                                  : bconfig.config_outputdir() + "/COs";
    fs::create_directory(output_path);

    kadath_config_boost<BCO_NS_INFO> nsconfig;
    nsconfig.set_defaults();

    // Tells the NS driver to initialize the numerical space and fields
    nsconfig.control(CONTROLS::SEQUENCES) = true;

    // copy parameters from binary configuration
    for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
        nsconfig.set(i) = bconfig.set(i, bco);

    for (int i = 0; i < EOS_PARAMS::NUM_EOS_PARAMS; ++i)
        nsconfig.set_eos(i) = bconfig.set_eos(i, bco);

    for (int i = 0; i < CONTROLS::NUM_CONTROLS; ++i)
        nsconfig.control(i) = bconfig.control(i);

    for (int i = 0; i < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++i)
        nsconfig.seq_setting(i) = bconfig.seq_setting(i);

    Parameter_sequence resolution("res", BCO_PARAMS::BCO_RES);
    resolution.set(9, 9, nsconfig(BCO_PARAMS::BCO_RES));

    nsconfig.set_filename("initns");
    nsconfig.set_outputdir(output_path);

    if (bconfig.control(CONTROLS::USE_BOOSTED_CO))
        nsconfig.set_stage(STAGES::BIN_BOOST) = true;

    const int ninshells = bconfig(BCO_PARAMS::NINSHELLS, bco);
    const int nshells = bconfig(BCO_PARAMS::NSHELLS, bco);

    nsconfig(BCO_PARAMS::NINSHELLS) = 0;
    nsconfig(BCO_PARAMS::NSHELLS) = (bconfig.control(CONTROLS::CO_USE_SHELLS)) ? nshells : 0;

    int err = ns_3d_xcts_binary_boost_driver(nsconfig,
                                             resolution,
                                             output_path,
                                             bconfig,
                                             bco);
    MPI_Barrier(MPI_COMM_WORLD);

    // update binary parameters
    for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
        bconfig.set(i, bco) = nsconfig.set(i);
    // put shells back for the binary
    bconfig.set(BCO_PARAMS::NINSHELLS, bco) = ninshells;
    bconfig.set(BCO_PARAMS::NSHELLS, bco) = nshells;

    return nsconfig.config_filename_abs();
}
/** @}*/
}  // namespace Kadath::FUKA_Solvers