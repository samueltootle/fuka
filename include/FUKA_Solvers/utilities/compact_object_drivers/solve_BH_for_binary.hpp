#pragma once
#include "FUKA_Solvers/utilities/fuka_path_tools.hpp"
#include "Solvers/bh_3d_xcts/bh_3d_xcts_driver.hpp"

namespace Kadath::FUKA_Solvers {
template <typename config_t>
std::string solve_BH_for_binary(config_t& bconfig, const size_t bco) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string output_path = (bconfig.control(CONTROLS::SAVE_COS))
                                  ? get_cos_path()
                                  : bconfig.config_outputdir() + "/COs";
    fs::create_directory(output_path);

    kadath_config_boost<BCO_BH_INFO> bhconfig;
    bhconfig.set_defaults();

    // Tells the NS driver to initialize the numerical space and fields
    bhconfig.control(CONTROLS::SEQUENCES) = true;

    // copy parameters from binary configuration
    for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
        bhconfig.set(i) = bconfig.set(i, bco);

    for (int i = 0; i < CONTROLS::NUM_CONTROLS; ++i)
        bhconfig.control(i) = bconfig.control(i);

    for (int i = 0; i < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++i)
        bhconfig.seq_setting(i) = bconfig.seq_setting(i);

    Parameter_sequence resolution("res", BCO_PARAMS::BCO_RES);
    resolution.set(9, 9, bhconfig(BCO_PARAMS::BCO_RES));

    bhconfig.set_filename("initbh");
    bhconfig.set_outputdir(output_path);

    if (bconfig.control(CONTROLS::USE_BOOSTED_CO))
        bhconfig.set_stage(STAGES::BIN_BOOST) = true;

    const int nshells = bconfig(NSHELLS, bco);
    bhconfig(BCO_PARAMS::NSHELLS) =
        (bconfig.control(CONTROLS::CO_USE_SHELLS)) ? nshells : 0;

    int err = bh_3d_xcts_binary_boost_driver(bhconfig,
                                             resolution,
                                             output_path,
                                             bconfig,
                                             bco);
    MPI_Barrier(MPI_COMM_WORLD);

    // update binary parameters
    for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
        bconfig.set(i, bco) = bhconfig.set(i);
    bconfig.set(BCO_PARAMS::NSHELLS, bco) = nshells;

    MPI_Barrier(MPI_COMM_WORLD);
    return bhconfig.config_filename_abs();
}
}  // namespace Kadath::FUKA_Solvers