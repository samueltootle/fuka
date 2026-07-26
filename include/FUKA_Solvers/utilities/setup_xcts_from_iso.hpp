#pragma once
#include "Solvers/ns_3d_xcts/NS_ISO_to_XCTS_convert.hpp"
#include "Solvers/ns_isotropic/ns_isotropic_driver.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "fuka_path_tools.hpp"

namespace Kadath::FUKA_Solvers {
using namespace Kadath::FUKA_EOS;

template <typename config_t>
std::string solve_NS_ISO_from_XCTS_config(config_t& bconfig,
                                          ns_sequence const& seq) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string output_path = (bconfig.control(CONTROLS::SAVE_COS))
                                  ? get_cos_path()
                                  : bconfig.config_outputdir() + "/COs";
    if (rank == 0) {
        fs::create_directories(output_path);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    kadath_config_boost<BCO_ISO_NS_INFO> nsconfig;
    nsconfig.set_defaults();

    // Activate all controls that are active in the seq config
    for (int i = 0; i < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++i)
        nsconfig.seq_setting(i) = bconfig.seq_setting(i);

    // copy parameters from binary configuration
    for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
        nsconfig.set(i) = bconfig.set(i);
    nsconfig.set(BCO_PARAMS::DIM) = 2;

    for (int i = 0; i < CONTROLS::NUM_CONTROLS; ++i)
        nsconfig.control(i) = bconfig.control(i);

    for (int i = 0; i < STAGES::NUM_STAGES; ++i) {
        nsconfig.set_stage(i) = bconfig.set_stage(i);
    }
    nsconfig.set_stage(STAGES::UNIFORM_ROT) =
        bconfig.set_stage(STAGES::TOTAL_BC) ||
        bconfig.set_stage(STAGES::UNIFORM_ROT);
    nsconfig.set_stage(STAGES::TOTAL_BC) = false;

    // Tells the NS driver to initialize the numerical space and fields
    nsconfig.control(CONTROLS::SEQUENCES) = true;
    update_eos_parameters(bconfig, nsconfig);
    update_diffrot_parameters(bconfig, nsconfig);

    Parameter_sequence resolution("res", BCO_PARAMS::BCO_RES);
    resolution.set(9, 9, nsconfig(BCO_PARAMS::BCO_RES));

    nsconfig.set_filename("initns");
    nsconfig.set_outputdir(output_path);

    auto ns_iso_sol_config =
        ns_isotropic_sequence(nsconfig, seq, resolution, output_path);

    const std::string eos_type = bconfig.template eos<std::string>(EOSTYPE);
    nsconfig.control(CONTROLS::SEQUENCES) = false;
    EOS_Function_Dispatcher::dispatch<NS_ISO_to_XCTS_convert>(ns_iso_sol_config,
                                                              eos_type,
                                                              ns_iso_sol_config,
                                                              output_path);
    ns_iso_sol_config.set_filename("initns_xcts.info");
    ns_iso_sol_config.set_outputdir(output_path);
    MPI_Barrier(MPI_COMM_WORLD);

    return ns_iso_sol_config.config_filename_abs();
}
}  // namespace Kadath::FUKA_Solvers