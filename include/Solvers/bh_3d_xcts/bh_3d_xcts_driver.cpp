/**
 * \addtogroup BH_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template <class config_t>
config_t bh_3d_xcts_sequence_setup(config_t& seqconfig, std::string outputdir) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    config_t bconfig = generate_sequence_config(seqconfig, outputdir);

    if (rank == 0)
        bconfig.write_config();
    return bconfig;
}

template <class Seq_t, class Res_t, class config_t>
config_t bh_3d_xcts_sequence(config_t& seqconfig,
                             Seq_t const& seq,
                             Res_t const& resolution,
                             std::string outputdir) {
    int rank = 0;
    int exit_status = EXIT_SUCCESS;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Initialize sequence variables
    auto sequence_var_indices = seq.get_indices();
    auto resolution_indices = resolution.get_indices();

    auto const& dx = seq.step_size();

    // Initialize full configurator
    config_t base_config = bh_3d_xcts_sequence_setup(seqconfig, outputdir);
    base_config.set(resolution_indices) = resolution.init();

    // Declare here such that it is within scope to be returned
    config_t bconfig{};

    std::array<bool, NUM_STAGES> stage_enabled = base_config.return_stages();
    auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

#ifdef DEBUG
    std::cout << seq << std::endl;
    std::cout << resolution << std::endl;
#endif
    auto single_seq = [&](auto val) {
        bconfig = base_config;
        if (seq.is_set())
            bconfig.set(sequence_var_indices) = val;

        if (bconfig.control(CONTROLS::SEQUENCES)) {
            if (rank == 0) {
                setup_3d_BH_xcts(bconfig);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // make sure all ranks have the same config
            bconfig.open_config();
            MPI_Barrier(MPI_COMM_WORLD);
        }

        exit_status = bh_3d_xcts_driver(bconfig, resolution, outputdir);
        return exit_status;
    };

    // Loop if a valid sequence is set
    if (seq.is_set()) {
        for (auto val = seq.init(); seq.loop_condition(val); val += dx) {
            exit_status = single_seq(val);
        }
    } else {
        exit_status = single_seq(0);
    }
    return bconfig;
}

template <typename config_t>
int bh_3d_xcts_stationary_driver(config_t& bconfig, std::string outputdir) {
    int exit_status = RELOAD_FILE;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (std::abs(bconfig(BCO_PARAMS::CHI)) > 0.84) {
        std::cerr << "Warning: Solutions above |chi| > 0.84 "
                     "are not necessarily stable solutions\n";
    }

    // In the event we wish to solve for a highly spinning solution
    // we need to do an initial slow rotating solution before going to
    // faster rotations otherwise the solution will diverge.
    bconfig.control(CONTROLS::ITERATIVE_CHI) =
        bconfig.control(CONTROLS::SEQUENCES) &&
        std::abs(bconfig(BCO_PARAMS::CHI)) > 0.5;

    // We need to store the final desired CHI in case of iterative chi
    bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI) = bconfig(BCO_PARAMS::CHI);
    // Lower chi in case of iterative chi
    bconfig(BCO_PARAMS::CHI) =
        (bconfig.control(CONTROLS::ITERATIVE_CHI))
            ? std::copysign(0.5, bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI))
            : bconfig(BCO_PARAMS::CHI);

    // make sure BH directory exists for outputs
    if (outputdir == "./") {
        fs::path cwd = fs::current_path();
        outputdir = cwd.string();
    }
    if (rank == 0)
        std::cout << "Solutions will be stored in: " << outputdir << "\n"
                  << "Directory will be created if it doesn't exist.\n";
    fs::create_directory(outputdir);

    std::string spacein = bconfig.space_filename();
    if (!fs::exists(spacein)) {
        // mainly for debugging MPI bugs
        if (rank == 0) {
            std::cerr << "File: " << spacein << " not found.\n\n";
        } else {
            std::cerr << "File: " << spacein
                      << " not found for another rank.\n\n";
        }
        std::_Exit(EXIT_FAILURE);
    }

    while (exit_status == RELOAD_FILE) {
        spacein = bconfig.space_filename();
        // just so you really know
        if (rank == 0) {
            std::cout << "Config File: " << bconfig.config_filename_abs()
                      << std::endl
                      << "Fields File: " << spacein << std::endl
                      << bconfig << std::endl;
        }
        FILE* ff1 = fopen(spacein.c_str(), "r");
        if (ff1 == NULL) {
            // mainly for debugging MPI bugs
            std::cerr << spacein.c_str() << " failed to open for rank " << rank
                      << "\n";
            std::_Exit(EXIT_FAILURE);
        }
        Space_adapted_bh space(ff1);
        Scalar conf(space, ff1);
        Scalar lapse(space, ff1);
        Vector shift(space, ff1);
        fclose(ff1);
        Base_tensor basis(space, CARTESIAN_BASIS);

        if (outputdir != "")
            bconfig.set_outputdir(outputdir);

        bh_3d_xcts_solver<decltype(bconfig), decltype(space)> bh_solver(bconfig,
                                                                        space,
                                                                        basis,
                                                                        conf,
                                                                        lapse,
                                                                        shift);

        exit_status = bh_solver.solve();

        MPI_Barrier(MPI_COMM_WORLD);
    }
    return exit_status;
}

template <typename config_t>
int bh_3d_xcts_linear_boost_driver(config_t& bconfig, std::string outputdir) {
    int exit_status = RELOAD_FILE;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string spacein = bconfig.space_filename();
    if (!fs::exists(spacein)) {
        // mainly for debugging MPI bugs
        if (rank == 0) {
            std::cerr << "File: " << spacein << " not found.\n\n";
        } else {
            std::cerr << "File: " << spacein
                      << " not found for another rank.\n\n";
        }
        std::_Exit(EXIT_FAILURE);
    }

    while (exit_status == RELOAD_FILE) {
        spacein = bconfig.space_filename();
        // just so you really know
        if (rank == 0) {
            std::cout << "Config File: " << bconfig.config_filename_abs()
                      << std::endl
                      << "Fields File: " << spacein << std::endl
                      << bconfig << std::endl;
        }
        FILE* ff1 = fopen(spacein.c_str(), "r");
        if (ff1 == NULL) {
            // mainly for debugging MPI bugs
            std::cerr << spacein.c_str() << " failed to open for rank " << rank
                      << "\n";
            std::_Exit(EXIT_FAILURE);
        }
        Space_adapted_bh space(ff1);
        Scalar conf(space, ff1);
        Scalar lapse(space, ff1);
        Vector shift(space, ff1);
        fclose(ff1);
        Base_tensor basis(space, CARTESIAN_BASIS);

        if (outputdir != "")
            bconfig.set_outputdir(outputdir);

        bh_3d_xcts_solver<decltype(bconfig), decltype(space)> bh_solver(bconfig,
                                                                        space,
                                                                        basis,
                                                                        conf,
                                                                        lapse,
                                                                        shift);
        bh_solver.set_solver_stage(STAGES::LINBOOST);
        exit_status = bh_solver.von_Neumann_stage("LINBOOST");

        MPI_Barrier(MPI_COMM_WORLD);
    }
    return exit_status;
}

template <class config_t, class Res_t>
inline int bh_3d_xcts_driver(config_t& bconfig,
                             Res_t& resolution,
                             std::string outputdir) {
    int exit_status = RELOAD_FILE;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string spacein = bconfig.space_filename();
    if (!fs::exists(spacein)) {
        // mainly for debugging MPI bugs
        if (rank == 0) {
            std::cerr << "File: " << spacein << " not found.\n\n";
        } else {
            std::cerr << "File: " << spacein
                      << " not found for another rank.\n\n";
        }
        std::_Exit(EXIT_FAILURE);
    }

    double const xboost = bconfig.set(BCO_PARAMS::BVELX);
    double const yboost = bconfig.set(BCO_PARAMS::BVELY);
    bool res_inc = (resolution.final() != resolution.init());
    auto resolution_indices = resolution.get_indices();
    auto const& final_res = resolution.final();

    std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
    auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

    // These need to be zero prior to computing
    // a stationary solution
    bconfig.set(BCO_PARAMS::BVELX) = 0.;
    bconfig.set(BCO_PARAMS::BVELY) = 0.;

    exit_status = bh_3d_xcts_stationary_driver(bconfig, outputdir);

    if (stage_enabled[STAGES::LINBOOST]) {
        bconfig.set(BCO_PARAMS::BVELX) = xboost;
        bconfig.set(BCO_PARAMS::BVELY) = yboost;
        exit_status = bh_3d_xcts_linear_boost_driver(bconfig, outputdir);
    }

    auto regrid = [&]() {
        std::string fname{"bh_regrid"};

        if (rank == 0)
            exit_status = bh_3d_xcts_regrid(bconfig, fname);
        MPI_Barrier(MPI_COMM_WORLD);
        bconfig.set_filename(fname);
        bconfig.open_config();

        // Deactivate to avoid iterative CHI
        bconfig.control(CONTROLS::SEQUENCES) = false;

        // Ensure last stage is activated to resolve
        stage_enabled[last_stage_idx] = true;
    };

    while (res_inc) {
        // iterative res increase
        if (bconfig(BCO_PARAMS::BCO_RES) + 2 >= final_res) {
            bconfig.set(BCO_PARAMS::BCO_RES) = final_res;
            res_inc = false;
        } else {
            bconfig.set(BCO_PARAMS::BCO_RES) += 2;
        }
        regrid();

        // Rerun with new grid
        if (stage_enabled[STAGES::LINBOOST]) {
            exit_status = bh_3d_xcts_linear_boost_driver(bconfig, outputdir);
        } else {
            exit_status = bh_3d_xcts_stationary_driver(bconfig, outputdir);
        }
    }
    return exit_status;
}

template <typename config_t, class Res_t>
inline int bh_3d_xcts_binary_boost_driver(
    config_t& bconfig,
    Res_t& resolution,
    std::string outputdir,
    kadath_config_boost<BIN_INFO> binconfig,
    const size_t bco) {
    int exit_status = RELOAD_FILE;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int final_res = bconfig(BCO_PARAMS::BCO_RES);
    bool res_inc = (bconfig.seq_setting(SEQ_SETTINGS::INIT_RES) < final_res);
    bconfig.set(BCO_PARAMS::BCO_RES) = resolution.init();

    // Obtain low res, stationary solution
    auto tmp_res(resolution);
    tmp_res.set(resolution.init(), resolution.init(), resolution.init());
    Parameter_sequence<BCO_PARAMS> tmp_seq{};
    bconfig = bh_3d_xcts_sequence(bconfig, tmp_seq, tmp_res, outputdir);

    while (exit_status == RELOAD_FILE || exit_status == RUN_BOOST) {
        auto spacein = bconfig.space_filename();
        // just so you really know
        if (rank == 0) {
            std::cout << "Config File: " << bconfig.config_filename_abs()
                      << std::endl
                      << "Fields File: " << spacein << std::endl
                      << bconfig << std::endl;
        }
        FILE* ff1 = fopen(spacein.c_str(), "r");
        if (ff1 == NULL) {
            // mainly for debugging MPI bugs
            std::cerr << spacein.c_str() << " failed to open for rank " << rank
                      << "\n";
            std::_Exit(EXIT_FAILURE);
        }
        Space_adapted_bh space(ff1);
        Scalar conf(space, ff1);
        Scalar lapse(space, ff1);
        Vector shift(space, ff1);
        fclose(ff1);
        Base_tensor basis(space, CARTESIAN_BASIS);

        bh_3d_xcts_solver<decltype(bconfig), decltype(space)> bh_solver(bconfig,
                                                                        space,
                                                                        basis,
                                                                        conf,
                                                                        lapse,
                                                                        shift);
        bh_solver.set_solver_stage(STAGES::BIN_BOOST);
        exit_status = bh_solver.binary_boost_stage(binconfig, bco);

        if (res_inc && exit_status != RELOAD_FILE) {
            res_inc = false;
            bconfig.set(BCO_PARAMS::BCO_RES) = final_res;
            if (rank == 0)
                exit_status = bh_3d_xcts_regrid(bconfig, "initbh");
            bconfig.set_filename("initbh");

            MPI_Barrier(MPI_COMM_WORLD);
            bconfig.open_config();
            exit_status = RELOAD_FILE;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return exit_status;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
