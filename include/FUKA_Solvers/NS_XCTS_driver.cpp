#include "FUKA_Solvers/utilities/compact_object_initializers/setup_3dns_xcts.hpp"
#include "FUKA_Solvers/utilities/setup_xcts_from_iso.hpp"
#include "NS_XCTS.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath::FUKA_Solvers {

template <class eos_t>
struct launch_ns_xcts_solver {
    template <class config_t>
    int operator()(const int rank,
                   config_t& bconfig,
                   std::string outputdir,
                   Parameter_sequence<BCO_PARAMS>& resolution,
                   ns_sequence const& seq) {

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
            std::stringstream ss;
            ss << spacein.c_str() << " failed to open for rank " << rank
               << "\n";
            throw std::runtime_error(ss.str().c_str());
        }
        int exit_status = EXIT_SUCCESS;
        auto launch = [](auto& solver,
                         bool ignore_resinc = false,
                         bool ignore_seq = false) {
            int exit_status = EXIT_SUCCESS;
            do {
                // initial solution
                solver.setup_syst();
                solver.do_newton();

                if (!ignore_resinc) {
                    // Make sure final solution uses optimal domain decomposition
                    solver.regrid();

                    // resolve at current resolution
                    solver.setup_syst();
                    solver.do_newton();
                }

                // Obtain final resolution for the desired
                // solution or the first solution in a sequence
                // All remaining sequences will be computed
                // at the final resolution only
                // Note: only occurs if the last stage is NOROT_BC
                while (!ignore_resinc && solver.increment_resolution() &&
                       (exit_status != EXIT_FAILURE)) {
                    // regrid to new resolution
                    solver.regrid();

                    // initial solution
                    solver.setup_syst();
                    exit_status = solver.do_newton();
                }
            } while (!ignore_seq && solver.increment_seq());
            return exit_status;
        };

        std::array<bool, NUM_STAGES> const stage_enabled =
            bconfig.return_stages();
        if (stage_enabled[STAGES::NOROT_BC] && (exit_status != EXIT_FAILURE)) {
            NS_XCTS_NOROT<eos_t> norot_solver(&bconfig,
                                              seq,
                                              resolution,
                                              outputdir,
                                              rank);
            launch(norot_solver);
        }

        if (stage_enabled[STAGES::UNIFORM_ROT] &&
            (exit_status != EXIT_FAILURE)) {
            NS_XCTS_UNIFORM_ROT<eos_t> uniformrot_solver(&bconfig,
                                                         seq,
                                                         resolution,
                                                         outputdir,
                                                         rank);
            auto const& spinup(uniformrot_solver.get_spinup());
            bool check = (spinup && spinup->is_set() && spinup->is_varying());
            do {
                launch(uniformrot_solver, check, check);
            } while (uniformrot_solver.increment_spin());
            launch(uniformrot_solver);
        } else if (stage_enabled[STAGES::DIFF_ROT] &&
                   bconfig.control(CONTROLS::SEQUENCES) &&
                   (exit_status != EXIT_FAILURE)) {
            auto tmp_seq(seq);
            tmp_seq.set_spin_idx(BCO_PARAMS::CHI);
            tmp_seq.set_spin_val(0.1);
            bconfig.set(BCO_PARAMS::CHI) = tmp_seq.spin_val();
            NS_XCTS_UNIFORM_ROT<eos_t> uniformrot_solver(&bconfig,
                                                         tmp_seq,
                                                         resolution,
                                                         outputdir,
                                                         rank);
            auto const& spinup(uniformrot_solver.get_spinup());
            bool check = (spinup && spinup->is_set() && spinup->is_varying());
            do {
                launch(uniformrot_solver, check, check);
            } while (uniformrot_solver.increment_spin());
        }

        if (stage_enabled[STAGES::DIFF_ROT] && (exit_status != EXIT_FAILURE)) {
            NS_XCTS_DIFF_ROT<eos_t> diffrot_solver(&bconfig,
                                                   seq,
                                                   resolution,
                                                   outputdir,
                                                   rank);
            auto const& spinup(diffrot_solver.get_spinup());
            bool check = (spinup && spinup->is_set() && spinup->is_varying());
            do {
                // launch(uniformrot_solver, spinup.is_set());
                // FIXME? launch(diffrot_solver, false, check); ????
                launch(diffrot_solver, true, true);
            } while (diffrot_solver.increment_spin());
            // FIXME? launch(diffrot_solver);
            launch(diffrot_solver, false, true);
        }
        return exit_status;
    };
};

inline int ns_xcts_driver(NS_XCTS_BASE::base_config_t& bconfig,
                          ns_sequence const& seq,
                          Parameter_sequence<BCO_PARAMS>& resolution,
                          std::string outputdir);

inline NS_XCTS_BASE::base_config_t ns_xcts_sequence_setup(
    NS_XCTS_BASE::base_config_t& seqconfig,
    std::string outputdir) {
    using config_t = NS_XCTS_BASE::base_config_t;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    config_t bconfig = generate_sequence_config(seqconfig, outputdir);

    update_eos_parameters(seqconfig, bconfig);
    update_diffrot_parameters(seqconfig, bconfig);

    if (rank == 0)
        bconfig.write_config();
    return bconfig;
}

inline int ns_xcts_solver_driver(NS_XCTS_BASE::base_config_t& bconfig,
                                 ns_sequence const& seq,
                                 Parameter_sequence<BCO_PARAMS>& resolution,
                                 std::string outputdir) {
    int exit_status = RELOAD_FILE;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // make sure NS directory exists for outputs
    if (outputdir == "./") {
        std::filesystem::path cwd = std::filesystem::current_path();
        outputdir = cwd.string();
    }

    const std::string eos_type =
        bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);
    using namespace Kadath::FUKA_EOS;
    exit_status =
        EOS_Function_Dispatcher::dispatch<launch_ns_xcts_solver>(bconfig,
                                                                 eos_type,
                                                                 rank,
                                                                 bconfig,
                                                                 outputdir,
                                                                 resolution,
                                                                 seq);

    MPI_Barrier(MPI_COMM_WORLD);
    return exit_status;
}

int ns_xcts_seq_driver(NS_XCTS_BASE::base_config_t& seqconfig,
                       ns_sequence const& seq,
                       Parameter_sequence<BCO_PARAMS>& resolution,
                       std::string const outputdir) {
    using config_t = NS_XCTS_BASE::base_config_t;
    int rank = 0, exit_status = EXIT_SUCCESS;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Ensure fixed values are initialized
    seqconfig.set(seq.mass_idx()) = seq.mass_val();

    // Initialize sequence variables
    auto sequence_var_indices = seq.get_indices();
    auto resolution_indices = resolution.get_indices();

    // Initialize full configurator
    config_t base_config = ns_xcts_sequence_setup(seqconfig, outputdir);
    base_config.set(resolution_indices) = resolution.init();

    // in the event the user wants an MADM > MTOV
    const double final_MADM = (seq.mass_idx() == BCO_PARAMS::MADM)
                                  ? base_config(BCO_PARAMS::MADM)
                                  : std::nan("1");
    base_config.control(CONTROLS::ITERATIVE_M) = false;
    config_t bconfig{base_config};

    if (!(base_config.control(CONTROLS::SEQUENCES) ||
          base_config.control(CONTROLS::RESOLVE))) {
        base_config = seqconfig;
        base_config.set(BCO_PARAMS::CHI) = 0;
    }

    auto mass_fixing = seq.mass_idx();

    // Get last enabled to make for safety check later
    std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
    auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

    if (bconfig.control(CONTROLS::SEQUENCES) ||
        bconfig.control(CONTROLS::RESOLVE)) {
        if (rank == 0) {
            // setup_co<NODES::NS>(bconfig);
            const std::string eos_type =
                bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);
            using namespace Kadath::FUKA_EOS;
            EOS_Function_Dispatcher::dispatch<setup_3dns_xcts_functor>(
                bconfig,
                eos_type,
                bconfig,
                mass_fixing);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // make sure all ranks have the same config
        bconfig.open_config();
        MPI_Barrier(MPI_COMM_WORLD);
        bconfig.control(CONTROLS::ITERATIVE_M) =
            !std::isnan(final_MADM) &&
            (std::fabs(1. - bconfig(BCO_PARAMS::MADM) / final_MADM) > 1e-3);

        if (bconfig.control(CONTROLS::ITERATIVE_M) &&
            last_stage_idx == STAGES::NOROT_BC) {
            std::stringstream ss;
            ss << "Cannot solve TOV for Madm = " << final_MADM
               << " without spin.\n";
            throw std::runtime_error(ss.str().c_str());
        }
    }
    // Get non-rotating solution for the given mass or TOV mass if
    // bconfig.control(CONTROLS::ITERATIVE_M)
    ns_xcts_solver_driver(bconfig, seq, resolution, outputdir);

    // Update config such that the next solving round uses
    // the final ADM mass and spin if applicable
    if (bconfig.control(CONTROLS::ITERATIVE_M)) {
        bconfig(BCO_PARAMS::MADM) = final_MADM;
    }
    bconfig.control(CONTROLS::SEQUENCES) = false;

    seqconfig = bconfig;

    return EXIT_SUCCESS;
}

inline int launch_final_stage_driver(NS_XCTS_BASE::base_config_t& bconfig,
                                     ns_sequence const& seq,
                                     Parameter_sequence<BCO_PARAMS>& resolution,
                                     std::string outputdir) {
    using config_t = NS_XCTS_BASE::base_config_t;
    using Res_t = Parameter_sequence<BCO_PARAMS>;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
    auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

    std::function<int(config_t&, ns_sequence const&, Res_t&, std::string)>
        final_stage_driver;
    if (rank == 0)
        std::cout << "Last stage: " << last_stage << '\n';
    if (bconfig.control(CONTROLS::SEQUENCES)) {
        final_stage_driver = &ns_xcts_seq_driver;
    } else {
        final_stage_driver = &ns_xcts_solver_driver;
    }
    return final_stage_driver(bconfig, seq, resolution, outputdir);
}

inline int ns_xcts_driver(NS_XCTS_BASE::base_config_t& bconfig,
                          ns_sequence const& seq,
                          Parameter_sequence<BCO_PARAMS>& resolution,
                          std::string outputdir) {
    int exit_status = EXIT_SUCCESS;
    auto resolution_indices = resolution.get_indices();
    bconfig.set(resolution_indices) = resolution.init();

    if (bconfig.control(CONTROLS::SEQUENCES)) {
        std::string initial_guess_filename =
            solve_NS_ISO_from_XCTS_config(bconfig, seq);
        bconfig = NS_XCTS_BASE::base_config_t(initial_guess_filename);
        bconfig.control(CONTROLS::SEQUENCES) = false;
    }
    exit_status =
        launch_final_stage_driver(bconfig, seq, resolution, outputdir);

    return exit_status;
}

/** @}*/
};  // namespace Kadath::FUKA_Solvers
