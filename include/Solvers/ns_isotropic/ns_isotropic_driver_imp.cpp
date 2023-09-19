#include<functional>
#include"Solvers/sequences/ns_sequence.hpp"
/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<class config_t>
config_t ns_isotropic_sequence_setup (config_t & seqconfig, std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig = generate_sequence_config(seqconfig, outputdir);

  update_eos_parameters(seqconfig, bconfig);
  update_diffrot_parameters(seqconfig, bconfig);

  if(rank == 0) bconfig.write_config();
  return bconfig;
}

template<class Res_t, class config_t>
config_t ns_isotropic_sequence (config_t & seqconfig, 
                          ns_sequence const & seq,
                          Res_t & resolution,
                          std::string outputdir) {
  
  int rank = 0, exit_status = EXIT_SUCCESS;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // Ensure fixed values are initialized
  seqconfig.set(seq.mass_idx()) = seq.mass_val();
  seqconfig.set(seq.spin_idx()) = seq.spin_val();

  // Initialize sequence variables
  auto sequence_idx = seq.get_sequence_idx();
  auto resolution_indices = resolution.get_indices();
  
  auto const & dx = seq.step_size();

  // Initialize full configurator
  config_t base_config = ns_isotropic_sequence_setup(seqconfig, outputdir);
  base_config.set(resolution_indices) = resolution.init();

  if(!(base_config.control(CONTROLS::SEQUENCES) || base_config.control(CONTROLS::RESOLVE))) {
    base_config = seqconfig;
    base_config.set(BCO_PARAMS::CHI) = 0;
  }
  auto mass_fixing = seq.mass_idx();

  // Should be deprecated...
  if(seq.is_set() && std::isnan(base_config.set(sequence_idx))) {
    base_config.set(sequence_idx) = seq.init();
    if(ns_seq_is_mass_fixing(seq))
      mass_fixing = sequence_idx;
  }

  // Save this in case an invalid ADM mass is given for the EOS used
  const double final_MADM = (seq.mass_idx() == BCO_PARAMS::MADM) ? base_config(BCO_PARAMS::MADM) : std::nan("1");
  base_config.control(CONTROLS::ITERATIVE_M) = false;
  config_t bconfig{base_config};
  
  // Not tested...
  if(bconfig.control(CONTROLS::SEQUENCES) || bconfig.control(CONTROLS::RESOLVE)) {
    int n_shells = (std::isnan(bconfig.set(BCO_PARAMS::NSHELLS))) ? 0 : bconfig.set(BCO_PARAMS::NSHELLS);
    // the 2D code is very sensitive to the initial domain decomposition.  Although polytropic
    // laws are very stable, tabulated solution are more sensitive.
    bconfig.set(BCO_PARAMS::NSHELLS) = 0.;
    if(rank == 0) {
      setup_2dns_isotropic(bconfig, mass_fixing);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // make sure all ranks have the same config
    bconfig.open_config();
    bconfig.set(BCO_PARAMS::NSHELLS) = n_shells;

    bconfig.control(CONTROLS::ITERATIVE_M) = !std::isnan(final_MADM) &&
      (std::fabs(1. - bconfig(BCO_PARAMS::MADM)/final_MADM) > 1e-3);

    if(bconfig.control(CONTROLS::ITERATIVE_M)) {
      if(rank == 0)
      std::cerr << "Cannot solve TOV for Madm = " << final_MADM
                << " without spin. " << bconfig(BCO_PARAMS::MADM) << '\n';
      std::_Exit(EXIT_FAILURE);
    }
    bconfig.control(CONTROLS::SEQUENCES) = true;
  }

  // Need to get a rotating solution before increasing the
  // NS mass up to final_MADM
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);
  if(bconfig.control(CONTROLS::ITERATIVE_M)){    

    // Only obtain the iterative solution at the initial_resolution
    auto const res_init{resolution.init()};
    Parameter_sequence tmp_res("res", BCO_PARAMS::BCO_RES);      
    tmp_res.set(res_init,res_init,res_init);

    // exit_status = ns_3d_xcts_driver(bconfig, tmp_res, outputdir);
    exit_status = ns_isotropic_base_solution_driver(bconfig, outputdir, &seq);

    // Update config such that the next solving round uses
    // the final ADM mass and spin
    bconfig(BCO_PARAMS::MADM) = final_MADM;
    bconfig.control(CONTROLS::SEQUENCES) = false;
  }
  exit_status = ns_isotropic_base_solution_driver(bconfig, outputdir, &seq);
  
  // Ensure only the final stage is used
  // e.g. avoid NOROT stage
  stage_enabled.fill(false);
  stage_enabled[last_stage_idx] = true;

  base_config = bconfig;

  #ifdef DEBUG
  std::cout << seq << std::endl;
  std::cout << resolution << std::endl;
  #endif
  auto single_seq = [&](auto val) {
    stage_enabled[last_stage_idx] = true;
    if(seq.is_set())
      bconfig.set(sequence_idx) = val;

    exit_status = ns_isotropic_driver(bconfig, resolution, outputdir, &seq);
    resolution.set(resolution.final(), resolution.final(), resolution.final());
    return exit_status;
  };

  // Loop if a valid sequence is set
  if(seq.is_set() && std::fabs(dx) > 0) {
    for(auto val = seq.init(); seq.loop_condition(val); val+=dx) {
      exit_status = single_seq(val);
    }
  } else if(seq.is_set()) {
    exit_status = single_seq(seq.init());
  } else {
    exit_status = single_seq(0);
  }
  return bconfig;
}

/**
 * @brief Driver to compute a stationary solution for a given resolution
 * 
 * @tparam config_t Config file type
 * @param bconfig NS config file
 * @param outputdir directory to store solutions in
 * @return int error code
 */
template<typename config_t>
int ns_isotropic_base_solution_driver (config_t& bconfig, std::string outputdir, ns_sequence const * seq=nullptr){
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(std::abs( bconfig(BCO_PARAMS::CHI) ) > 0.7) {
    std::cerr << "Unable to handle chi > 0.7\n";
    return EXIT_FAILURE;
  }
  // We need to store the final desired CHI in case of iterative chi
  bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI) = bconfig(BCO_PARAMS::CHI);

  // In the event we wish to solve for a highly spinning solution
  // we need to do an initial slow rotating solution before going to 
  // faster rotations otherwise the solution will diverge.
  bconfig.control(CONTROLS::ITERATIVE_CHI) = 
    bconfig.control(CONTROLS::SEQUENCES) 
      && std::fabs(bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI)) > 0.2;

  // Lower chi in case of iterative chi
  bconfig(BCO_PARAMS::CHI) = (bconfig.control(CONTROLS::ITERATIVE_CHI)) ? 
    std::copysign(0.1, bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI)) : 
    bconfig(BCO_PARAMS::CHI);

  // make sure NS directory exists for outputs
  if(outputdir == "./") {
    std::filesystem::path cwd = std::filesystem::current_path();
    outputdir = cwd.string();
  }
  if(rank == 0)
    std::cout << "Solutions will be stored in: " << outputdir << "\n" \
              << "Directory will be created if it doesn't exist.\n";
  fs::create_directory(outputdir);

  // if(std::isnan(bconfig.set(BCO_PARAMS::MADM)) && std::isnan(bconfig.set(BCO_PARAMS::MB))){
  //   if(rank == 0)
  //     std::cout << "Config error.  No madm nor mb found. \n\n";
  //   std::_Exit(EXIT_FAILURE);
  // }
  
  std::string spacein = bconfig.space_filename();
  if(!fs::exists(spacein)) {
    // mainly for debugging MPI bugs
    if(rank == 0) {
      std::cerr << "File: " << spacein << " not found.\n\n";
    } else {
      std::cerr << "File: " << spacein << " not found for another rank.\n\n";
    }
    std::_Exit(EXIT_FAILURE);
  }
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);
  while(exit_status == RELOAD_FILE) { 
    exit_status == EXIT_FAILURE;
    // exit_status = ns_isotropic_stationary_driver(bconfig, outputdir);
    if(stage_enabled[STAGES::NOROT_BC]) {
      double const omega = bconfig(BCO_PARAMS::OMEGA);
      double const chi = bconfig(BCO_PARAMS::CHI);
      bconfig(BCO_PARAMS::OMEGA) = 0.;
      bconfig(BCO_PARAMS::CHI) = 0.;

      exit_status = ns_isotropic_norot_stationary_driver(bconfig, outputdir, seq);
      stage_enabled[STAGES::NOROT_BC] = false;
      bconfig(BCO_PARAMS::OMEGA) = omega;
      bconfig(BCO_PARAMS::CHI) = chi;
      if(last_stage_idx != STAGES::NOROT_BC) {
        stage_enabled[STAGES::NOROT_BC] = false;
        exit_status = RELOAD_FILE;
      }
      bconfig.return_stages() = stage_enabled;
    } else if(stage_enabled[STAGES::UNIFORM_ROT]) {
      exit_status = ns_isotropic_uniform_rot_stationary_driver(bconfig, outputdir, seq);
      // exit_status = EXIT_FAILURE;
      if(last_stage_idx != STAGES::UNIFORM_ROT) {
        stage_enabled[STAGES::UNIFORM_ROT] = false;
        exit_status = RELOAD_FILE;
      }
      bconfig.return_stages() = stage_enabled;
    } else if(stage_enabled[STAGES::DIFF_ROT]) {
      exit_status = ns_isotropic_diff_rot_stationary_driver(bconfig, outputdir, seq);
      // exit_status = EXIT_FAILURE;
      if(last_stage_idx != STAGES::DIFF_ROT) {
        stage_enabled[STAGES::DIFF_ROT] = false;
        exit_status = RELOAD_FILE;
      }
      bconfig.return_stages() = stage_enabled;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return exit_status;
}

template<class config_t, class Res_t>
inline int ns_isotropic_driver (config_t& bconfig, Res_t& resolution, 
  std::string outputdir, ns_sequence const * seq) {
  
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string spacein = bconfig.space_filename();
  if(!fs::exists(spacein)) {
    // mainly for debugging MPI bugs
    if(rank == 0) {
      std::cerr << "File: " << spacein << " not found.\n\n";
    } else {
      std::cerr << "File: " << spacein << " not found for another rank.\n\n";
    }
    std::_Exit(EXIT_FAILURE);
  }

  bool res_inc = (resolution.final() > resolution.init());
  auto resolution_indices = resolution.get_indices();
  auto const & final_res = resolution.final();
  bconfig.set(resolution_indices) = resolution.init();
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  std::function<int(config_t&, Res_t&, std::string, ns_sequence const *)> final_stage_driver;
  if(rank == 0)
    std::cout << "Last stage: " << last_stage << '\n';
  switch(last_stage_idx) {
    case STAGES::NOROT_BC:
      final_stage_driver = &ns_isotropic_norot_driver<config_t, Res_t>;
      break;
    case STAGES::UNIFORM_ROT:
      final_stage_driver = &ns_isotropic_uniform_rot_driver<config_t, Res_t>;
      break;
    case STAGES::DIFF_ROT:
      final_stage_driver = &ns_isotropic_diff_rot_driver<config_t, Res_t>;
      break;
  }
  
  exit_status = final_stage_driver(bconfig, resolution, outputdir, seq);
  
  return exit_status;
}
/** @}*/
}}
