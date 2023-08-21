/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<class config_t>
config_t ns_3d_xcts_sequence_setup (config_t & seqconfig, std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig = generate_sequence_config(seqconfig, outputdir);

  update_eos_parameters(seqconfig, bconfig);

  if(rank == 0) bconfig.write_config();
  return bconfig;
}

template<class Seq_t, class Res_t, class config_t>
config_t ns_3d_xcts_sequence (config_t & seqconfig, 
                          Seq_t const & seq,
                          Res_t const & resolution,
                          std::string outputdir) {
  
  int rank = 0, exit_status = EXIT_SUCCESS;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Initialize sequence variables
  auto sequence_var_indices = seq.get_indices();
  auto resolution_indices   = resolution.get_indices();
  
  auto const & dx = seq.step_size();

  // // Initialize full configurator
  config_t base_config = ns_3d_xcts_sequence_setup(seqconfig, outputdir);
  base_config.set(resolution_indices) = resolution.init();

  // in the event the user wants an MADM > MTOV
  const double final_MADM = base_config(BCO_PARAMS::MADM);
  base_config.control(CONTROLS::ITERATIVE_M) = false;
  config_t bconfig{base_config};

  if(bconfig.control(CONTROLS::SEQUENCES) || bconfig.control(CONTROLS::RESOLVE)) {
    if(rank == 0) {
      setup_co<NODES::NS>(bconfig);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // make sure all ranks have the same config
    bconfig.open_config();
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.control(CONTROLS::ITERATIVE_M) = 
      (bconfig(BCO_PARAMS::MADM) < final_MADM);

    if(bconfig.control(CONTROLS::ITERATIVE_M) 
        && std::fabs(bconfig(BCO_PARAMS::CHI)) < 1e-5) {
      if(rank == 0)
      std::cerr << "Cannot solve TOV for Madm = " << final_MADM
                << " without spin.\n";
      std::_Exit(EXIT_FAILURE);
    }
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
    exit_status = ns_3d_xcts_base_solution_driver(bconfig, outputdir);

    // Update config such that the next solving round uses
    // the final ADM mass and spin
    bconfig(BCO_PARAMS::MADM) = final_MADM;
    bconfig.control(CONTROLS::SEQUENCES) = false;
    

  }
  exit_status = ns_3d_xcts_base_solution_driver(bconfig, outputdir);
  // Ensure only the final stage is used
  // e.g. avoid NOROT stage
  stage_enabled.fill(false);
  stage_enabled[last_stage_idx] = true;

  base_config = bconfig;
  if(seq.is_set() && std::isnan(bconfig.set(sequence_var_indices)))
    bconfig.set(sequence_var_indices) = seq.init();
  #ifdef DEBUG
  std::cout << seq << std::endl;
  std::cout << resolution << std::endl;
  #endif
  auto single_seq = [&](auto val) {
    // bconfig = base_config;
    stage_enabled[last_stage_idx] = true;
    auto old_val = bconfig(sequence_var_indices);
    if(seq.is_set())
      bconfig.set(sequence_var_indices) = val;
    if(std::get<0>(sequence_var_indices) == BCO_PARAMS::CHI) {
      bconfig.seq_setting(SEQ_SETTINGS::FINAL_CHI) = bconfig(BCO_PARAMS::CHI);
      if(std::fabs(old_val) < 0.1 && val > 0.2) {
        bconfig.control(CONTROLS::SEQUENCES) = true;
        exit_status = ns_3d_xcts_base_solution_driver(bconfig, outputdir);
      }
    }
    exit_status = ns_3d_xcts_driver(bconfig, resolution, outputdir, &seq); 
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
int ns_3d_xcts_stationary_driver (config_t& bconfig, std::string outputdir, Parameter_sequence<BCO_PARAMS> const * seq){
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // make sure NS directory exists for outputs
  if(outputdir == "./") {
    std::filesystem::path cwd = std::filesystem::current_path();
    outputdir = cwd.string();
  }
  
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

  spacein = bconfig.space_filename();
  // just so you really know
  if(rank == 0) {
    std::cout << "Config File: " 
              << bconfig.config_filename_abs() << std::endl
              << "Fields File: " << spacein << std::endl
              << bconfig << std::endl;
  }
  FILE* ff1 = fopen (spacein.c_str(), "r") ;
  if(ff1 == NULL){
    // mainly for debugging MPI bugs
    std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
    std::_Exit(EXIT_FAILURE);
  }
  Space_spheric_adapted space (ff1) ;
  Scalar conf   (space, ff1) ;
  Scalar lapse  (space, ff1) ;
  Vector shift  (space, ff1) ;
  Scalar logh   (space, ff1) ;
  fclose(ff1) ;
  Base_tensor basis(space, CARTESIAN_BASIS);
  
  if(outputdir != "") bconfig.set_outputdir(outputdir) ;

  // load and setup the EOS
  const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT);
  const std::string eos_file = bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE);
  const std::string eos_type = bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t, eos_var_t::PRESSURE>::init(eos_file, h_cut);
    ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
      ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
    exit_status = ns_solver.solve(seq);

  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS) == 0) ? 
      2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
      ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
    
    exit_status = ns_solver.solve(seq);
  } else { 
    std::cerr << "Unknown EOSTYPE." << endl;
    std::_Exit(EXIT_FAILURE);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  return exit_status;
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
int ns_3d_xcts_base_solution_driver (config_t& bconfig, std::string outputdir){
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

  if(std::isnan(bconfig.set(BCO_PARAMS::MADM)) && std::isnan(bconfig.set(BCO_PARAMS::MB))){
    if(rank == 0)
      std::cout << "Config error.  No madm nor mb found. \n\n";
    std::_Exit(EXIT_FAILURE);
  }
  
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

  while(exit_status == RELOAD_FILE) { 
    exit_status = ns_3d_xcts_stationary_driver(bconfig, outputdir);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return exit_status;
}

template<class config_t, class Res_t>
inline int ns_3d_xcts_driver (config_t& bconfig, 
  Res_t& resolution, std::string outputdir, Parameter_sequence<BCO_PARAMS> const * seq) {
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

  double const xboost = bconfig.set(BCO_PARAMS::BVELX);
  double const yboost = bconfig.set(BCO_PARAMS::BVELY);
  bool res_inc = (resolution.final() > resolution.init());
  auto resolution_indices = resolution.get_indices();
  auto const & final_res = resolution.final();
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  // Placeholder
  // These need to be zero prior to computing
  // a stationary solution
  bconfig.set(BCO_PARAMS::BVELX) = 0.;
  bconfig.set(BCO_PARAMS::BVELY) = 0.;

  exit_status = ns_3d_xcts_stationary_driver(bconfig, outputdir, seq);
  // We now have a "low" resolution solution for the NS of interest
  // Set this to false to avoid iterative M and CHI
  bconfig.control(CONTROLS::SEQUENCES) = false;

  // Placeholder
  // if(stage_enabled[STAGES::LINBOOST]) {
  //   bconfig.set(BCO_PARAMS::BVELX) = xboost;
  //   bconfig.set(BCO_PARAMS::BVELY) = yboost;
  //   exit_status = bh_3d_xcts_linear_boost_driver(bconfig, outputdir);
  // }

  auto regrid = [&]() {
    std::string fname{"ns_regrid"};

    if(rank == 0)
      exit_status = ns_3d_xcts_regrid(bconfig, fname);
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.set_filename(fname);
    bconfig.open_config();
    
    stage_enabled.fill(false);
    stage_enabled[last_stage_idx] = true;
  };

  while(res_inc) {        

    // iterative res increase
    if(bconfig(BCO_PARAMS::BCO_RES) + 2 >= final_res) {
      bconfig.set(BCO_PARAMS::BCO_RES) = final_res;
      res_inc = false;
    } else {
      bconfig.set(BCO_PARAMS::BCO_RES) += 2;
    }
    regrid();
  
    // Placeholder
    // Rerun with new grid
    // if(stage_enabled[STAGES::LINBOOST]) {
    //   exit_status = bh_3d_xcts_linear_boost_driver(bconfig, outputdir);
    // } else {
    //   exit_status = bh_3d_xcts_stationary_driver(bconfig, outputdir);
    // }
    exit_status = ns_3d_xcts_stationary_driver(bconfig, outputdir, seq);
  }
  return exit_status;
}

template<typename config_t, class Res_t>
inline int ns_3d_xcts_binary_boost_driver (config_t& bconfig, 
  Res_t& resolution, std::string outputdir,
    kadath_config_boost<BIN_INFO> binconfig, const size_t bco) {

  int exit_status = RUN_BOOST;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int final_res = bconfig(BCO_PARAMS::BCO_RES);
  bool res_inc = (bconfig.seq_setting(SEQ_SETTINGS::INIT_RES) < final_res);
  bconfig.set(BCO_PARAMS::BCO_RES) = resolution.init();

  // Obtain stationary solution
  Parameter_sequence<BCO_PARAMS> tmp_seq{};
  bconfig = ns_3d_xcts_sequence(bconfig, tmp_seq, resolution, outputdir);
  
  while(exit_status == RUN_BOOST) { 
    auto spacein = bconfig.space_filename();
    // just so you really know
    if(rank == 0) {
      std::cout << "Config File: " 
                << bconfig.config_filename_abs() << std::endl
                << "Fields File: " << spacein << std::endl
                << bconfig << std::endl;
    }
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    if(ff1 == NULL){
      // mainly for debugging MPI bugs
      std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
      std::_Exit(EXIT_FAILURE);
    }
    Space_spheric_adapted space (ff1) ;
    Scalar conf   (space, ff1) ;
    Scalar lapse  (space, ff1) ;
    Vector shift  (space, ff1) ;
    Scalar logh   (space, ff1) ;
    fclose(ff1) ;
    Base_tensor basis(space, CARTESIAN_BASIS);

    if(outputdir != "") bconfig.set_outputdir(outputdir) ;
        // load and setup the EOS
    const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT);
    const std::string eos_file = bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE);
    const std::string eos_type = bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
      
      exit_status = ns_solver.binary_boost_stage(binconfig, bco);
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS) == 0) ? \
                              2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      ns_3d_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, basis, conf, lapse, logh, shift);
      
      exit_status = ns_solver.binary_boost_stage(binconfig, bco);
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return exit_status;
}
/** @}*/
}}
