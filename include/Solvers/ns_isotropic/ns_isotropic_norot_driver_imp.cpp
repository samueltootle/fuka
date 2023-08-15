/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<class config_t>
config_t ns_isotropic_norot_sequence_setup (config_t & seqconfig, std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig = generate_sequence_config(seqconfig, outputdir);

  update_eos_parameters(seqconfig, bconfig);

  if(rank == 0) bconfig.write_config();
  return bconfig;
}

template<class Seq_t, class Res_t, class config_t>
config_t ns_isotropic_norot_sequence (config_t & seqconfig, 
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
  config_t base_config = ns_isotropic_norot_sequence_setup(seqconfig, outputdir);
  base_config.set(resolution_indices) = resolution.init();

  // Save this in case an invalid ADM mass is given for the EOS used
  const double final_MADM = base_config(BCO_PARAMS::MADM);
  config_t bconfig{};
  #ifdef DEBUG
  std::cout << seq << std::endl;
  std::cout << resolution << std::endl;
  #endif
  auto single_seq = [&](auto val) {
    bconfig = base_config;
    if(seq.is_set())
      bconfig.set(sequence_var_indices) = val;

    if(bconfig.control(CONTROLS::SEQUENCES)) {
      if(rank == 0) {
        setup_2dns_isotropic(bconfig);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      // make sure all ranks have the same config
      bconfig.open_config();
      MPI_Barrier(MPI_COMM_WORLD);
      bconfig.control(CONTROLS::ITERATIVE_M) = 
        (bconfig(BCO_PARAMS::MADM) < final_MADM);

      if(bconfig.control(CONTROLS::ITERATIVE_M)) {
        if(rank == 0)
        std::cerr << "Cannot solve TOV for Madm = " << final_MADM
                  << " without spin.\n";
        std::_Exit(EXIT_FAILURE);
      }
    }

    exit_status = ns_isotropic_norot_driver(bconfig, resolution, outputdir); 
    return exit_status;
  };

  // Loop if a valid sequence is set
  if(seq.is_set()) {
    for(auto val = seq.init(); seq.loop_condition(val); val+=dx) {
      exit_status = single_seq(val);
    }
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
int ns_isotropic_norot_stationary_driver (config_t& bconfig, std::string outputdir){
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // make sure NS directory exists for outputs
  if(outputdir == "./") {
    std::filesystem::path cwd = std::filesystem::current_path();
    outputdir = cwd.string();
  }
  if(rank == 0)
    std::cout << "Solutions will be stored in: " << outputdir << "\n" \
              << "Directory will be created if it doesn't exist.\n";
  fs::create_directory(outputdir);

  // Not important atm
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

  while(exit_status == RELOAD_FILE) { 
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
    Space_polar_adapted space (ff1) ;

    // load the fields defined on the space
    Scalar lap_Aterm   (space, ff1) ;
    Scalar nu  (space, ff1) ;
    Scalar logh   (space, ff1) ;
    fclose(ff1) ;
    
    if(outputdir != "") bconfig.set_outputdir(outputdir) ;

    // load and setup the EOS
    const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT);
    const std::string eos_file = bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE);
    const std::string eos_type = bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t, eos_var_t::PRESSURE>::init(eos_file, h_cut);
      ns_isotropic_norot_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, nu, lap_Aterm, logh);
      exit_status = ns_solver.solve();

    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS) == 0) ? 
        2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      ns_isotropic_norot_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, nu, lap_Aterm, logh);
      
      exit_status = ns_solver.solve();
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return exit_status;
}

template<class config_t, class Res_t>
inline int ns_isotropic_norot_driver (config_t& bconfig, Res_t& resolution, std::string outputdir) {
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
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  exit_status = ns_isotropic_norot_stationary_driver(bconfig, outputdir);
  // We now have a "low" resolution solution for the NS of interest
  // Set this to false to avoid iterative M and CHI
  bconfig.control(CONTROLS::SEQUENCES) = false;

  // auto regrid = [&]() {
  //   std::string fname{"ns_regrid"};

  //   if(rank == 0)
  //     exit_status = ns_isotropic_norot_regrid(bconfig, fname);
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   bconfig.set_filename(fname);
  //   bconfig.open_config();
    
  //   stage_enabled.fill(false);
  //   stage_enabled[last_stage_idx] = true;
  // };

  // while(res_inc) {        

  //   // iterative res increase
  //   if(bconfig(BCO_PARAMS::BCO_RES) + 2 >= final_res) {
  //     bconfig.set(BCO_PARAMS::BCO_RES) = final_res;
  //     res_inc = false;
  //   } else {
  //     bconfig.set(BCO_PARAMS::BCO_RES) += 2;
  //   }
  //   regrid();
  
  //   // Placeholder
  //   // Rerun with new grid
  //   // if(stage_enabled[STAGES::LINBOOST]) {
  //   //   exit_status = bh_3d_xcts_linear_boost_driver(bconfig, outputdir);
  //   // } else {
  //   //   exit_status = bh_3d_xcts_stationary_driver(bconfig, outputdir);
  //   // }
  //   exit_status = ns_isotropic_norot_stationary_driver(bconfig, outputdir);
  // }
  return exit_status;
}
/** @}*/
}}
