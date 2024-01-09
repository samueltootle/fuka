/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<typename config_t>
int ns_isotropic_norot_stationary_driver (config_t& bconfig, 
  std::string outputdir, ns_sequence const * seq) {
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // make sure NS directory exists for outputs
  if(outputdir == "./") {
    std::filesystem::path cwd = std::filesystem::current_path();
    outputdir = cwd.string();
  }

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
    Scalar lap_Aterm (space, ff1) ;
    Scalar nu        (space, ff1) ;
    Scalar logh      (space, ff1) ;
    Scalar lap_Bterm (space, ff1) ;
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
        ns_solver(bconfig, space, nu, lap_Aterm, logh, lap_Bterm);
      exit_status = ns_solver.solve(seq);

    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS) == 0) ? 
        2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      ns_isotropic_norot_solver<eos_t, decltype(bconfig), decltype(space)> 
        ns_solver(bconfig, space, nu, lap_Aterm, logh, lap_Bterm);
      
      exit_status = ns_solver.solve(seq);
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return exit_status;
}

template<class config_t, class Res_t>
inline int ns_isotropic_norot_driver (config_t& bconfig, 
  Res_t& resolution, std::string outputdir, ns_sequence const * seq) {
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

  bool res_inc = (resolution.final() > resolution.init() || bconfig.control(CONTROLS::REGRID));
  auto resolution_indices = resolution.get_indices();
  auto const & final_res = resolution.final();
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  exit_status = ns_isotropic_norot_stationary_driver(bconfig, outputdir, seq);
  // We now have a "low" resolution solution for the NS of interest
  // Set this to false to avoid iterative M and CHI
  bconfig.control(CONTROLS::SEQUENCES) = false;

  auto regrid = [&]() {
    std::string fname{"ns_regrid"};

    if(rank == 0)
      exit_status = ns_isotropic_norot_regrid(bconfig, fname);
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.set_filename(fname);
    bconfig.open_config();
    
    stage_enabled[STAGES::NOROT_BC] = true;
  };
  // Since the 2D code focuses on sequences, we always regrid to make sure
  // we start/end on an optimal grid structure.
  regrid();
  exit_status = ns_isotropic_norot_stationary_driver(bconfig, outputdir, seq);

  while(res_inc) {        
    int next_res = bco_utils::next_resolution(bconfig(BCO_PARAMS::BCO_RES));
    // iterative res increase
    if(bconfig.control(CONTROLS::REGRID)) {
      bconfig.control(CONTROLS::REGRID) = false;
      res_inc = (resolution.final() > resolution.init());
    } else if(next_res >= final_res) {
      bconfig.set(BCO_PARAMS::BCO_RES) = next_res;
      res_inc = false;
    } else {
      bconfig.set(BCO_PARAMS::BCO_RES) = next_res;
    }
    regrid();

    exit_status = ns_isotropic_norot_stationary_driver(bconfig, outputdir, seq);
  }
  return exit_status;
}
/** @}*/
}}
