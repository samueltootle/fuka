#include "bns_xcts_driver.hpp"
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

/**
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template<class config_t>
int bns_xcts_solution_driver (config_t& bconfig, std::string outputdir) {
  int exit_status = 0;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // make sure outputdir directory exists for outputs
  if(rank == 0)
    std::cout << "Solutions will be stored in: " << outputdir << "\n" \
              << "Directory will be created if it doesn't exist.\n";

  std::string spacein = bconfig.space_filename();
  if(!fs::exists(spacein)) {
    // For debugging MPI bugs
    if(rank == 0) {
      std::cerr << "File: " << spacein << " not found.\n\n";
    } else {
      std::cerr << "File: " << spacein << " not found for another rank.\n\n";
    }
    std::_Exit(EXIT_FAILURE);
  }

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
  Space_bin_ns space  (ff1) ;
	Scalar conf  (space, ff1) ;
	Scalar lapse (space, ff1) ;
  Vector shift (space, ff1) ;
  Scalar logh  (space, ff1) ;
  Scalar phi   (space, ff1) ;
	fclose(ff1) ;
  Base_tensor basis(space, CARTESIAN_BASIS);
  
  if(outputdir != "") {
    fs::create_directory(outputdir);
    bconfig.set_outputdir(outputdir) ;
  }
  if(bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  // load and setup the EOS
  const double h_cut = 
    bconfig.template eos<double>(EOS_PARAMS::HCUT, NODES::BCO1);
  const std::string eos_file = 
    bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE, NODES::BCO1);
  const std::string eos_type = 
    bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE, NODES::BCO1);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = ::Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,eos_var_t::PRESSURE>::init(eos_file, h_cut);
    bns_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
      bns_solver(bconfig, space, basis, conf, lapse, shift, logh, phi);
    bns_solver.solve();
  } else if(eos_type == "Cold_Table") {
    using eos_t = ::Kadath::Margherita::Cold_Table;

    const int interp_pts = 
      (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS, NODES::BCO1) == 0) ? \
        2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS, NODES::BCO1);

    EOS<eos_t,eos_var_t::PRESSURE>::init(eos_file, h_cut, interp_pts);
    bns_xcts_solver<eos_t, decltype(bconfig), decltype(space)> 
      bns_solver(bconfig, space, basis, conf, lapse, shift, logh, phi);
    bns_solver.solve();
  } else { 
    std::cerr << "Unknown EOSTYPE." << endl;
    std::_Exit(EXIT_FAILURE);
  }
  
  return exit_status;
}

template<class Res_t, class config_t>
int bns_xcts_driver (config_t & bconfig, 
                          Res_t const & resolution,
                          std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<bool, NUM_STAGES> const final_stages = bconfig.return_stages();
  auto & stages = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, final_stages);
  
  // solve at low res first    
  int const final_res = resolution.final();
  int const init_res = resolution.init();
  auto resolution_indices = resolution.get_indices();
    
  bool res_inc = (final_res > init_res);
  if(res_inc) {
    // Solve only TOTAL - hydrostatic equilibrium -
    // before increasing resolution
    // All other stages only rescale the matter
    stages.fill(false);
    stages[TOTAL] = true;
  }

  auto regrid = [&]() {
    std::string fname{"bns_regrid"};

    if(rank == 0)
      bns_xcts_regrid(bconfig, fname);
    bconfig.set_filename(fname);
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.open_config();

    // The resolution is iteratively increased
    // for the hydrostatic equilibrium solution.
    // The solutions where hydrostatic equilibrium
    // are not rigorously enforced, i.e. the matter is 
    // rescaled instead of solved for consistently, can only
    // be solved for at the end
    if(!res_inc) {
      // Ensure input stages are activated for the final solution
      // e.g. TOTAL_BC, ECC_RED
      stages = final_stages;
    }
  };

  int err = bns_xcts_solution_driver(bconfig, outputdir);
  while(res_inc) {

    // iterative res increase
    if(bconfig(resolution_indices) + 2 >= final_res) {
      bconfig.set(resolution_indices) = final_res;
      res_inc = false;
    } else {
      bconfig.set(resolution_indices) += 2;
    }
    regrid();
    err = bns_xcts_solution_driver(bconfig, outputdir);
  }
  return err;
}

template<class config_t>
config_t bns_xcts_sequence_setup (config_t & seqconfig, std::string outputdir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig = binary_generate_sequence_config(seqconfig, {"ns","ns"}, outputdir);

  update_eos_parameters(seqconfig, bconfig, NODES::BCO1);
  update_eos_parameters(seqconfig, bconfig, NODES::BCO2);

  if(rank == 0) bconfig.write_config();
  return bconfig;
}

template<class Seq_t, class Res_t, class config_t>
int bns_xcts_sequence (config_t & seqconfig, 
                          Seq_t const & seq,
                          Res_t const & resolution,
                          std::string outputdir) {

  int rank = 0, exit_status = EXIT_SUCCESS;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Initialize sequence variables
  auto sequence_var_indices = seq.get_indices();
  auto resolution_indices   = resolution.get_indices();

  auto const & dx = seq.step_size();

  // Initialize full configurator
  config_t base_config = bns_xcts_sequence_setup(seqconfig, outputdir);
  // Needed to obtain hydrostatic equilibrium solution
  base_config.set_stage(STAGES::TOTAL) = true;

  base_config.set(resolution_indices) = resolution.init();

  #ifdef DEBUG
  if(rank == 0) {
    std::cout << seq << std::endl;
    std::cout << resolution << std::endl;
  }
  #endif

  auto single_seq = [&](auto val) {
    config_t bconfig = base_config;
    if(seq.is_set())
      bconfig.set(sequence_var_indices) = val;

    if(bconfig.control(CONTROLS::SEQUENCES)) {
      // Retain adot in case it is set from the start manually
      auto const adot = bconfig.set(BIN_PARAMS::ADOT);

      // Adot is deleted here
      bns_xcts_setup_bin_config(bconfig);

      // Reset ADOT
      bconfig.set(BIN_PARAMS::ADOT) = adot;

      bns_xcts_setup_space(bconfig);
    }
    
    // Execute ID solution driver
    exit_status = bns_xcts_driver(bconfig, resolution, outputdir);
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
  return exit_status;
}
/** @}*/
}}
