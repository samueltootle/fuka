#include "NS_XCTS.hpp"
#include "Solvers/co_solver_utils.hpp"
/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath::FUKA_Solvers {

inline NS_XCTS_BASE::base_config_t ns_xcts_sequence_setup (NS_XCTS_BASE::base_config_t& seqconfig, std::string outputdir) {
  using config_t = NS_XCTS_BASE::base_config_t;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  config_t bconfig = generate_sequence_config(seqconfig, outputdir);

  update_eos_parameters(seqconfig, bconfig);
  update_diffrot_parameters(seqconfig, bconfig);

  if(rank == 0) bconfig.write_config();
  return bconfig;
}

template<class eos_t>
NS_XCTS_BASE::base_config_t ns_xcts_norot_seq_driver(NS_XCTS_BASE::base_config_t& seqconfig, ns_sequence& seq, 
  Parameter_sequence<BCO_PARAMS>& resolution, std::string const outputdir) {
  
  using config_t = NS_XCTS_BASE::base_config_t;
  int rank = 0, exit_status = EXIT_SUCCESS;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Ensure fixed values are initialized
  seqconfig.set(seq.mass_idx()) = seq.mass_val();

  // Initialize sequence variables
  auto sequence_var_indices = seq.get_indices();
  auto resolution_indices   = resolution.get_indices();

  // Initialize full configurator
  config_t base_config = ns_xcts_sequence_setup(seqconfig, outputdir);
  base_config.set(resolution_indices) = resolution.init();
  
  // in the event the user wants an MADM > MTOV
  const double final_MADM = (seq.mass_idx() == BCO_PARAMS::MADM) ? base_config(BCO_PARAMS::MADM) : std::nan("1");
  base_config.control(CONTROLS::ITERATIVE_M) = false;
  config_t bconfig{base_config};

  if(!(base_config.control(CONTROLS::SEQUENCES) || base_config.control(CONTROLS::RESOLVE))) {
    base_config = seqconfig;
    base_config.set(BCO_PARAMS::CHI) = 0;
  }
  
  auto mass_fixing = seq.mass_idx();

  // Get last enabled to make for safety check later
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  if(bconfig.control(CONTROLS::SEQUENCES) || bconfig.control(CONTROLS::RESOLVE)) {
    if(rank == 0) {
      // setup_co<NODES::NS>(bconfig);
      setup_ns_3d_xcts(bconfig, mass_fixing);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // make sure all ranks have the same config
    bconfig.open_config();
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.control(CONTROLS::ITERATIVE_M) = !std::isnan(final_MADM) &&
      (std::fabs(1. - bconfig(BCO_PARAMS::MADM)/final_MADM) > 1e-3);

    if(bconfig.control(CONTROLS::ITERATIVE_M)
      && last_stage_idx == STAGES::NOROT_BC) {
      std::stringstream ss;
      ss << "Cannot solve TOV for Madm = " << final_MADM << " without spin.\n";
      throw std::runtime_error(ss.str().c_str());
    }
  }
  if(last_stage_idx != STAGES::NOROT_BC){

    // Only obtain the iterative solution at the initial_resolution
    auto const res_init{resolution.init()};
    Parameter_sequence tmp_res("res", BCO_PARAMS::BCO_RES);      
    tmp_res.set(res_init,res_init,res_init);

    NS_XCTS_NOROT<eos_t> solver(bconfig, seq, tmp_res, outputdir);
    // solver.solve();
  } else {
    NS_XCTS_NOROT<eos_t> solver(bconfig, seq, resolution, outputdir);
    // sequence...
    // solver.solve();
  }
  
  // Update config such that the next solving round uses
  // the final ADM mass and spin if applicable
  bconfig(BCO_PARAMS::MADM) = final_MADM;
  bconfig.control(CONTROLS::SEQUENCES) = false;

  return bconfig;
}
/** @}*/
};