#include "NS_XCTS.hpp"
#include "Solvers/co_solver_utils.hpp"
#include "name_tools.hpp"
/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath::FUKA_Solvers {

template<class eos_t>
inline int ns_isotropic_norot_driver (NS_XCTS_BASE::base_config_t& bconfig, ns_sequence const & seq, 
  Parameter_sequence<BCO_PARAMS> & resolution, std::string outputdir);

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
int ns_xcts_norot_seq_driver(NS_XCTS_BASE::base_config_t& seqconfig, ns_sequence const & seq, 
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

    ns_isotropic_norot_driver<eos_t>(bconfig, seq, tmp_res, outputdir);
    // solver.solve();
  } else {
    ns_isotropic_norot_driver<eos_t>(bconfig, seq, resolution, outputdir);
    // sequence...
    // solver.solve();
  }
  
  // Update config such that the next solving round uses
  // the final ADM mass and spin if applicable
  bconfig(BCO_PARAMS::MADM) = final_MADM;
  bconfig.control(CONTROLS::SEQUENCES) = false;

  return EXIT_SUCCESS;
}

template<class eos_t>
inline int ns_isotropic_norot_driver (NS_XCTS_BASE::base_config_t& bconfig, ns_sequence const & seq, 
  Parameter_sequence<BCO_PARAMS> & resolution, std::string outputdir) {
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // make sure NS directory exists for outputs
  if(outputdir == "./") {
    std::filesystem::path cwd = std::filesystem::current_path();
    outputdir = cwd.string();
  }

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
    std::stringstream ss;
    ss << spacein.c_str() << " failed to open for rank " << rank << "\n";
    throw std::runtime_error(ss.str().c_str());
  }

  NS_XCTS_NOROT<eos_t> solver(bconfig, seq, resolution, outputdir);
  solver.setup_syst();
  solver.do_newton();
  
  MPI_Barrier(MPI_COMM_WORLD);
  return exit_status;
}

template<class eos_t>
inline int launch_final_stage_driver(NS_XCTS_BASE::base_config_t& bconfig, ns_sequence const & seq, 
  Parameter_sequence<BCO_PARAMS> & resolution, std::string outputdir) {
  
  using config_t = NS_XCTS_BASE::base_config_t;
  using Res_t = Parameter_sequence<BCO_PARAMS>;

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  std::function<int(config_t&, ns_sequence const &, Res_t&, std::string)> final_stage_driver;
  if(rank == 0)
    std::cout << "Last stage: " << last_stage << '\n';
  if(seq.is_set() || bconfig.control(CONTROLS::SEQUENCES)) {
    final_stage_driver = &ns_xcts_norot_seq_driver<eos_t>;
  } else {
    switch(last_stage_idx) {
      case STAGES::NOROT_BC:
        final_stage_driver = &ns_isotropic_norot_driver<eos_t>;
        break;
      // case STAGES::UNIFORM_ROT:
      //   final_stage_driver = &ns_isotropic_uniform_rot_driver<config_t, Res_t>;
      //   break;
      // case STAGES::DIFF_ROT:
      //   final_stage_driver = &ns_isotropic_diff_rot_driver<config_t, Res_t>;
      //   break;
      default:
        throw std::runtime_error("No valid stages enabled.\n");
    }
  }
  return final_stage_driver(bconfig, seq, resolution, outputdir);
}

inline int ns_3d_xcts_driver (NS_XCTS_BASE::base_config_t& bconfig, ns_sequence const & seq, 
  Parameter_sequence<BCO_PARAMS> & resolution, std::string outputdir) {
  
  int exit_status = EXIT_SUCCESS;
  auto resolution_indices = resolution.get_indices();
  bconfig.set(resolution_indices) = resolution.init();  

  // load and setup the EOS
  const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT);
  const std::string eos_file = bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE);
  const std::string eos_type_ = bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);

  const std::string eos_type = str_tolower(eos_type_);
  if(eos_type == "cold_pwpoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    exit_status = launch_final_stage_driver<eos_t>(bconfig, seq, resolution, outputdir);
  } else if(eos_type == "cold_table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS) == 0) ? \
                            2000 : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    exit_status = launch_final_stage_driver<eos_t>(bconfig, seq, resolution, outputdir);
  } else { 
    throw std::runtime_error("Unknown EOS type \n");
  }
  return exit_status;
}
/** @}*/
};