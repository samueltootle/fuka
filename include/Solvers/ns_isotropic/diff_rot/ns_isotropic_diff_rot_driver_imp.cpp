/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template <typename config_t>
void initialize_fields_diff(config_t& bconfig) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string spacein = bconfig.space_filename();
  FILE* ff1 = fopen(spacein.c_str(), "r");
  if (ff1 == NULL) {
    // mainly for debugging MPI bugs
    std::cerr << spacein.c_str() << " failed to open for rank " << rank << "\n";
    std::_Exit(EXIT_FAILURE);
  }
  Space_polar_adapted space(ff1);

  // load the fields from non-rotating solution
  Scalar lap_Aterm(space, ff1);
  Scalar nu(space, ff1);
  Scalar logh(space, ff1);
  Scalar lap_Bterm(space, ff1);
  Scalar lap_wterm(space, ff1);
  fclose(ff1);

  // Construct initial guess for Omega field
  Scalar Omega(space);
  Omega = (std::fabs(bconfig(BCO_PARAMS::OMEGA)) < 1e-7)
              ? 1e-7
              : bconfig(BCO_PARAMS::OMEGA);
  Omega.std_base();

  bconfig.set_filename("initns_diffrot");
  bconfig.set_field(BCO_FIELDS::DIFF_OMEGA) = true;
  bconfig.control(CONTROLS::SEQUENCES) = true;
  if (rank == 0)
    bco_utils::save_to_file(space, bconfig, lap_Aterm, nu, logh, lap_Bterm,
                            lap_wterm, Omega);
  MPI_Barrier(MPI_COMM_WORLD);
}

template <class eos_t>
class launch_diffrot_solver {
  int exit_status{RELOAD_FILE};

  template <class config_t>
  void regrid(const int rank, config_t& bconfig, bool& regridded) {
    std::string fname{"ns_regrid"};
    std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();

    if (rank == 0)
      exit_status = ns_isotropic_diff_rot_regrid(bconfig, fname);
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.set_filename(fname);
    bconfig.open_config();

    stage_enabled.fill(false);
    stage_enabled[STAGES::DIFF_ROT] = true;
    exit_status = RELOAD_FILE;
    regridded = true;
  }

 public:
  template <class config_t>
  int operator()(const int rank,
                 config_t& bconfig,
                 std::string outputdir,
                 bool& regridded,
                 ns_sequence const* seq) {

    std::string spacein = bconfig.space_filename();

    if (!fs::exists(spacein)) {
      // mainly for debugging MPI bugs
      if (rank == 0) {
        std::cerr << "File: " << spacein << " not found.\n\n";
      } else {
        std::cerr << "File: " << spacein << " not found for another rank.\n\n";
      }
      std::_Exit(EXIT_FAILURE);
    }

    // just so you really know
    if (rank == 0) {
      std::cout << "Config File: " << bconfig.config_filename_abs() << std::endl
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
    Space_polar_adapted space(ff1);

    // load the fields defined on the space
    Scalar lap_Aterm(space, ff1);
    Scalar nu(space, ff1);
    Scalar logh(space, ff1);
    Scalar lap_Bterm(space, ff1);
    Scalar lap_wterm(space, ff1);
    Scalar Omega(space, ff1);
    fclose(ff1);

    if (outputdir != "") {
      bconfig.set_outputdir(outputdir);
    }

    ns_isotropic_diff_rot_solver<eos_t, decltype(bconfig), decltype(space)>
        ns_solver(bconfig, space, nu, lap_Aterm, logh, lap_Bterm, lap_wterm,
                  Omega);
    exit_status = ns_solver.solve(seq);

    if (exit_status == EXIT_SUCCESS && !regridded &&
        bconfig.control(CONTROLS::ITERATIVE_RRATIO)) {
      regrid(rank, bconfig, regridded);
    } else if (bconfig.control(CONTROLS::ITERATIVE_RRATIO)) {
      double rr = bconfig.template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO);
      rr -= 0.05;
      bconfig.set_diffrot(DIFFROT_PARAMS::DIFF_RRATIO) =
          (rr < bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO))
              ? bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO)
              : rr;
      exit_status = RELOAD_FILE;
      bconfig.control(CONTROLS::ITERATIVE_RRATIO) =
          (bconfig.template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO) !=
           bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO));
      regridded = false;
    }
    return exit_status;
  };
};

/**
 * @brief Driver to compute a stationary solution for a given resolution
 *
 * @tparam config_t Config file type
 * @param bconfig NS config file
 * @param outputdir directory to store solutions in
 * @return int error code
 */
template <typename config_t>
int ns_isotropic_diff_rot_stationary_driver(config_t& bconfig,
                                            std::string outputdir,
                                            ns_sequence const* seq) {
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  using namespace Kadath::FUKA_EOS;

  // make sure NS directory exists for outputs
  if (outputdir == "./") {
    fs::path cwd = fs::current_path();
    outputdir = cwd.string();
  }

  // Make sure fields needed for rotating solution are initialized before
  // opening files
  if (!bconfig.field(BCO_FIELDS::DIFF_OMEGA))
    initialize_fields_diff(bconfig);

  bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO) =
      bconfig.template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO);

  bconfig.control(CONTROLS::ITERATIVE_RRATIO) =
      bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO) < 0.9 &&
      bconfig.control(CONTROLS::SEQUENCES);

  bconfig.set_diffrot(DIFFROT_PARAMS::DIFF_RRATIO) =
      (bconfig.control(CONTROLS::ITERATIVE_RRATIO))
          ? 0.9
          : bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO);
  // cout << "RATIOS: " << (0 == 0) << " = " <<
  // (bconfig.seq_setting(SEQ_SETTINGS::FINAL_RRATIO) >= 0.9) << " - " <<
  // !bconfig.control(CONTROLS::SEQUENCES) << ", " << bconfig.template
  // diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO) << '\n'; Not important atm
  // if(std::isnan(bconfig.set(BCO_PARAMS::MADM)) &&
  // std::isnan(bconfig.set(BCO_PARAMS::MB))){
  //   if(rank == 0)
  //     std::cout << "Config error.  No madm nor mb found. \n\n";
  //   std::_Exit(EXIT_FAILURE);
  // }

  const std::string eos_type =
      bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE);
  bool regridded = false;
  while (exit_status == RELOAD_FILE) {
    exit_status = EOS_Function_Dispatcher::dispatch<launch_diffrot_solver>(
        bconfig, eos_type, rank, bconfig, outputdir, regridded, seq);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  bconfig.control(CONTROLS::SEQUENCES) = false;
  return exit_status;
}

template <class config_t, class Res_t>
inline int ns_isotropic_diff_rot_driver(config_t& bconfig,
                                        Res_t& resolution,
                                        std::string outputdir,
                                        ns_sequence const* seq) {
  int exit_status = RELOAD_FILE;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bool res_inc = (resolution.final() > resolution.init());
  auto resolution_indices = resolution.get_indices();
  auto const& final_res = resolution.final();

  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

  exit_status =
      ns_isotropic_diff_rot_stationary_driver(bconfig, outputdir, seq);

  auto regrid = [&]() {
    std::string fname{"ns_regrid"};

    if (rank == 0)
      exit_status = ns_isotropic_diff_rot_regrid(bconfig, fname);
    MPI_Barrier(MPI_COMM_WORLD);
    bconfig.set_filename(fname);
    bconfig.open_config();

    stage_enabled.fill(false);
    stage_enabled[STAGES::DIFF_ROT] = true;
  };
  // Since the 2D code focuses on sequences, we always regrid to make sure
  // we start/end on an optimal grid structure.
  regrid();
  exit_status =
      ns_isotropic_diff_rot_stationary_driver(bconfig, outputdir, seq);

  while (res_inc) {
    int next_res = bco_utils::next_resolution(bconfig(resolution_indices));
    // iterative res increase
    if (next_res >= final_res) {
      bconfig.set(resolution_indices) = next_res;
      res_inc = false;
    } else {
      bconfig.set(resolution_indices) = next_res;
    }
    regrid();

    exit_status =
        ns_isotropic_diff_rot_stationary_driver(bconfig, outputdir, seq);
  }
  return exit_status;
}

/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath
