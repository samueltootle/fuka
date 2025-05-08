/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author:
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "Solvers/bns_xcts/bns_xcts_driver.hpp"
#include "Solvers/bns_xcts/bns_xcts_setup.hpp"
#include "Solvers/bns_xcts/bns_xcts_solver.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/solver_startup.hpp"
#include "mpi.h"
using namespace Kadath::FUKA_EOS;

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int err = EXIT_SUCCESS;

  namespace solvers = ::Kadath::FUKA_Solvers;
  using config_t = kadath_config_boost<BIN_INFO>;
  using InitSolver = solvers::Initialize_Solver<config_t>;

  // Initialize static member variables
  InitSolver::input_configname = "initial_bns.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;

  if (InitSolver::example_setup) {
    if (rank == 0) {
      bconfig.initialize_binary({"ns", "ns"});
      if (InitSolver::minimal_config) {
        bconfig.set_minimal_defaults();
      } else {
        bconfig.set_defaults();
      }

      bconfig.set(BIN_PARAMS::DIST) = 30.2;  // 45km
      bconfig.set_stage(TOTAL) = true;
      bconfig.control(CONTROLS::SEQUENCES) = InitSolver::setup_first;

      if (InitSolver::minimal_config) {
        bconfig.write_minimal_config();
      } else {
        solvers::bns_xcts_setup_bin_config(bconfig);
        bconfig.write_config();
      }
    }
  } else {
    bconfig.open_config();
    bconfig.control(CONTROLS::SEQUENCES) = InitSolver::setup_first;
    EOS_initialize::init(bconfig, NODES::BCO1);

    auto tree = bconfig.get_config_tree();
    auto N = solvers::number_of_sequences_binary(tree);
    if (N > 1) {
      if (rank == 0) {
        std::cerr << "Only sequences along one component is allowed.\n";
        std::_Exit(EXIT_FAILURE);
      }
    }

    auto resolution =
        solvers::parse_seq_tree(tree, "binary", "res", BIN_PARAMS::BIN_RES);
    solvers::verify_resolution_sequence(bconfig, resolution);

    auto seq = solvers::find_sequence_binary(tree);
    auto seq_bin = solvers::find_sequence(tree, MBIN_PARAMS, "binary");

    if (!(seq.is_set() || seq_bin.is_set()) &&
        !bconfig.control(CONTROLS::SEQUENCES))
      err = bns_xcts_driver(bconfig, resolution, InitSolver::outputdir);
    else {

      auto [branch_name, key, val] = find_leaf(tree, "N");
      if (!key.empty()) {
        seq.set_N(std::stoi(val));
        seq_bin.set_N(std::stoi(val));
      }
      if (seq.is_set()) {
        if (rank == 0) {
          std::cout << stdio_header("Companion Sequence", '-') << std::endl
                    << seq << std::endl;
        }
        err = solvers::bns_xcts_sequence(bconfig, seq, resolution,
                                         InitSolver::outputdir);
      } else {
        if (rank == 0) {
          std::cout << stdio_header("Binary Sequence", '-') << std::endl
                    << seq_bin << std::endl;
        }
        err = solvers::bns_xcts_sequence(bconfig, seq_bin, resolution,
                                         InitSolver::outputdir);
      }
    }
  }

  MPI_Finalize();
  return err;
}
