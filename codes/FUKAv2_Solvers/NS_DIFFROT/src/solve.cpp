/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
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
#include "mpi.h"
// #include "FUKA_Solvers/NS_XCTS.hpp"
// #include "Solvers/ns_3d_xcts/ns_3d_xcts_driver.hpp"
#include "FUKA_Solvers/NS_XCTS_driver.cpp"
#include "Solvers/solver_startup.hpp"
#include "Solvers/sequences/ns_sequence.hpp"
#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/sequences/sequence_utilities.hpp"

using namespace Kadath::FUKA_Config;
using namespace Kadath::FUKA_Solvers;

int main(int argc, char** argv) {
  int rc = MPI_Init(&argc, &argv) ;
  if (rc!=MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl ;
    MPI_Abort(MPI_COMM_WORLD, rc) ;
  }
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

  using config_t = kadath_config_boost<BCO_NS_INFO>;
  using InitSolver = Initialize_Solver<config_t>;
  
  // Initialize static member variables
  InitSolver::input_configname = "initial_ns.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;

  if(InitSolver::example_setup) {
    if(rank == 0) {
      // Generate <example name>.info and terminate      
      if(InitSolver::minimal_config) {
        bconfig.set_minimal_defaults();
        bconfig.write_minimal_config();
      } else {
        bconfig.set_defaults();
        bconfig.control(CONTROLS::SEQUENCES) = InitSolver::setup_first;
        bconfig.write_config();
      }
    }
  } else {
    // We now have to assume bconfig is a minimal config
    // that contains sequences _init/_final
    bconfig.open_config();
    bconfig.control(CONTROLS::SEQUENCES) = InitSolver::setup_first;

    Tree tree;
    pt::read_info(bconfig.config_filename_abs(), tree);
    auto N = number_of_sequences(tree, MBCO_PARAMS, "ns");
    if(N > 1) {
      if(rank == 0) {
        std::cerr << "Only sequences along one component is allowed.\n";
        std::_Exit(EXIT_FAILURE);
      }
    }

    Parameter_sequence<BCO_PARAMS> resolution = parse_seq_tree(tree, "ns", "res", BCO_PARAMS::BCO_RES);
    verify_resolution_sequence(bconfig, resolution);

    ns_sequence seq = find_ns_sequence(tree);
    if(seq.is_set()) {
      auto [ branch_name, key, val ] = find_leaf(tree, "N");
      if(!key.empty()) seq.set_N(std::stoi(val));
    }

    verify_ns_fixing_values(bconfig, seq);

    if(!seq.is_set() && !bconfig.control(CONTROLS::SEQUENCES)) {
      initialize_config_from_fixing_values(bconfig, seq);
    }
    ns_3d_xcts_driver(bconfig, seq, resolution, InitSolver::outputdir);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
