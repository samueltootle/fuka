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
#pragma once
#include <string>
#include "Solvers/solvers.hpp"
namespace Kadath {
namespace FUKA_Solvers {
/**
 * \addtogroup Solver_startup
 * \ingroup FUKA
 * @{*/
/**
 * Initialize_Solver
 *
 * Static class to handle the boiler plate code at the start of all
 * the v2 solvers.
 *
 */
template <class config_t>
class Initialize_Solver {
  // static member variables
 public:
  static std::string outputdir;
  static std::string input_configname;
  static bool example_setup;
  static bool setup_first;
  static bool minimal_config;
  static config_t bconfig;
  static int rank;

 public:
  static void init_solver(int argc, char** argv);
};

// static member initialization
template <class config_t>
std::string Initialize_Solver<config_t>::outputdir = "./";

template <class config_t>
std::string Initialize_Solver<config_t>::input_configname;

template <class config_t>
bool Initialize_Solver<config_t>::example_setup = false;

template <class config_t>
bool Initialize_Solver<config_t>::setup_first = false;

template <class config_t>
bool Initialize_Solver<config_t>::minimal_config = true;

template <class config_t>
config_t Initialize_Solver<config_t>::bconfig;

template <class config_t>
int Initialize_Solver<config_t>::rank = 0;
// end static member initialization

template <class config_t>
void Initialize_Solver<config_t>::init_solver(int argc, char** argv) {
  // Necessary to get default construction
  Initialize_Solver<config_t>::bconfig = config_t();
  if (argc < 2) {
    if (Initialize_Solver::rank == 0) {
      std::cout << "INFO config file missing - generating a minimal setup"
                << std::endl
                << "Modify as needed before rerunning `solve "
                << Initialize_Solver::input_configname << " <outputdir>`\n"
                << "For a full setup, rerun with `solve full`\n";
    }
    Initialize_Solver::example_setup = Initialize_Solver::setup_first = true;
  }

  // if the input config is called something other than the default filename
  // check if we have a dat file associated with it.  If not, create a dat
  // file based on the info file before solving.
  else if (Initialize_Solver::input_configname != std::string{argv[1]}) {
    std::string arg1 = std::string{argv[1]};
    if (arg1 == "full") {
      Initialize_Solver<config_t>::minimal_config = false;
      if (Initialize_Solver::rank == 0) {
        std::cout << "INFO config file missing - generating a full setup"
                  << std::endl
                  << "Modify as needed before rerunning `solve "
                  << Initialize_Solver::input_configname << " <outputdir>`\n";
      }
      Initialize_Solver::example_setup = Initialize_Solver::setup_first = true;
    } else {
      Initialize_Solver::input_configname = arg1;
      Initialize_Solver::bconfig.set_filename(
          Initialize_Solver::input_configname);

      std::string input_setupname = Initialize_Solver::bconfig.space_filename();
      std::string file_path = Initialize_Solver::bconfig.config_outputdir();

      if (fs::exists(Initialize_Solver::input_configname) &&
          fs::exists(input_setupname)) {
        if (Initialize_Solver::rank == 0)
          std::cout << "Solving based on previous solution: "
                    << Initialize_Solver::input_configname << std::endl;
      } else if (fs::exists(Initialize_Solver::input_configname) &&
                 !fs::exists(input_setupname)) {
        if (Initialize_Solver::rank == 0) {
          std::cout << "No dat file associated with : "
                    << Initialize_Solver::input_configname << std::endl;
          std::cout << "Creating a setup based on config file..." << std::endl;
        }
        Initialize_Solver::setup_first = true;
      } else {
        if (Initialize_Solver::rank == 0)
          std::cerr << "Config file input is invalid.  Check for typos."
                    << std::endl;
        std::_Exit(EXIT_FAILURE);
      }
    }
  } else {
    // if the default <example name>.info file is used, we always start over.
    Initialize_Solver::setup_first = true;
  }

  if (argc > 2)
    Initialize_Solver::outputdir = std::string{argv[2]};

  Initialize_Solver::bconfig.set_filename(Initialize_Solver::input_configname);
}
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath