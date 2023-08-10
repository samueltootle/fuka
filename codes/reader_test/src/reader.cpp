/*
 * Copyright 2021
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
// #include "kadath_adapted_polar.hpp"
// #include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
// #include "bco_utilities.hpp"
// #include "EOS/EOS.hh"
// #include "coord_fields.hpp"
// #include "mpi.h"
#include "Solvers/reader.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_bh.hpp"

using namespace Kadath;
using namespace Kadath::FUKA_Config;

int main(int argc, char **argv) {
  // // initialize MPI
  // int rc = MPI_Init(&argc, &argv);

  // if (rc != MPI_SUCCESS) {
  //   cerr << "Error starting MPI" << endl;
  //   MPI_Abort(MPI_COMM_WORLD, rc);
  // }

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::string ifilename{argv[1]};
  kadath_config_boost<BCO_BH_INFO> bconfig(ifilename);
  Kadath::FUKA_Solvers::CFMS_BH_Reader<decltype(bconfig), Space_adapted_bh> reader(ifilename);

  // Space_adapted_bh bh(*reader.get_space().get());
  auto reader_copy(reader);
  
  
  // MPI_Finalize();
  return EXIT_SUCCESS;
}
