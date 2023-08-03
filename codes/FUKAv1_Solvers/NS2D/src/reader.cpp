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
// FUKA includes
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
// #include "EOS/standalone/tov.hh"

// Kadath includes
#include "kadath.hpp"
#include "kadath_adapted_polar.hpp"
#include "kadath_adapted.hpp"

// C++ includes
#include <fstream>
#include <string>
#include <type_traits>
#include <memory>
#include <algorithm>

using namespace Kadath;
using namespace Kadath::FUKA_Config;

template<class eos_t, typename config_t>
void reader_2d(config_t bconfig);

// conversion from solar mass to km
constexpr double M2km = 1.4769994423016508;

int main(int argc, char **argv) {
  // initialize MPI
  //int rc = MPI_Init(&argc, &argv);

  //if (rc != MPI_SUCCESS) {
  //  cerr << "Error starting MPI" << endl;
  //  MPI_Abort(MPI_COMM_WORLD, rc);
  //}

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  // load the configuration
  std::string ifilename{argv[1]};
  kadath_config_boost<BCO_NS_INFO> bconfig(ifilename);

  // setup the EOS
  const double h_cut = bconfig.eos<double>(HCUT);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);

    // call reader to output diagnostics
    reader_2d<eos_t>(bconfig);
  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                            2000 : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);

    // call reader to output diagnostics
    reader_2d<eos_t>(bconfig);
  } else {
    std::cerr << "Unknown EOSTYPE." << endl;
    std::_Exit(EXIT_FAILURE);
  }

  //MPI_Finalize();
  return EXIT_SUCCESS;
}

template<class eos_t, typename config_t>
void reader_2d(config_t bconfig) {
  // load the space (and thus the domain setup)
  auto spacein = bconfig.space_filename();
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // load the fields defined on the space
	Scalar bigA   (space, ff1) ;
	Scalar lapse  (space, ff1) ;
  Scalar logh   (space, ff1) ;
	fclose(ff1) ;
  std::cout << space;
}
