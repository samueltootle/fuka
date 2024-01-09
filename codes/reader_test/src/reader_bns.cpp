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
#include "Configurator/config_binary.hpp"
#include "Solvers/bns_xcts/bns_exporter.hpp"
#include "kadath_adapted.hpp"
#include "kadath_bin_bh.hpp"

#include "./reader_test_tools.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif
#include<thread>

using namespace Kadath;
using namespace Kadath::FUKA_Config;
  
  using config_t = kadath_config_boost<BIN_INFO>;
  
  // template<class eos_t>
  using reader_t = Kadath::FUKA_Solvers::CFMS_BNS_Exporter;

constexpr unsigned int Npts = 1e4;  
constexpr double range = 100;
constexpr double dx = range / Npts;

int main(int argc, char **argv) {

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.BBH.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::string ifilename{argv[1]};


  std::vector<double> xx(Npts);
  std::vector<double> yy(Npts);
  std::vector<double> zz(Npts);
  
  #pragma omp parallel for
  for(auto i = 0; i < Npts; ++i) {
    xx[i]+= i * dx;
    yy[i]+= i * dx;
    zz[i]+= i * dx;
  }

  config_t bconfig(ifilename);
  reader_t input_reader(ifilename);
  interp_data(input_reader, xx, yy, zz);

  return EXIT_SUCCESS;
}
