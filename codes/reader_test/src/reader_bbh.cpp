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
#include "Configurator/config_bco.hpp"
#include "Solvers/bbh_xcts/bbh_exporter.hpp"
#include "kadath_adapted.hpp"
#include "kadath_bin_bh.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif
#include<thread>

using namespace Kadath;
using namespace Kadath::FUKA_Config;
  
  using config_t = kadath_config_boost<BIN_INFO>;
  using reader_t = Kadath::FUKA_Solvers::CFMS_BBH_Exporter;
  using ary_t = std::vector<reader_t::output_ary_t>;

constexpr unsigned int Npts = 256;  
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
  ary_t all_data(Npts);  

  #pragma omp parallel for firstprivate(input_reader)
  for(auto i = 0; i < Npts; ++i) {
    all_data[i] = input_reader.export_pointwise(xx[i], yy[i], zz[i]);
  }

  for(auto i = 0; i < Npts; ++i) {
    std::cout << "(" << xx[i] << ", " << yy[i] << ", " << zz[i] << ") - ";
    
    for(auto& e : all_data[i])
      cout << e << " - ";
    cout << endl;
  }
  cout << xx[200] << ", " << yy[200] << ", " << zz[200] << ", " 
    << all_data[200][reader_t::OUTPUT_VARS::ALPHA] << "\n\t"
    << all_data[200][reader_t::OUTPUT_VARS::KXY] << ", "
    << all_data[200][reader_t::OUTPUT_VARS::KXZ] << ", "
    << all_data[200][reader_t::OUTPUT_VARS::KYZ] << "\n";

  return EXIT_SUCCESS;
}
