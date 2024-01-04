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
#include "Solvers/reader.hpp"
// #include "Solvers/interpolator.hpp"
// #include "Solvers/bh_3d_xcts/bh_reader.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_bh.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif
#include<thread>

using namespace Kadath;
using namespace Kadath::FUKA_Config;
  
  using config_t = kadath_config_boost<BCO_BH_INFO>;
  using reader_t = Kadath::FUKA_Solvers::CFMS_BH_Reader<config_t, Space_adapted_bh>;
  using ary_t = std::array<double, reader_t::OUTPUT_VARS::NUM_OUTPUT_VARS>;

constexpr unsigned int Npts = 256;  
constexpr double range = 2;
constexpr double dx = range / Npts;

struct get_space_ptr {
  using space_t = Space_adapted_bh;
  using ptr_t = std::unique_ptr<space_t>;
  static ptr_t&& get(std::string filename) {
    
    FILE* ff1 = fopen (filename.c_str(), "r") ;
    space_t _space{ff1};
    fclose(ff1);
    static ptr_t p(new space_t(_space));
    return std::move(p);
  }
};

int main(int argc, char **argv) {

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
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
  auto space_ptr = get_space_ptr::get(bconfig.space_filename());
  std::vector<reader_t::pointwise_ary_t> all_data(Npts);
  // auto x = input_reader.export_pointwise(0.5, 0., 0.);
  // Kadath::FUKA_Solvers::Interpolator interp{input_reader, space_ptr};
  

  #pragma omp parallel for firstprivate(input_reader)
  for(auto i = 0; i < Npts; ++i) {
    all_data[i] = input_reader.export_pointwise(xx[i], yy[i], zz[i]);
  }

  for(auto i = 0; i < Npts; ++i)
    std::cout << "(" << xx[i] << ", " << yy[i] << ", " << zz[i] << ") - " << all_data[i][reader_t::OUTPUT_VARS::ALPHA] << " - "
              // << all_data[i][reader_t::OUTPUT_VARS::GXX] << ", "
              // << all_data[i][reader_t::OUTPUT_VARS::BETAX] << ", "
              // << all_data[i][reader_t::OUTPUT_VARS::BETAY] << ", "
              // << all_data[i][reader_t::OUTPUT_VARS::BETAZ] << "\n";
              << all_data[i][reader_t::OUTPUT_VARS::KXY] << ", "
              << all_data[i][reader_t::OUTPUT_VARS::KXZ] << ", "
              << all_data[i][reader_t::OUTPUT_VARS::KYZ] << "\n";
  cout << xx[200] << ", " << yy[200] << ", " << zz[200] << ", " 
    << all_data[200][reader_t::OUTPUT_VARS::ALPHA] << "\n\t"
    << all_data[200][reader_t::OUTPUT_VARS::KXY] << ", "
    << all_data[200][reader_t::OUTPUT_VARS::KXZ] << ", "
    << all_data[200][reader_t::OUTPUT_VARS::KYZ] << "\n";

  return EXIT_SUCCESS;
}
