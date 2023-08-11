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

constexpr unsigned int Npts = 128;  
constexpr double range = 3;
constexpr double dx = range / Npts;

void fill_coords(std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {
  unsigned int nthreads = std::thread::hardware_concurrency();
  unsigned int chunksize = Npts / nthreads;

  auto fill_points = [&](auto start, auto stop) {
    for(auto i = start; i < stop; ++i) {
      xx[i]+= i * dx;
      yy[i]+= i * dx;
      zz[i]+= i * dx;
    }
  };

  std::vector<std::thread> threads;
  for(auto j = 0; j < nthreads; ++j) {
    auto start = chunksize * j;
    auto stop = (chunksize * (j+1) > Npts) ? Npts : chunksize * (j+1);
    threads.push_back(std::thread(fill_points, start, stop));
  }
  for(auto j = 0; j < nthreads; ++j) {
    threads[j].join();
  }
}

ary_t interpolate(std::string fn, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {

  
  config_t bconfig(fn);  
  reader_t input_reader(fn);

  unsigned int nthreads = std::thread::hardware_concurrency();
  unsigned int chunksize = Npts / nthreads;

  std::vector<ary_t> all;

  auto interp_points = [&](auto start, auto stop) {
    for(auto i = start; i < stop; ++i) {
      reader_t reader(input_reader);
      auto vals = reader.export_pointwise(xx[i], yy[i], zz[i]);
      all.push_back(vals);
      // xx[i]+= i * dx;
      // yy[i]+= i * dx;
      // zz[i]+= i * dx;
    }
  };

  std::vector<std::thread> threads;
  for(auto j = 0; j < nthreads; ++j) {
    auto start = chunksize * j;
    auto stop = (chunksize * (j+1) > Npts) ? Npts : chunksize * (j+1);
    threads.push_back(std::thread(interp_points, start, stop));
  }
  for(auto j = 0; j < nthreads; ++j) {
    threads[j].join();
  }
}

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

  fill_coords(xx,yy,zz);
  auto v = interpolate(ifilename, xx, yy, zz);


  
  // #pragma omp parallel for
  // for(auto i = 0; i < Npts; ++i) {
  //   xx[i]+= i * dx;
  //   yy[i]+= i * dx;
  //   zz[i]+= i * dx;
  // }

  // #pragma omp parallel for firstprivate(reader)
  // for(auto i = 0; i < Npts; ++i) {
  //   reader_t r(reader);

  //   // r = reader;
  //   // #pragma omp critical 
  //   // {
  //   //   r = reader;
  //   // std::cout << i << "\n";
  //   std::cout << r.get_space() << "\n";
  //   auto interp = r.interpolate_pointwise(xx[i], yy[i], zz[i]);
  //   // std::cout << "Psi(" << xx[i] << ", " << yy[i] << ", " << zz[i] << ") = " << interp[reader_t::XCTS_VARS::XCTS_PSI] << '\n';
    
    

  // }

  return EXIT_SUCCESS;
}
