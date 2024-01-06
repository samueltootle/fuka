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
#include "Solvers/ns_3d_xcts/ns_reader.hpp"
#include "kadath_adapted.hpp"
#ifdef _OPENMP
  #include <omp.h>
#endif
#include<thread>

using namespace Kadath;
using namespace Kadath::FUKA_Config;
  
  using config_t = kadath_config_boost<BCO_NS_INFO>;
  
  template<class eos_t>
  using reader_t = Kadath::FUKA_Solvers::CFMS_NS_Reader<eos_t, config_t, Space_spheric_adapted>;

constexpr unsigned int Npts = 1e5;
constexpr double range = 10;
constexpr double dx = range / Npts;

template<class reader_t>
void interp_data(reader_t& input_reader, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {
  std::vector<typename reader_t::pointwise_ary_t> all_data(Npts);

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

}

int main(int argc, char **argv) {

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::string ifilename{argv[1]};

  // ID Coords
  std::vector<double> xx(Npts);
  std::vector<double> yy(Npts);
  std::vector<double> zz(Npts);
  
  #pragma omp parallel for
  for(auto i = 0; i < Npts; ++i) {
    xx[i]+= i * dx;
    yy[i]+= i * dx;
    zz[i]+= i * dx;
  }
  // END ID Coords

  config_t bconfig(ifilename);
  
  // get const EOS information - used for initializing EOS later
  const double h_cut = bconfig.eos<double>(HCUT);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);
  if(eos_type == "Cold_Table") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                            2000 : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    reader_t<eos_t> input_reader(ifilename);
    interp_data(input_reader, xx, yy, zz);
  }

  if(eos_type == "Cold_PWPoly") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    reader_t<eos_t> input_reader(ifilename);
    interp_data(input_reader, xx, yy, zz);
  } // end adding EOS OPEs

  return EXIT_SUCCESS;
}
