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
  
  template<class eos_t>
  using reader_t = Kadath::FUKA_Solvers::CFMS_BNS_Exporter<eos_t>;

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
  
  // get const EOS information - used for initializing EOS later
  const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT, NODES::BCO1);
  const std::string eos_file = bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE, NODES::BCO1);
  const std::string eos_type = bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE, NODES::BCO1);
  if(eos_type == "Cold_Table") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(EOS_PARAMS::INTERP_PTS, NODES::BCO1) == 0) ? \
                            2000 : bconfig.eos<int>(EOS_PARAMS::INTERP_PTS, NODES::BCO1);

    EOS<eos_t, eos_var_t::PRESSURE>::init(eos_file, h_cut, interp_pts);
    reader_t<eos_t> input_reader(ifilename);
    interp_data(input_reader, xx, yy, zz);
  }

  if(eos_type == "Cold_PWPoly") {
    using namespace Kadath::Margherita;
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,eos_var_t::PRESSURE>::init(eos_file, h_cut);
    reader_t<eos_t> input_reader(ifilename);
    interp_data(input_reader, xx, yy, zz);
  } // end adding EOS OPEs

  return EXIT_SUCCESS;
}
