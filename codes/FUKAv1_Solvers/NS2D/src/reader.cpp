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
	Scalar nulogA   (space, ff1) ;
	Scalar nu (space, ff1) ;
  Scalar logh   (space, ff1) ;
  Scalar bet   (space, ff1) ;
  Scalar wrsint(space,ff1);
	fclose(ff1) ;

  // central values of the matter fields
  double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);
  double hc = std::exp(loghc);
  double nc = EOS<eos_t,DENSITY>::get(hc);
  double pc = EOS<eos_t,PRESSURE>::get(hc);

  // minimal and maximal radius of the adapted surface domain
  auto [ rmin, rmax ] = bco_utils::get_rmin_rmax(space, 1);
  // inner radius of the nucleus
  double rin1 = bco_utils::get_radius(space.get_domain(0), OUTER_BC);

  int ndom = space.get_nbr_domains();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  syst.add_cst("H", logh);
  syst.add_cst("nu", nu);
  syst.add_cst("nulogA", nulogA);
  syst.add_cst("bet", bet);
  syst.add_cst("wrsint", wrsint);

  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(nulogA - nu)");
  syst.add_def("B = (divrsint(bet) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");

  syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
  syst.add_def(ndom - 1, "intMadm2 = - (dr(A)) / 4piG ");
  syst.add_def(ndom - 1, "intMk = B * (dr(N) - multrsint(multrsint(B^2) / 2 / N * w * dr(w)))  / 4piG");
  syst.add_def(ndom - 1, "intJ = -multr(multrsint(dr(w)))  / 4/4piG");
 
 Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
  Val_domain integMadm2(syst.give_val_def("intMadm2")()(ndom - 1));
  double Madm2 = space.get_domain(ndom - 1)->integ(integMadm2, OUTER_BC);

  // Komar mass at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // ADM angular momentum at infinity 
  Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
  double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);
  
  #ifdef FORMAT
    #undef FORMAT
  #endif
  #define FORMAT std::setw(25) << std::right << std::setprecision(5) << std::fixed << std::showpos
  auto print_shells = [&](int dom_min, int dom_max)
  {
    int cnt = 1;
    for(int i = dom_min; i < dom_max; ++i) {
      std::string shell{"SHELL"+std::to_string(cnt)+" = "};
      std::cout << FORMAT << shell << bco_utils::get_radius(space.get_domain(i), OUTER_BC) << std::endl;
      cnt++;
    }
  };

  auto res_r = space.get_domain(0)->get_nbr_points()(0);
  auto res_t = space.get_domain(0)->get_nbr_points()(1);

  // output to stdout
  std::cout << FORMAT << "RES = "  << "[" << res_r << "," << res_t << "]\n"
            << FORMAT << "Coord R_IN = "  << rin1 << std::endl
            << FORMAT << "Coord R = "     << "[" << rmin << ", " << rmax << "]\n";
  std::cout << FORMAT << "Coord R_OUT = " << bco_utils::get_radius(space.get_domain(2), OUTER_BC) << "\n";
  print_shells(3, ndom-1); cout << endl;

  // std::cout << FORMAT << "Areal R = "    << AR << " [" << AR * M2km << "km]\n"
            // << FORMAT << "Baryonic Mass = " << baryonic_mass << std::endl
    std::cout \
            << FORMAT << "ADM Mass = " << Madm << " [" << Madm2 << "]\n"
            << FORMAT << "ADM Momentum = " << J << std::endl
            // << FORMAT << "Chi = " << J / Madm / Madm << " [" << bconfig(CHI) << "]\n"
            << FORMAT << "Omega = "<< bconfig(OMEGA) << std::endl
            << FORMAT << std::scientific << "Central Density = " << nc  << std::endl
            << FORMAT << std::scientific << "Central h = " << hc << std::endl
            << FORMAT << std::scientific << "Central log(h) = " << loghc << std::endl
            << FORMAT << std::scientific << "Central Pressure = " << pc << "\n\n";
            // << FORMAT << std::scientific << "Central dlog(h)/dx = " << central_dHdx << std::endl
            // << FORMAT << std::scientific << "Central Euler Constant = "<< central_euler << std::endl
            // << FORMAT << "Integrated log(h) = "    << H_integral << "\n\n";

  std::cout << FORMAT << "Mk = "   << Mk << std::scientific
            << ", Diff: " << 2. * fabs(Madm-Mk)/(Madm+Mk) << std::endl;
            // << FORMAT << "Px = "   << Px   << std::endl
            // << FORMAT << "Py = "   << Py   << std::endl
            // << FORMAT << "Pz = "   << Pz   << std::endl;
}
