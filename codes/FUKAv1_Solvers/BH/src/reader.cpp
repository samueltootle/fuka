/*
 * Copyright 2021
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
#include "kadath_adapted_bh.hpp"
#include "kadath.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include <iostream>
#include <memory>
#include "mpi.h"
using namespace Kadath;
using namespace Kadath::FUKA_Config;

int main(int argc, char **argv) {
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  if(argc < 2) {
    std::cerr << "Usage: ./reeader /<path>/<ID base name>.info " << std::endl;
    std::cerr << "Ex: ./reader converged.BH.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }
  cout.precision(14);

  //Name of config.info file
  std::string input_filename = argv[1];
  kadath_config_boost<BCO_BH_INFO> bconfig(input_filename);  
  std::string kadath_filename = 
    input_filename.substr(0,input_filename.size()-5)+".dat";
  
  FILE *fich = fopen(kadath_filename.c_str(), "r");

  Space_adapted_bh space(fich);
  Scalar conf(space, fich);
  Scalar lapse(space, fich);
  Vector shift(space, fich);
  fclose(fich);

	int ndom = space.get_nbr_domains() ;
	Base_tensor basis (space, CARTESIAN_BASIS) ;
	Metric_flat fmet (space, basis) ;

  double r = bco_utils::get_radius(space.get_domain(1), OUTER_BC);
 	double xo = bco_utils::get_center(space,0);

	CoordFields<Space_adapted_bh> cf_generator(space);
  vec_ary_t coord_vectors {default_co_vector_ary(space)};

  update_fields_co   (cf_generator, coord_vectors, {}, xo);
	
  System_of_eqs syst  (space  , 0, ndom-1) ;

	fmet.set_system     (syst   , "f") ;

  syst.add_cst        ("PI"   , M_PI) ;

	syst.add_cst        ("M"    , bconfig(MIRR)) ;
	syst.add_cst        ("chi"  , bconfig(CHI)) ;
  syst.add_cst        ("CM"   , bconfig(MCH));

  syst.add_cst        ("mbh"  , *coord_vectors[GLOBAL_ROT]) ;
	syst.add_cst        ("sbh"  , *coord_vectors[S_BCO1]) ;

	syst.add_cst        ("ex"   , *coord_vectors[EX])  ;
	syst.add_cst        ("ey"   , *coord_vectors[EY])  ;
	syst.add_cst        ("ez"   , *coord_vectors[EZ])  ;
	syst.add_cst        ("einf" , *coord_vectors[S_INF]) ;

  syst.add_cst        ("ome"  , bconfig(OMEGA)) ;
	syst.add_cst        ("P"    , conf) ;
	syst.add_cst        ("N"    , lapse) ;
	syst.add_cst        ("bet"  , shift) ;

  syst.add_def        ("NP = P*N");
  syst.add_def        ("Ntilde = N / P^6");
  syst.add_def        ("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

  syst.add_def        ("intPx = A_ij * ex^j * einf^i / 8 / PI") ;
  syst.add_def        ("intPy = A_ij * ey^j * einf^i / 8 / PI") ;
  syst.add_def        ("intPz = A_ij * ez^j * einf^i / 8 / PI") ;

  syst.add_def        ("intS = A_ij * mbh^i * sbh^j / 8. / PI") ;
  syst.add_def        ("intA = P^4") ;
  syst.add_def        ("intMsq = intA / 16. / PI") ;

  syst.add_def        (ndom - 1, "intJ = multr(A_ij * mbh^j * einf^i) / 8 / PI");
  syst.add_def        (ndom - 1, "Madm = -dr(P) / 2 / PI");
  syst.add_def        (ndom - 1, "Mk   =  dr(N) / 4 / PI");

  // ADM mass :
  double adm_inf  = space.get_domain(ndom-1)->integ(syst.give_val_def("Madm")()(ndom-1) , OUTER_BC);

  // Komar mass :
  double komar    = space.get_domain(ndom-1)->integ(syst.give_val_def("Mk")()(ndom-1)   , OUTER_BC);

  // Computation of the irreducible mass
  double mirrsq   = space.get_domain(2)->integ(syst.give_val_def("intMsq")()(2)    , INNER_BC);
  double mirr     = std::sqrt(mirrsq);
  
  // Computation of local angular momentum and Mch
  Val_domain integS(syst.give_val_def("intS")()(2));
  double S        = space.get_domain(2)->integ(integS, INNER_BC);
  double Mch      = std::sqrt( mirrsq + S * S / 4. / mirrsq );
  double Chi      = S / Mch / Mch;

  double A = space.get_domain(2)->integ(syst.give_val_def("intA")()(2), INNER_BC);
  double AR = sqrt(A / 4. / acos(-1.));

  // J at infinity :
  double Jinf     = space.get_domain(ndom-1)->integ(syst.give_val_def("intJ")()(ndom-1) , OUTER_BC);
  
  auto [ lapsemin, lapsemax ] = bco_utils::get_field_min_max(lapse, 2, INNER_BC);
  auto [ confmin , confmax  ] = bco_utils::get_field_min_max(conf , 2, INNER_BC);

  #define FORMAT1 std::setw(25) << std::right << std::setprecision(5) << std::fixed << std::showpos
  auto print_shells = [&](int dom_min, int dom_max)
  {
    int cnt = 1;
    for(int i = dom_min; i < dom_max; ++i) {
      std::string shell{"SHELL"+std::to_string(cnt)+" = "};
      std::cout << FORMAT1 << shell << bco_utils::get_radius(space.get_domain(i), OUTER_BC) << std::endl;
      cnt++;
    }
  };
  auto res_r = space.get_domain(0)->get_nbr_points()(0);
  auto res_t = space.get_domain(0)->get_nbr_points()(1);
  auto res_p = space.get_domain(0)->get_nbr_points()(2);

  std::cout << FORMAT1 << "RES = " << "[" << res_r << "," << res_t << "," << res_p << "]\n";
  std::cout << FORMAT1 << "Coord R_IN = "      << bco_utils::get_radius(space.get_domain(0), OUTER_BC) << std::endl
            << FORMAT1 << "Coord R = " << bco_utils::get_radius(space.get_domain(1), OUTER_BC) << std::endl;
  print_shells(2, ndom-2);
  std::cout << FORMAT1 << "Coord R_OUT = "     << bco_utils::get_radius(space.get_domain(ndom-2), OUTER_BC) << "\n\n"
            << FORMAT1 << "Areal R = " << AR << std::endl
            << FORMAT1 <<" LAPSE = "       << "[" << lapsemin << ", " << lapsemax  <<"]\n"
            << FORMAT1 <<" PSI = "         << "[" << confmin  << ", " << confmax   <<"]\n\n";
  
  #ifdef FORMAT
    #undef FORMAT
  #endif

  #define FORMAT std::setw(25) << std::right << std::setprecision(5) << std::scientific << std::showpos
  std::cout << FORMAT1 << "Mirr = "        << mirr << "[" << bconfig(MIRR) << "]\n"
            << FORMAT1 << "Mch = "         << Mch  << "[" << bconfig(MCH) << "]\n"
            << FORMAT1 << "Chi = "         << Chi  << "[" << bconfig(CHI) << "]\n\n"
            << FORMAT << "Omega = "       << bconfig(OMEGA) << std::endl
            << FORMAT << "Komar mass = "   << komar          << std::endl
            << FORMAT << "Adm mass = "   << adm_inf
                      << ", Diff: "        << 2. * fabs(adm_inf-komar)/(adm_inf+komar) << std::endl
            << FORMAT << "Adm moment. = "  << Jinf << std::endl
            << FORMAT << "Px = "           << space.get_domain(ndom - 1)->integ(syst.give_val_def("intPx")()(ndom-1), OUTER_BC)  << std::endl
            << FORMAT << "Py = "           << space.get_domain(ndom - 1)->integ(syst.give_val_def("intPy")()(ndom-1), OUTER_BC)  << std::endl
            << FORMAT << "Pz = "           << space.get_domain(ndom - 1)->integ(syst.give_val_def("intPz")()(ndom-1), OUTER_BC)  << std::endl;

  MPI_Finalize();
  return EXIT_SUCCESS;
}
