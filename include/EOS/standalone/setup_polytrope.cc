//
// This file is part of Margherita, the light-weight EOS framework
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

/********************************
 * Setup piecewise polytropes
 *
 * Written in June 2016 by Elias R. Most
 * <emost@th.physik.uni-frankfurt.de>
 *
 * Revised in May 2021 by Samuel D. Tootle
 * <tootle@itp.uni-frankfurt.de>
 *
 * -revision <28.May.21>
 *  Updated to be purely standalone
 *
 * We do things here as in Whisky_Exp
 * (See Takami et al. https://arxiv.org/pdf/1412.3240v2.pdf)
 * (also see Read et al. https://arxiv.org/pdf/0812.2163v1.pdf)
 ********************************/
#pragma once
#define PWPOLY_SETUP

#include "cold_pwpoly.hh"
#include "cold_pwpoly_implementation.hh"
#include "polytrope_file_parser.hpp"

#include "Margherita_EOS.h"
#include "margherita.hh"
#include <fstream>
#include <string>
#include <vector>

namespace Kadath {
namespace Margherita {

// From a parsed polytrope file, populate the Margherita
// piecewise polytrope implementation
inline void populate_Margherita_polytrope(std::string polytrope_file) {
  using namespace Kadath::FUKA_EOS;
  parse_polytrope_file parser;
  parser(polytrope_file);
  Cold_PWPoly::num_pieces = parser.num_pieces;
  Cold_PWPoly::rhomin = parser.rhomin;
  Cold_PWPoly::rhomax = parser.rhomax;
  Cold_PWPoly::k_tab[0] = parser.K0;
  Cold_PWPoly::P_tab[0] = parser.Pmin;
  for (int i = 0; i < Cold_PWPoly::num_pieces; ++i) {
    Cold_PWPoly::gamma_tab[i] = parser.gamma_tab[i];
    Cold_PWPoly::rho_tab[i] = parser.rho_tab[i];
  }
}

inline void Margherita_setup_polytrope(std::string polytrope_file) {
  using namespace Margherita_constants;
  populate_Margherita_polytrope(polytrope_file);

  //eps_tab == continuity coefficients on the bounds between
  //pieces.  The first is always 0.
  Cold_PWPoly::eps_tab[0] = 0.0;

  // Setup piecewise polytrope
  for (int i = 1; i < Cold_PWPoly::num_pieces; ++i) {
    // Consistency check
    if (Cold_PWPoly::rho_tab[i] <= Cold_PWPoly::rho_tab[i-1]) {
      std::cerr << "rho_tab must increase monotonically!\n";
			std::_Exit(EXIT_FAILURE);
    }
  	const double gam_im1 = Cold_PWPoly::gamma_tab[i] - 1.0;

    Cold_PWPoly::k_tab[i] =
        Cold_PWPoly::k_tab[i - 1] *
        pow(Cold_PWPoly::rho_tab[i],
            Cold_PWPoly::gamma_tab[i - 1] - Cold_PWPoly::gamma_tab[i]);

    Cold_PWPoly::eps_tab[i] =
        Cold_PWPoly::eps_tab[i - 1] +
        Cold_PWPoly::k_tab[i - 1] *
            pow(Cold_PWPoly::rho_tab[i], Cold_PWPoly::gamma_tab[i - 1] - 1.0) /
            (Cold_PWPoly::gamma_tab[i - 1] - 1.0) -
        Cold_PWPoly::k_tab[i] *
            pow(Cold_PWPoly::rho_tab[i], gam_im1) /  (gam_im1);
    Cold_PWPoly::P_tab[i] =
        Cold_PWPoly::k_tab[i] *
        pow(Cold_PWPoly::rho_tab[i], Cold_PWPoly::gamma_tab[i]);

    double eps = Cold_PWPoly::eps_tab[i] +
        Cold_PWPoly::P_tab[i] / Cold_PWPoly::rho_tab[i] / gam_im1;

    Cold_PWPoly::h_tab[i] = 1. + eps + Cold_PWPoly::P_tab[i] / Cold_PWPoly::rho_tab[i];
  }
  #ifdef DEBUG
  for (int nn = 0; nn < Cold_PWPoly::num_pieces; ++nn) {
    printf("nn= %d : K=%.15e , rho= %.15e, gamma= %.15e, eps= %.15e, P=.%15e, h=.%15e \n", nn,
        Cold_PWPoly::k_tab[nn], Cold_PWPoly::rho_tab[nn],
        Cold_PWPoly::gamma_tab[nn], Cold_PWPoly::eps_tab[nn],
        Cold_PWPoly::P_tab[nn], Cold_PWPoly::h_tab[nn]);
  }
  #endif

  return;
}
}
}
