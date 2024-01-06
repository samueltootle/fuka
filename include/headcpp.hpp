/*
    Copyright 2017 Philippe Grandclement

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __HEADCPP_HPP_
#define __HEADCPP_HPP_

#ifdef USE_CXX_STANDARD_17_OR_HIGHER
#define CXX_17_ATTRIBUTES(...) [[__VA_ARGS__]]
#else
#define CXX_17_ATTRIBUTES(...)
#endif

#include <iomanip>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cmath>
#include "stdlib.h"
#define PRECISION 1e-14

using namespace std ;

namespace Kadath {

  using std::cos ;
  using std::sin ;
  using std::sqrt ;
  using std::pow ;

  using std::log10 ;
  using std::log ;
  using std::exp ;
  using std::acos ;
  using std::asin ;

  using std::abs ;
  using std::fabs ;
  using std::max ;
  using std::min ;

  using std::cos;
  using std::sin;
  using std::tan ;
  using std::atan ;
  using std::atanh ; 
  using std::cosh;
  using std::sinh;

}
#endif
