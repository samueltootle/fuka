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

#include "headcpp.hpp"
#include "bin_ns.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"

namespace Kadath {

void Space_bin_ns::add_bc_sphere_one (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	  sys.add_eq_bc (ADAPTED1+1, INNER_BC, name, nused, pused) ;
}

void Space_bin_ns::add_bc_sphere_two (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	sys.add_eq_bc (ADAPTED2+1, INNER_BC, name, nused, pused) ;
}

void Space_bin_ns::add_bc_outer (System_of_eqs& sys, const char* name, int nused, Array<int>** pused)  {
	sys.add_eq_bc (OUTER+n_shells_outer+5, OUTER_BC, name, nused, pused) ;
}

void Space_bin_ns::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
	// First NS
	sys.add_eq_inside (NS1, eq, nused, pused) ;
	sys.add_eq_matching (NS1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (NS1, OUTER_BC, rac_der, nused, pused) ;

	for(int i = 0; i < n_inner_shells1; ++i) {
		sys.add_eq_inside (NS1+1+i, eq, nused, pused) ;
		sys.add_eq_matching (NS1+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (NS1+1+i, OUTER_BC, rac_der, nused, pused) ;
	}

	sys.add_eq_inside   (ADAPTED1  , eq, nused, pused) ;
	sys.add_eq_matching (ADAPTED1  , OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (ADAPTED1  , OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside   (ADAPTED1+1, eq, nused, pused) ;

	for(int i = 0; i < n_shells1; ++i) {
		sys.add_eq_matching (ADAPTED1+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (ADAPTED1+1+i, OUTER_BC, rac_der, nused, pused) ;
		sys.add_eq_inside   (ADAPTED1+2+i, eq, nused, pused) ;
	}

	// Matching with bispheric :
	sys.add_eq_matching_import (ADAPTED1+1+n_shells1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (OUTER     , INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (OUTER+1   , INNER_BC, rac_der, nused, pused) ;

	// Second NS :
	sys.add_eq_inside   (NS2, eq, nused, pused) ;
	sys.add_eq_matching (NS2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (NS2, OUTER_BC, rac_der, nused, pused) ;

	for(int i = 0; i < n_inner_shells2; ++i) {
		sys.add_eq_inside (NS2+1+i, eq, nused, pused) ;
		sys.add_eq_matching (NS2+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (NS2+1+i, OUTER_BC, rac_der, nused, pused) ;
	}

	sys.add_eq_inside   (ADAPTED2  , eq, nused, pused) ;
	sys.add_eq_matching (ADAPTED2  , OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (ADAPTED2  , OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside   (ADAPTED2+1, eq, nused, pused) ;

	for(int i = 0; i < n_shells2; ++i) {
		sys.add_eq_matching (ADAPTED2+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (ADAPTED2+1+i, OUTER_BC, rac_der, nused, pused) ;
		sys.add_eq_inside   (ADAPTED2+2+i, eq, nused, pused) ;
	}

	// Matching with bispheric :
	sys.add_eq_matching_import (ADAPTED2+1+n_shells2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (OUTER+3   , INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (OUTER+4   , INNER_BC, rac_der, nused, pused) ;

	// Chi first
	sys.add_eq_inside   (OUTER, eq, nused, pused) ;
	sys.add_eq_matching (OUTER, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside   (OUTER+1, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+1, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+1, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside   (OUTER+2, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+2, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+2, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect
	sys.add_eq_inside   (OUTER+3, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+3, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+3, CHI_ONE_BC, rac_der, nused, pused) ;


	// chi first :
	sys.add_eq_inside   (OUTER+4, eq, nused, pused) ;

	// Matching outer domain :
	for (int d=OUTER ; d<=OUTER+4 ; d++)
		sys.add_eq_matching_import (d, OUTER_BC, rac, nused, pused) ;

  //Matching for first shell or compactified domain
	sys.add_eq_matching_import (OUTER+5, INNER_BC, rac_der, nused, pused) ;

  // Optional spherical shells between bi-spherical and compactified domain
  for (int d=0 ; d<n_shells_outer ; d++) {
      sys.add_eq_inside (OUTER+5+d, eq, nused, pused) ;
      sys.add_eq_matching (OUTER+5+d, OUTER_BC, rac, nused, pused) ;
      sys.add_eq_matching (OUTER+5+d, OUTER_BC, rac_der, nused, pused) ;
  }

  //Compactified domain
  sys.add_eq_inside (OUTER+5+n_shells_outer, eq, nused, pused) ;
}

void Space_bin_ns::add_eq_nozec (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {
  // FIXME not fixed yet

	// First NS
	sys.add_eq_inside (0, eq, nused, pused) ;
	sys.add_eq_matching (0, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (0, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (1, eq, nused, pused) ;
	sys.add_eq_matching (1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (1, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (2, eq, nused, pused) ;


	// Matching with bispheric :
	sys.add_eq_matching_import (2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (6, INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (7, INNER_BC, rac_der, nused, pused) ;

	// Second NS :
	sys.add_eq_inside (3, eq, nused, pused) ;
	sys.add_eq_matching (3, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (3, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (4, eq, nused, pused) ;
	sys.add_eq_matching (4, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (4, OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside (5, eq, nused, pused) ;

	// Matching with bispheric :
	sys.add_eq_matching_import (5, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (9, INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (10, INNER_BC, rac_der, nused, pused) ;

	// Chi first
	sys.add_eq_inside (6, eq, nused, pused) ;
	sys.add_eq_matching (6, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (6, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside (7, eq, nused, pused) ;
	sys.add_eq_matching (7, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (7, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside (8, eq, nused, pused) ;
	sys.add_eq_matching (8, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (8, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect
	sys.add_eq_inside (9, eq, nused, pused) ;
	sys.add_eq_matching (9, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (9, CHI_ONE_BC, rac_der, nused, pused) ;


	// chi first :
	sys.add_eq_inside (10, eq, nused, pused) ;

	if (n_shells_outer !=0) {
	// Matching outer domain :
	for (int d=6 ; d<=10 ; d++)
		sys.add_eq_matching_import (d, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (11, INNER_BC, rac_der, nused, pused) ;

	for (int d=0 ; d<n_shells_outer ; d++) {
		sys.add_eq_inside (11+d, eq, nused, pused) ;
		if (d!=n_shells_outer-1) {
			sys.add_eq_matching (11+d, OUTER_BC, rac, nused, pused) ;
			sys.add_eq_matching (11+d, OUTER_BC, rac_der, nused, pused) ;
		}
	}
	}
}

void Space_bin_ns::add_eq_nozec (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, const List_comp& used) {
	add_eq_nozec (sys, eq, rac,rac_der, used.get_ncomp(), used.get_pcomp()) ;
}

void Space_bin_ns::add_eq_noshell (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused)  {

	// First NS
	sys.add_eq_inside (NS1, eq, nused, pused) ;
	sys.add_eq_matching (NS1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (NS1, OUTER_BC, rac_der, nused, pused) ;

	for(int i = 0; i < n_inner_shells1; ++i) {
		sys.add_eq_inside (NS1+1+i, eq, nused, pused) ;
		sys.add_eq_matching (NS1+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (NS1+1+i, OUTER_BC, rac_der, nused, pused) ;
	}

	sys.add_eq_inside   (ADAPTED1  , eq, nused, pused) ;
	sys.add_eq_matching (ADAPTED1  , OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (ADAPTED1  , OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside   (ADAPTED1+1, eq, nused, pused) ;

	for(int i = 0; i < n_shells1; ++i) {
		sys.add_eq_matching (ADAPTED1+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (ADAPTED1+1+i, OUTER_BC, rac_der, nused, pused) ;
		sys.add_eq_inside   (ADAPTED1+2+i, eq, nused, pused) ;
	}

	// Matching with bispheric :
	sys.add_eq_matching_import (ADAPTED1+1+n_shells1, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (OUTER     , INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (OUTER+1   , INNER_BC, rac_der, nused, pused) ;

	// Second NS :
	sys.add_eq_inside   (NS2, eq, nused, pused) ;
	sys.add_eq_matching (NS2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (NS2, OUTER_BC, rac_der, nused, pused) ;

	for(int i = 0; i < n_inner_shells2; ++i) {
		sys.add_eq_inside (NS2+1+i, eq, nused, pused) ;
		sys.add_eq_matching (NS2+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (NS2+1+i, OUTER_BC, rac_der, nused, pused) ;
	}

	sys.add_eq_inside   (ADAPTED2  , eq, nused, pused) ;
	sys.add_eq_matching (ADAPTED2  , OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching (ADAPTED2  , OUTER_BC, rac_der, nused, pused) ;
	sys.add_eq_inside   (ADAPTED2+1, eq, nused, pused) ;

	for(int i = 0; i < n_shells2; ++i) {
		sys.add_eq_matching (ADAPTED2+1+i, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (ADAPTED2+1+i, OUTER_BC, rac_der, nused, pused) ;
		sys.add_eq_inside   (ADAPTED2+2+i, eq, nused, pused) ;
	}

	// Matching with bispheric :
	sys.add_eq_matching_import (ADAPTED2+1+n_shells2, OUTER_BC, rac, nused, pused) ;
	sys.add_eq_matching_import (OUTER+3   , INNER_BC, rac_der, nused, pused) ;
	sys.add_eq_matching_import (OUTER+4   , INNER_BC, rac_der, nused, pused) ;

	// Chi first
	sys.add_eq_inside (OUTER, eq, nused, pused) ;
	sys.add_eq_matching (OUTER, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER, CHI_ONE_BC, rac_der, nused, pused) ;

	// Rect :
	sys.add_eq_inside (OUTER+1, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+1, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+1, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Eta first
	sys.add_eq_inside (OUTER+2, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+2, ETA_PLUS_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+2, ETA_PLUS_BC, rac_der, nused, pused) ;

	// Rect
	sys.add_eq_inside (OUTER+3, eq, nused, pused) ;
	sys.add_eq_matching (OUTER+3, CHI_ONE_BC, rac, nused, pused) ;
	sys.add_eq_matching (OUTER+3, CHI_ONE_BC, rac_der, nused, pused) ;


	// chi first :
	sys.add_eq_inside (OUTER+4, eq, nused, pused) ;
}

void Space_bin_ns::add_eq_noshell (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, const List_comp& used) {
	add_eq_noshell (sys, eq, rac,rac_der, used.get_ncomp(), used.get_pcomp()) ;
}

void Space_bin_ns::add_eq_int_inf (System_of_eqs& sys, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, nbr_domains-1, OUTER_BC));

	// Check the last domain is of the right type :
	const Domain_compact* pcomp = dynamic_cast <const Domain_compact*> (domains[nbr_domains-1]) ;
	if (pcomp==0x0) {
		cerr << "add_eq_int_inf requires a compactified domain" << endl ;
		abort() ;
	}
	int dom = nbr_domains-1 ;

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

		// Construction of the equation
		sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		// Affectation :
		// no lhs :
		if (indic)
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(dom, p1, OUTER_BC)) ;

		else
			sys.eq_int[sys.neq_int]->set_part(0, new Ope_sub(&sys, sys.give_ope(dom, p1, OUTER_BC), sys.give_ope(dom, p2, OUTER_BC))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}

void Space_bin_ns::add_eq_int_volume (System_of_eqs& sys, int dmin, int dmax, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, dmin, -1));

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;

	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {
		// Construction of the equation
		int nz = dmax - dmin + 1;
		sys.eq_int[sys.neq_int] = new Eq_int(nz+1) ;

		// Affectation of the intregrale parts
		for (int d=0 ; d<nz ; d++)
		  sys.eq_int[sys.neq_int]->set_part(d, sys.give_ope(d+dmin, p1)) ;
		// Affectation of the second member (constant value)
		sys.eq_int[sys.neq_int]->set_part(nz, new Ope_minus(&sys, sys.give_ope(dmin, p2))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}

void Space_bin_ns::add_eq_ori_one (System_of_eqs& sys, const char* name) {

	Index pos (domains[NS1]->get_nbr_points()) ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val (NS1, auxi, pos) ;
}

void Space_bin_ns::add_eq_ori_two (System_of_eqs& sys, const char* name) {

	Index pos (domains[NS2]->get_nbr_points()) ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val (NS2, auxi, pos) ;
}

void Space_bin_ns::add_eq_int_outer_sphere_one (System_of_eqs& sys, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, ADAPTED1+1, OUTER_BC));

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {

		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

	  // Construction of the equation
	  sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		  // Affectation :
		  // no lhs :
		  if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(ADAPTED1+1, p1, OUTER_BC)) ;
			}
		  else {
			sys.eq_int[sys.neq_int]->set_part(0,
				new Ope_sub(&sys, sys.give_ope(ADAPTED1+1, p1, OUTER_BC), sys.give_ope(ADAPTED1+1, p2, OUTER_BC))) ;
		  }
		}
		sys.neq_int ++ ;
	sys.nbr_conditions = -1 ;
}

void Space_bin_ns::add_eq_int_outer_sphere_two (System_of_eqs& sys, const char* nom) {
  sys.eq_int_list.push_back(std::make_tuple(nom, ADAPTED2+1, OUTER_BC));

	// Get the lhs and rhs
	char p1[LMAX] ;
	char p2[LMAX] ;
	bool indic = sys.is_ope_bin(nom, p1, p2, '=') ;
	if (!indic) {
		cerr << "= needed for equations" << endl ;
		abort() ;
	}
	else {

		// Verif lhs = 0 ?
		indic = ((p2[0]=='0') && (p2[1]==' ') && (p2[2]=='\0')) ?
			true : false ;

	  // Construction of the equation
	  sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		  // Affectation :
		  // no lhs :
		  if (indic) {
			sys.eq_int[sys.neq_int]->set_part(0, sys.give_ope(ADAPTED2+1, p1, OUTER_BC)) ;
			}
		  else {
			sys.eq_int[sys.neq_int]->set_part(0,
				new Ope_sub(&sys, sys.give_ope(ADAPTED2+1, p1, OUTER_BC), sys.give_ope(ADAPTED2+1, p2, OUTER_BC))) ;
		  }
		}
		sys.neq_int ++ ;
	sys.nbr_conditions = -1 ;
}

}

