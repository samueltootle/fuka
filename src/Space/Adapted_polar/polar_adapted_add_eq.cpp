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

#include "adapted_polar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {

void Space_polar_adapted::add_eq (System_of_eqs& sys, const char* eq, const char* rac, const char* rac_der, int nused, Array<int>** pused) const  {
	for (int dd=sys.get_dom_min() ; dd<sys.get_dom_max(); dd++) {
		sys.add_eq_inside (dd, eq, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac, nused, pused) ;
		sys.add_eq_matching (dd, OUTER_BC, rac_der, nused, pused) ;
	}
	sys.add_eq_inside (sys.get_dom_max(), eq, nused, pused) ;
}

void Space_polar_adapted::add_eq_int_inf (System_of_eqs& sys, const char* nom) {

	// Check the last domain is of the right type :
	const Domain_polar_compact* pcomp = dynamic_cast <const Domain_polar_compact*> (domains[nbr_domains-1]) ;
	if (pcomp==0x0) {
		cerr << "add_eq_int_inf requires a compactified domain" << endl ;
		abort() ;
	}
	int dom = nbr_domains-1 ;
    sys.eq_int_list.push_back(std::make_tuple(nom,dom,OUTER_BC));
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
		sys.eq_int[sys.neq_int] = new Eq_int(1) ;

		// Affectation :
		sys.eq_int[sys.neq_int]->set_part(0, new Ope_sub(&sys, sys.give_ope(dom, p1, OUTER_BC), sys.give_ope(dom, p2, OUTER_BC))) ;
		sys.neq_int ++ ;
	}
	sys.nbr_conditions = -1 ;
}
}
