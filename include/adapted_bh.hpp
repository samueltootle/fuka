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

#ifndef __ADAPTED_BH_HPP_
#define __ADAPTED_BH_HPP_

#include "space.hpp"
#include "term_eq.hpp"
#include "spheric.hpp"
#include "metric.hpp"
#include "adapted.hpp"
#include "homothetic.hpp"
#include <vector>

namespace Kadath {


/**
 * The \c Space_spheric_adapted class fills the space with one shell adapted on the inside, several standard shells and a compactified domain, all centered on the same point.
 * \ingroup domain
 */
class Space_adapted_bh : public Space {
     public:
	/**
     	* Standard constructor ; all the shells are initially not deformed.
     	* @param ttype [input] : the type of basis.
	* @param cr [input] : absolute coordinates of the center.
	* @param nbr [input] : number of points in each domain.
	* @param bounds [input] : radii of the various shells (and also determines the total number of domains).
	*/
	Space_adapted_bh (int ttype, const Point& cr, const Dim_array& nbr, const std::vector<double>& BH_bounds) ;
	Space_adapted_bh (FILE*) ; ///< Constructor from a file
  Space_adapted_bh (Space_adapted_bh const &) ;
	void add_eq (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ;
  void add_bc_bh  (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
  void add_bc_inf (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
  void add_eq_int_inf (System_of_eqs&, const char*);
  void add_eq_int_bh  (System_of_eqs&, const char*);
  void add_eq_zero_mode_inf (System_of_eqs&, const char*, int, int);

	virtual ~Space_adapted_bh() ; ///< Destructor        
	virtual void save(FILE*) const ;
	
	virtual int  nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;
	virtual Array<int> get_indices_matching_non_std(int, int) const ;
} ;
}
#endif
