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
#include "adapted_polar.hpp"
#include "polar.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"

namespace Kadath {
Space_polar_adapted::Space_polar_adapted (int ttype, const Point& center, const Dim_array& res, const Array<double>& bounds) {

    // Verif :
    assert (bounds.get_ndim()== 1) ;

    ndim = 2 ;

    n_shells = bounds.get_size(0) - 3;
    ADAPTED_OUTER = 1;
    ADAPTED_INNER = 2;

    nbr_domains = bounds.get_size(0)+1 ;
    assert (nbr_domains>=4) ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;

    // Nucleus
    domains[0] = new Domain_polar_nucleus (0, ttype, bounds(0), center, res) ;
    domains[ADAPTED_OUTER] = new Domain_polar_shell_outer_adapted (*this, ADAPTED_OUTER, ttype, bounds(ADAPTED_OUTER-1), bounds(ADAPTED_OUTER), center, res) ;
    domains[ADAPTED_INNER] = new Domain_polar_shell_inner_adapted (*this, ADAPTED_INNER, ttype, bounds(ADAPTED_OUTER), bounds(ADAPTED_INNER), center, res) ;
    for (int i=3 ; i<nbr_domains-1 ; i++)
       domains[i] = new Domain_polar_shell(i, ttype, bounds(i-1), bounds(i), center, res) ;
   domains[nbr_domains-1] = new Domain_polar_compact (nbr_domains-1, ttype, bounds(nbr_domains-2), center, res) ;

   const Domain_polar_shell_outer_adapted* pshell_outer = dynamic_cast<const Domain_polar_shell_outer_adapted*> (domains[ADAPTED_OUTER]) ;
   pshell_outer->vars_to_terms() ;
   const Domain_polar_shell_inner_adapted* pshell_inner = dynamic_cast<const Domain_polar_shell_inner_adapted*> (domains[ADAPTED_INNER]) ;
   pshell_inner->vars_to_terms() ;
}


Space_polar_adapted::Space_polar_adapted (FILE* fd) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;

  ADAPTED_OUTER = 1;
  ADAPTED_INNER = 2;

	domains = new Domain* [nbr_domains] ;
	// interior
	domains[0] = new Domain_polar_nucleus (0, fd) ;
	domains[ADAPTED_OUTER] = new Domain_polar_shell_outer_adapted (*this, ADAPTED_OUTER, fd) ;	
	domains[ADAPTED_INNER] = new Domain_polar_shell_inner_adapted (*this, ADAPTED_INNER, fd) ;
  
  //Shells :
	for (int i=3 ; i<nbr_domains-1 ; i++)
		domains[i] = new Domain_polar_shell(i, fd) ;
	
	// Compactified
	domains[nbr_domains-1] = new Domain_polar_compact(nbr_domains-1, fd) ;  
	
	const Domain_polar_shell_outer_adapted* pshell_outer = dynamic_cast<const Domain_polar_shell_outer_adapted*> (domains[ADAPTED_OUTER]) ;
	pshell_outer->vars_to_terms() ;
	const Domain_polar_shell_inner_adapted* pshell_inner = dynamic_cast<const Domain_polar_shell_inner_adapted*> (domains[ADAPTED_INNER]) ;
	pshell_inner->vars_to_terms() ;
}

Space_polar_adapted::~Space_polar_adapted() {
  const Domain_polar_shell_outer_adapted* pshell_outer = dynamic_cast<const Domain_polar_shell_outer_adapted*> (domains[ADAPTED_OUTER]) ;
  pshell_outer->del_deriv() ;
  const Domain_polar_shell_inner_adapted* pshell_inner = dynamic_cast<const Domain_polar_shell_inner_adapted*> (domains[ADAPTED_INNER]) ;
  pshell_inner->del_deriv() ;
}

void Space_polar_adapted::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;	
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

int Space_polar_adapted::nbr_unknowns_from_variable_domains () const {
  
    return domains[ADAPTED_OUTER]->nbr_unknowns_from_adapted() ;
}

void Space_polar_adapted::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {
  
  bool found_outer = false ;
  int old_pos = pos ;
  domains[ADAPTED_OUTER]->affecte_coef(pos, cc, found_outer) ;
  pos = old_pos ;
  bool found_inner = false ;
  domains[ADAPTED_INNER]->affecte_coef(pos, cc, found_inner) ;
  assert (found_outer == found_inner) ;
  if (found_outer) {
    zedoms.set(0) = ADAPTED_OUTER ;
    zedoms.set(1) = ADAPTED_INNER ;
  }
  else {
    zedoms.set(0) = -1 ;
    zedoms.set(1) = -1 ;
  }
}

void Space_polar_adapted::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {
  
  int old_conte = conte ;
  domains[ADAPTED_OUTER]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[ADAPTED_INNER]->xx_to_ders_from_adapted(xx, conte) ;
}

void Space_polar_adapted::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {
  
  // First get the corrections :
  int old_pos = pos ;
  Val_domain cor_outer (domains[ADAPTED_OUTER]) ;
  domains[ADAPTED_OUTER]->xx_to_vars_from_adapted(cor_outer, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner (domains[ADAPTED_INNER]) ;
  domains[ADAPTED_INNER]->xx_to_vars_from_adapted(cor_inner, xx, pos) ;
      
  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) { 
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
	domains[ADAPTED_OUTER]->update_variable(cor_outer, *sys->var[i]->cmp[n], res) ;
	domains[ADAPTED_INNER]->update_variable(cor_inner, *sys->var[i]->cmp[n], res) ;
	
	sys->var[i]->cmp[n]->set_domain(ADAPTED_OUTER) = res(ADAPTED_OUTER) ;
	sys->var[i]->cmp[n]->set_domain(ADAPTED_INNER) = res(ADAPTED_INNER) ; 
	
    }
  }
   
  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++) 
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {
	
	Scalar so (*this) ;
	so = 0 ;
  so.set_domain(ADAPTED_OUTER) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED_OUTER-sys->dom_min]->val_t->cmp[n])(ADAPTED_OUTER) ;
  so.set_domain(ADAPTED_INNER) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED_INNER-sys->dom_min]->val_t->cmp[n])(ADAPTED_INNER) ;
	
	Scalar res(*this) ;
	res = 0 ;
	domains[ADAPTED_OUTER]->update_constante(cor_outer, so, res) ;
	domains[ADAPTED_INNER]->update_constante(cor_inner, so, res) ;
	
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED_OUTER-sys->dom_min]->val_t->cmp[n]->set_domain(ADAPTED_OUTER) = res(ADAPTED_OUTER) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED_INNER-sys->dom_min]->val_t->cmp[n]->set_domain(ADAPTED_INNER) = res(ADAPTED_INNER) ;
	
      }
    
      // Update the mapping :
      domains[ADAPTED_OUTER]->update_mapping(cor_outer) ;
      domains[ADAPTED_INNER]->update_mapping(cor_inner) ;
}

void Space_polar_adapted::add_eq_ori (System_of_eqs& sys, const char* name) {

	Index pos (domains[0]->get_nbr_points()) ;
	char auxi[LMAX] ;
	trim_spaces (auxi, name) ;
	sys.add_eq_val (0, auxi, pos) ;
}
}
