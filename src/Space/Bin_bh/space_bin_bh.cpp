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
#include "bin_bh.hpp"
#include "utilities.hpp"
#include "point.hpp"
#include "scalar.hpp"
#include "system_of_eqs.hpp"
#include "name_tools.hpp"
namespace Kadath {
double eta_lim_chi (double chi, double rext, double a, double eta_c) ;
double chi_lim_eta (double chi, double rext, double a, double chi_c) ;

double zerosec(double (*f)(double, const Param&), const Param& parf,
    double x1, double x2, double precis, int nitermax, int& niter) ;


inline double func_ab (double aa, const Param& par) {
	double r1 = par.get_double(0) ;
	double r2 = par.get_double(1) ;
	double d = par.get_double(2) ;
	return (sqrt(aa*aa+r1*r1)+sqrt(aa*aa+r2*r2)-d) ;
}

Space_bin_bh::Space_bin_bh (int ttype, double dist, double rbh1, double rbh2, double rext, int nr) {

    ndim = 3;
    nbr_domains = 12;

    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;

    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;

     // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(2*rbh1,0) ;
    par_a.add_double(2*rbh2,1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_ab, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/2/rbh2) ;
    double eta_minus = -asinh(aa/2/rbh1) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;

    // First BH
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus (0, ttype, rbh1/4., center_minus, res) ;
    domains[1] = new Domain_shell_outer_homothetic (*this, 1, ttype, rbh1/4., rbh1, center_minus, res) ;
    domains[2] = new Domain_shell_inner_homothetic (*this, 2, ttype, rbh1, rbh1*2., center_minus, res) ;

    // Second bh
    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[3] = new Domain_nucleus (3, ttype, rbh2/4., center_plus, res) ;
    domains[4] = new Domain_shell_outer_homothetic (*this, 4, ttype, rbh2/4., rbh2, center_plus, res) ;
    domains[5] = new Domain_shell_inner_homothetic (*this, 5, ttype, rbh2, rbh2*2., center_plus, res) ;

    // Bispheric part
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[7] = new Domain_bispheric_rect(7, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[8] = new Domain_bispheric_eta_first(8, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[9] = new Domain_bispheric_rect(9, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[10] = new Domain_bispheric_chi_first(10, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;

    // Compactified domain
    Point center(3) ;
    domains[11] = new Domain_compact (11, ttype, rext, center, res) ;

    const Domain_shell_outer_homothetic* pouter_1 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_homothetic* pinner_1 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[2]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_homothetic* pouter_2 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[4]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_homothetic* pinner_2 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[5]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}

Space_bin_bh::Space_bin_bh (Space_bin_bh const & sp) {

	nbr_domains = sp.nbr_domains;
  ndim = sp.ndim;
  type_base = sp.type_base;

  n_shells1 = sp.n_shells1;
  n_shells2 = sp.n_shells2;

  // We calculate domain indicies and ensure they match
  // the imported space - should be consistent!
  this->BH1 = 0;
  assert(BH1 == sp.BH1);

  this->BH2 = 3 + n_shells1;
  assert(BH2 == sp.BH2);

  this->OUTER = 6 + n_shells1 + n_shells2;
  assert(OUTER == sp.OUTER);

  this->n_shells_outer = nbr_domains - 12 - n_shells1 - n_shells2;
  assert(n_shells_outer >= 0 && n_shells_outer == sp.n_shells_outer);

	domains = new Domain* [nbr_domains] ;
	//First BH :
  const Domain_nucleus* d_nuc1 = dynamic_cast<const Domain_nucleus*> (sp.get_domain(BH1)) ;
  domains[BH1] = new Domain_nucleus(*d_nuc1, true);

  const Domain_shell_outer_homothetic* sp_pouter_1 = dynamic_cast<const Domain_shell_outer_homothetic*> (sp.get_domain(BH1+1)) ;
  domains[BH1+1] = new Domain_shell_outer_homothetic(*this, *sp_pouter_1) ;
  
  const Domain_shell_inner_homothetic* sp_pinner_1 = dynamic_cast<const Domain_shell_inner_homothetic*> (sp.get_domain(BH1+2)) ;
  domains[BH1+2] = new Domain_shell_inner_homothetic(*this, *sp_pinner_1) ;

  for(int i = 0; i < n_shells1; ++i) {
    const Domain_shell* d_shell = dynamic_cast<const Domain_shell*> (sp.get_domain(BH1+3+i)) ;
    domains[BH1+3+i] = new Domain_shell(*d_shell, true);
  }    

	//second BH :
  const Domain_nucleus* d_nuc2 = dynamic_cast<const Domain_nucleus*> (sp.get_domain(BH2)) ;
  domains[BH2] = new Domain_nucleus(*d_nuc2, true);

  const Domain_shell_outer_homothetic* sp_pouter_2 = dynamic_cast<const Domain_shell_outer_homothetic*> (sp.get_domain(BH2+1)) ;
  domains[BH2+1] = new Domain_shell_outer_homothetic(*this, *sp_pouter_2) ;
  
  const Domain_shell_inner_homothetic* sp_pinner_2 = dynamic_cast<const Domain_shell_inner_homothetic*> (sp.get_domain(BH2+2)) ;
  domains[BH2+2] = new Domain_shell_inner_homothetic(*this, *sp_pinner_2) ;

  for(int i = 0; i < n_shells1; ++i) {
    const Domain_shell* d_shell = dynamic_cast<const Domain_shell*> (sp.get_domain(BH2+3+i)) ;
    domains[BH2+3+i] = new Domain_shell(*d_shell, true);
  } 

	// Bispheric
  const Domain_bispheric_chi_first* sp_chi_1 = dynamic_cast<const Domain_bispheric_chi_first*> (sp.get_domain(OUTER)) ;
	domains[OUTER] = new Domain_bispheric_chi_first(*this, *sp_chi_1) ;

  const Domain_bispheric_rect* sp_rec_1 = dynamic_cast<const Domain_bispheric_rect*> (sp.get_domain(OUTER+1)) ;
	domains[OUTER+1] = new Domain_bispheric_rect(*this, *sp_rec_1) ;

  const Domain_bispheric_eta_first* sp_eta = dynamic_cast<const Domain_bispheric_eta_first*> (sp.get_domain(OUTER+2)) ;
  domains[OUTER+2] = new Domain_bispheric_eta_first(*this, *sp_eta) ;

  const Domain_bispheric_rect* sp_rec_2 = dynamic_cast<const Domain_bispheric_rect*> (sp.get_domain(OUTER+3)) ;
	domains[OUTER+3] = new Domain_bispheric_rect(*this, *sp_rec_2) ;

  const Domain_bispheric_chi_first* sp_chi_2 = dynamic_cast<const Domain_bispheric_chi_first*> (sp.get_domain(OUTER+4)) ;
	domains[OUTER+4] = new Domain_bispheric_chi_first(*this, *sp_chi_2) ;

	for (int i=0 ; i<n_shells_outer ; i++) {
		const Domain_shell* d_shell = dynamic_cast<const Domain_shell*> (sp.get_domain(OUTER+5+i)) ;
    domains[OUTER+5+i] = new Domain_shell(*d_shell, true);
  }

	// Compactified
  const Domain_compact* d_compact = dynamic_cast<const Domain_compact*> (sp.get_domain(nbr_domains-1)) ;
	domains[nbr_domains-1] = new Domain_compact(*d_compact, true) ;

  const Domain_shell_outer_homothetic* pouter_1 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH1+1]) ;
  pouter_1->vars_to_terms() ;
  pouter_1->update() ;
  const Domain_shell_inner_homothetic* pinner_1 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH1+2]) ;
  pinner_1->vars_to_terms() ;
  pinner_1->update() ;
  const Domain_shell_outer_homothetic* pouter_2 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH2+1]) ;
  pouter_2->vars_to_terms() ;
  pouter_2->update() ;
  const Domain_shell_inner_homothetic* pinner_2 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH2+2]) ;
  pinner_2->vars_to_terms() ;
  pinner_2->update() ;
}

Space_bin_bh::Space_bin_bh (FILE* fd) {

	fread_be (&nbr_domains, sizeof(int), 1, fd) ;

	fread_be (&n_shells1, sizeof(int), 1, fd) ;
	fread_be (&n_shells2, sizeof(int), 1, fd) ;

	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;

  this->BH1 = 0;
  this->BH2 = 3 + n_shells1;
  this->OUTER = 6 + n_shells1 + n_shells2;
  this->n_shells_outer = nbr_domains - 12 - n_shells1 - n_shells2;
  assert(n_shells_outer >= 0);

	domains = new Domain* [nbr_domains] ;
	//First BH :
	domains[BH1] = new Domain_nucleus (BH1, fd) ;
	domains[BH1+1] = new Domain_shell_outer_homothetic (*this, BH1+1, fd) ;
	domains[BH1+2] = new Domain_shell_inner_homothetic (*this, BH1+2, fd) ;

  for(int i = 0; i < n_shells1; ++i)
    domains[BH1+3+i] = new Domain_shell(BH1+3+i, fd);

	//second BH :
	domains[BH2] = new Domain_nucleus (BH2, fd) ;
	domains[BH2+1] = new Domain_shell_outer_homothetic (*this, BH2+1, fd) ;
	domains[BH2+2] = new Domain_shell_inner_homothetic (*this, BH2+2, fd) ;

  for(int i = 0; i < n_shells2; ++i)
    domains[BH2+3+i] = new Domain_shell(BH2+3+i, fd);

	// Bispheric
	domains[OUTER] = new Domain_bispheric_chi_first(OUTER, fd) ;
  domains[OUTER+1] = new Domain_bispheric_rect(OUTER+1, fd) ;
  domains[OUTER+2] = new Domain_bispheric_eta_first(OUTER+2, fd) ;
  domains[OUTER+3] = new Domain_bispheric_rect(OUTER+3, fd) ;
  domains[OUTER+4] = new Domain_bispheric_chi_first(OUTER+4, fd) ;

	for (int i=0 ; i<n_shells_outer ; i++)
		domains[OUTER+5+i] = new Domain_shell (OUTER+5+i, fd) ;

	// Compactified
	domains[OUTER+5+n_shells_outer] = new Domain_compact(OUTER+5+n_shells_outer, fd) ;

  const Domain_shell_outer_homothetic* pouter_1 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH1+1]) ;
  pouter_1->vars_to_terms() ;
  pouter_1->update() ;
  const Domain_shell_inner_homothetic* pinner_1 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH1+2]) ;
  pinner_1->vars_to_terms() ;
  pinner_1->update() ;
  const Domain_shell_outer_homothetic* pouter_2 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH2+1]) ;
  pouter_2->vars_to_terms() ;
  pouter_2->update() ;
  const Domain_shell_inner_homothetic* pinner_2 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH2+2]) ;
  pinner_2->vars_to_terms() ;
  pinner_2->update() ;
}

Space_bin_bh::~Space_bin_bh() {
  const Domain_shell_outer_homothetic* pouter_1 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH1+1]) ;
  pouter_1->del_deriv() ;
  const Domain_shell_inner_homothetic* pinner_1 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH1+2]) ;
  pinner_1->del_deriv() ;
  const Domain_shell_outer_homothetic* pouter_2 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH2+1]) ;
  pouter_2->del_deriv() ;
  const Domain_shell_inner_homothetic* pinner_2 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH2+2]) ;
  pinner_2->del_deriv() ;
}

void Space_bin_bh::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;

	fwrite_be (&n_shells1, sizeof(int), 1, fd) ;
	fwrite_be (&n_shells2, sizeof(int), 1, fd) ;

	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

int Space_bin_bh::nbr_unknowns_from_variable_domains () const {
  return (domains[BH1+2]->nbr_unknowns_from_adapted() + domains[BH2+2]->nbr_unknowns_from_adapted()) ;
}

void Space_bin_bh::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {
   // In BH 1 ?
  bool found_outer_1 = false ;
  int old_pos = pos ;
  domains[BH1+1]->affecte_coef(pos, cc, found_outer_1) ;
  pos = old_pos ;
  bool found_inner_1 = false ;
  domains[BH1+2]->affecte_coef(pos, cc, found_inner_1) ;
  assert (found_outer_1 == found_inner_1) ;
  if (found_outer_1) {
    zedoms.set(0) = BH1+1 ;
    zedoms.set(1) = BH1+2 ;
  }
  else {
    zedoms.set(0) = -1 ;
    zedoms.set(1) = -1 ;
  }

  // In BH 2 ?
  bool found_outer_2 = false ;
  old_pos = pos ;
  domains[BH2+1]->affecte_coef(pos, cc, found_outer_2) ;
  pos = old_pos ;
  bool found_inner_2 = false ;
  domains[BH2+2]->affecte_coef(pos, cc, found_inner_2) ;
  assert (found_outer_2 == found_inner_2) ;
  if (found_outer_2) {
    zedoms.set(0) = BH2+1 ;
    zedoms.set(1) = BH2+2 ;
  }
}

void Space_bin_bh::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {
  // BH 1
  int old_conte = conte ;
  domains[BH1+1]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[BH1+2]->xx_to_ders_from_adapted(xx, conte) ;
  // BH 2
  old_conte = conte ;
  domains[BH2+1]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[BH2+2]->xx_to_ders_from_adapted(xx, conte) ;
}

void Space_bin_bh::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {
  // First get the corrections :
  // BH 1
  int old_pos = pos ;
  Val_domain cor_outer_1 (domains[BH1+1]) ;
  domains[BH1+1]->xx_to_vars_from_adapted(cor_outer_1, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner_1 (domains[BH1+2]) ;
  domains[BH1+2]->xx_to_vars_from_adapted(cor_inner_1, xx, pos) ;

  // BH 2
  old_pos = pos ;
  Val_domain cor_outer_2 (domains[BH2+1]) ;
  domains[BH2+1]->xx_to_vars_from_adapted(cor_outer_2, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner_2 (domains[BH2+2]) ;
  domains[BH2+2]->xx_to_vars_from_adapted(cor_inner_2, xx, pos) ;

  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) {
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
	domains[BH1+1]->update_variable(cor_outer_1, *sys->var[i]->cmp[n], res) ;
	domains[BH1+2]->update_variable(cor_inner_1, *sys->var[i]->cmp[n], res) ;
	domains[BH2+1]->update_variable(cor_outer_2, *sys->var[i]->cmp[n], res) ;
	domains[BH2+2]->update_variable(cor_inner_2, *sys->var[i]->cmp[n], res) ;


	sys->var[i]->cmp[n]->set_domain(BH1+1) = res(BH1+1) ;
	sys->var[i]->cmp[n]->set_domain(BH1+2) = res(BH1+2) ;
	sys->var[i]->cmp[n]->set_domain(BH2+1) = res(BH2+1) ;
	sys->var[i]->cmp[n]->set_domain(BH2+2) = res(BH2+2) ;
    }
  }

  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++)
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {

	Scalar so (*this) ;
	so = 0 ;
  so.set_domain(BH1+1) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + BH1+1]->val_t->cmp[n])(BH1+1) ;
  so.set_domain(BH1+2) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + BH1+2]->val_t->cmp[n])(BH1+2) ;
  so.set_domain(BH2+1) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + BH2+1]->val_t->cmp[n])(BH2+1) ;
  so.set_domain(BH2+2) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + BH2+2]->val_t->cmp[n])(BH2+2) ;

	Scalar res(*this) ;
	res = 0 ;
	domains[BH1+1]->update_constante(cor_outer_1, so, res) ;
	domains[BH1+2]->update_constante(cor_inner_1, so, res) ;
	domains[BH2+1]->update_constante(cor_outer_2, so, res) ;
	domains[BH2+2]->update_constante(cor_inner_2, so, res) ;

	sys->cst[i*(sys->dom_max-sys->dom_min+1)+BH1+1]->val_t->cmp[n]->set_domain(BH1+1) = res(BH1+1) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+BH1+2]->val_t->cmp[n]->set_domain(BH1+2) = res(BH1+2) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+BH2+1]->val_t->cmp[n]->set_domain(BH2+1) = res(BH2+1) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+BH2+2]->val_t->cmp[n]->set_domain(BH2+2) = res(BH2+2) ;
      }

  // Update the mapping :
  domains[BH1+1]->update_mapping(cor_outer_1) ;
  domains[BH1+2]->update_mapping(cor_inner_1) ;
  domains[BH2+1]->update_mapping(cor_outer_2) ;
  domains[BH2+2]->update_mapping(cor_inner_2) ;

}


Array<int> Space_bin_bh::get_indices_matching_non_std(int dom, int bound) const {
	if (dom==BH1 + n_shells1 + 2) {
	  // First BH ;
	  Array<int> res (2,2) ;
	  switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = OUTER ; // Matching with chi first
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = OUTER+1 ; // Matching with rect
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;
		  std::_Exit(EXIT_FAILURE) ;
	  }
	return res ;
	}

	if (dom== BH2 + n_shells2 + 2) {
		// Second BH ;
		Array<int> res(2, 2) ;
		switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = OUTER+3 ; // Matching with rect
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = OUTER+4 ; // Matching with chi_first
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;			std::_Exit(EXIT_FAILURE) ;
		}
	return res ;
	}

	if (dom == OUTER)	{
	  // first chi first :
	  Array<int> res(2,1) ;
	  switch (bound) {
	    case INNER_BC :
	      res.set(0,0) = BH1 + n_shells1 + 2 ; // First bh
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;			std::_Exit(EXIT_FAILURE) ;
	}
	return res ;
	}

	if (dom == OUTER+1)	{
	  // first rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC :
	      res.set(0,0) = BH1 + n_shells1 + 2 ; // First bh
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;			std::_Exit(EXIT_FAILURE) ;
	}
	return res ;
	}


	if (dom == OUTER+2) {
	  // eta first
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case OUTER_BC :
	      res.set(0, 0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1, 0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;			std::_Exit(EXIT_FAILURE) ;
	}
	return res ;
	}

	if (dom == OUTER+3) {
	  // second rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
		case INNER_BC :
		  res.set(0, 0) = BH2 + n_shells2 + 2 ; // Second bh
		  res.set(1, 0) = OUTER_BC ;
		  break ;
		case OUTER_BC :
		  res.set(0, 0) = OUTER+5 ; // Compactified domain or first shell
		  res.set(1, 0) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;			std::_Exit(EXIT_FAILURE) ;
	}
	return res ;
	}

	if (dom == OUTER+4) {
	  // second chi first :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC :
	      res.set(0,0) = BH2 + n_shells2 + 2 ; // second BH
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
		cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;
		std::_Exit(EXIT_FAILURE) ;
	  }
	return res ;
	}

	if (dom==OUTER+5) {
	    // compactified domain or first shell :
	    Array<int> res(2,5) ;
	    switch (bound) {
		case INNER_BC :
		  res.set(0,0) = OUTER ; // first chi first
		  res.set(0,1) = OUTER+1 ; // first rect
		  res.set(0,2) = OUTER+2 ; // eta first
		  res.set(0,3) = OUTER+3 ; // second rect
		  res.set(0,4) = OUTER+4 ; // second chi first
		  // Outer BC for all :
		  for (int i=0 ; i<5 ; i++)
		      res.set(1,i) = OUTER_BC ;
		  break ;
	  default :
	      cerr << "Bad bound in Space_bin_bh::get_indices_matching_non_std" << endl ;
	      std::_Exit(EXIT_FAILURE) ;
	}
	return res ;
	}

	cerr << "Bad domain in Space_bin_bh::get_indices_matching_non_std" << endl ;
	std::_Exit(EXIT_FAILURE) ;
}

Space_bin_bh::Space_bin_bh (int ttype, double dist, const std::vector<double>& BH1_bounds, 
                            const std::vector<double>& BH2_bounds, 
                            const std::vector<double>& outer_bounds, int nr) {

    ndim = 3;
    enum bound_index {RIN=0, RMID=1, NUM_IND};
    double rext = outer_bounds[0];
    
    n_shells1 = BH1_bounds.size()-3;
    n_shells2 = BH2_bounds.size()-3;
    n_shells_outer = outer_bounds.size()-1;

    BH2 = 3 + n_shells1;
    OUTER = 6 + n_shells1 + n_shells2;

    auto router1 = BH1_bounds[n_shells1+2];
    auto router2 = BH2_bounds[n_shells2+2];

    nbr_domains = 12 + n_shells1 + n_shells2 + n_shells_outer;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;

    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;

    // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(router1,0) ;
    par_a.add_double(router2,1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_ab, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/router2) ;
    double eta_minus = -asinh(aa/router1) ;

    double chi_c = 2.*atan(aa/rext) ;
    double eta_c = log((1.+rext/aa)/(rext/aa-1.)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;

    auto gen_bh_domains = [&](const int nuc_i, auto& bounds, const double& eta) {
        int shells = bounds.size()-3;
        Point center (ndim) ;
        center.set(1) = aa*cosh(eta)/sinh(eta) ;
        domains[nuc_i] = new Domain_nucleus (nuc_i, ttype, bounds[RIN], center, res) ;
        domains[nuc_i+1] = 
            new Domain_shell_outer_homothetic (*this, nuc_i+1, ttype, bounds[RIN], bounds[RMID], center, res) ;
        domains[nuc_i+2] = 
            new Domain_shell_inner_homothetic (*this, nuc_i+2, ttype, bounds[RMID], bounds[RMID+1], center, res) ;
        
        for(int i = 0; i < shells; ++i)
            domains[nuc_i+3+i] = new Domain_shell(nuc_i+3+i, ttype, bounds[RMID+1+i], bounds[RMID+1+i+1], center, res);
    };
    gen_bh_domains(BH1, BH1_bounds, eta_minus);
    gen_bh_domains(BH2, BH2_bounds, eta_plus);

    // Bispheric part
    domains[OUTER  ] = new Domain_bispheric_chi_first(OUTER, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[OUTER+1] = new Domain_bispheric_rect(OUTER+1, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[OUTER+2] = new Domain_bispheric_eta_first(OUTER+2, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[OUTER+3] = new Domain_bispheric_rect(OUTER+3, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[OUTER+4] = new Domain_bispheric_chi_first(OUTER+4, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;

    Point center(3) ;
  	for (int i=0 ; i<n_shells_outer ; i++)
	  	domains[OUTER+5+i] = new Domain_shell (OUTER+5+i, ttype, outer_bounds[i], outer_bounds[i+1], center, res) ;

	  // Compactified
	  domains[OUTER+5+n_shells_outer] = 
        new Domain_compact(OUTER+5+n_shells_outer, ttype, outer_bounds[n_shells_outer], center, res) ;

    const Domain_shell_outer_homothetic* pouter_1 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH1+1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_homothetic* pinner_1 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH1+2]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_homothetic* pouter_2 = dynamic_cast<const Domain_shell_outer_homothetic*> (domains[BH2+1]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_homothetic* pinner_2 = dynamic_cast<const Domain_shell_inner_homothetic*> (domains[BH2+2]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}
}
