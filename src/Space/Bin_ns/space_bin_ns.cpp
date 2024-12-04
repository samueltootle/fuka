
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
double eta_lim_chi (double chi, double rext, double a, double eta_c) ;
double chi_lim_eta (double chi, double rext, double a, double chi_c) ;

double zerosec(double (*f)(double, const Param&), const Param& parf,
    double x1, double x2, double precis, int nitermax, int& niter) ;


double func_abns (double aa, const Param& par) {
	double r1 = par.get_double(0) ;
	double r2 = par.get_double(1) ;
	double d = par.get_double(2) ;
	return (sqrt(aa*aa+r1*r1)+sqrt(aa*aa+r2*r2)-d) ;
}

Space_bin_ns::Space_bin_ns (int ttype, double dist, double rinstar1, double rstar1, double routstar1,
                            double rinstar2, double rstar2, double routstar2, double rext, int nr, int nshells1, int nshells2) {

    ndim = 3 ;

    n_shells1 = nshells1;
    n_shells2 = nshells2;

    NS1 = 0;
    NS2 = 3 + n_shells1;
    ADAPTED1 = NS1 + n_shells1 + 1;
    ADAPTED2 = NS2 + n_shells2 + 1;
    OUTER = 6 + n_shells1 + n_shells2;

    n_shells_outer = 0;
    nbr_domains = 12 + n_shells1 + n_shells2;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;

    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;

     // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(routstar1,0) ;
    par_a.add_double(routstar2,1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_abns, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/routstar2) ;
    double eta_minus = -asinh(aa/routstar1) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;
    // First NS
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[NS1] = new Domain_nucleus (NS1, ttype, rinstar1, center_minus, res) ;

    double delta_r1 = (rstar1 - rinstar1);
    std::vector<double> r_shells1(n_shells1+1);
    //r_shells1[0] = rinstar1;
    for(int i = 0; i <= n_shells1; ++i) {
      //r_shells1[i] = (delta_r1/(M_PI/2))*atan(i+1)+rstar1;  //puts shells closer to router
      r_shells1[i] = rinstar1 + delta_r1 / (nshells1 + 1) * i;
      //std::cout << r_shells1[i] << std::endl;
    }

    for(int i = 0; i < n_shells1; ++i) {
      domains[NS1+1+i] = new Domain_shell(NS1+1+i, ttype, r_shells1[i], r_shells1[i+1], center_minus, res);
    }

    domains[ADAPTED1] = new Domain_shell_outer_adapted (*this, ADAPTED1, ttype, r_shells1[n_shells1], rstar1, center_minus, res) ;
    domains[ADAPTED1+1] = new Domain_shell_inner_adapted (*this, ADAPTED1+1, ttype, rstar1, routstar1, center_minus, res) ;

    // Second NS
    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[NS2] = new Domain_nucleus (NS2, ttype, rinstar2, center_plus, res) ;

    double delta_r2 = (rstar2 - rinstar2);
    std::vector<double> r_shells2(n_shells2+1);
    //r_shells2[0] = rinstar2;
    for(int i = 0; i <= n_shells2; ++i) {
      //r_shells2[i] = (delta_r1/(M_PI/2))*atan(i+1)+rstar2;  //puts shells closer to router
      r_shells2[i] = rinstar2 + delta_r2 / (nshells2 + 1) * i;
      //std::cout << r_shells2[i] << std::endl;
    }

    for(int i = 0; i < n_shells2; ++i) {
      domains[NS2+1+i] = new Domain_shell(NS2+1+i, ttype, r_shells2[i], r_shells2[i+1], center_plus, res);
    }

    domains[ADAPTED2] = new Domain_shell_outer_adapted (*this, ADAPTED2, ttype, r_shells2[n_shells2], rstar2, center_plus, res) ;
    domains[ADAPTED2+1] = new Domain_shell_inner_adapted (*this, ADAPTED2+1, ttype, rstar2, routstar2, center_plus, res) ;

    // Bispheric part
    domains[OUTER] = new Domain_bispheric_chi_first(OUTER, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[OUTER+1] = new Domain_bispheric_rect(OUTER+1, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[OUTER+2] = new Domain_bispheric_eta_first(OUTER+2, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[OUTER+3] = new Domain_bispheric_rect(OUTER+3, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[OUTER+4] = new Domain_bispheric_chi_first(OUTER+4, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;

    // Compactified domain
    Point center(3) ;
    domains[OUTER+5] = new Domain_compact (OUTER+5, ttype, rext, center, res) ;

    const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED1+1]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED2]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED2+1]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}

Space_bin_ns::Space_bin_ns (int ttype, double dist, double rinstar1, double rstar1, double routstar1,
			double rinstar2, double rstar2, double routstar2, double rext, double rshell, int nr) {

    ndim = 3 ;

    n_shells1 = 0;
    n_shells2 = 0;

    NS1 = 0;
    NS2 = 3;
    ADAPTED1 = NS1 + n_shells1 + 1;
    ADAPTED2 = NS2 + n_shells2 + 1;
    OUTER = 6;

    n_shells_outer = 1 ;
    nbr_domains = 13 ;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;

    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;

     // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(routstar1,0) ;
    par_a.add_double(routstar2,1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_abns, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus = asinh(aa/routstar2) ;
    double eta_minus = -asinh(aa/routstar1) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;

    // First NS
    Point center_minus (ndim) ;
    center_minus.set(1) = aa*cosh(eta_minus)/sinh(eta_minus) ;
    domains[0] = new Domain_nucleus (0, ttype, rinstar1, center_minus, res) ;
    domains[1] = new Domain_shell_outer_adapted (*this, 1, ttype, rinstar1, rstar1, center_minus, res) ;
    domains[2] = new Domain_shell_inner_adapted (*this, 2, ttype, rstar1, routstar1, center_minus, res) ;

    // Second NS
    Point center_plus (ndim) ;
    center_plus.set(1) = aa*cosh(eta_plus)/sinh(eta_plus) ;
    domains[3] = new Domain_nucleus (3, ttype, rinstar2, center_plus, res) ;
    domains[4] = new Domain_shell_outer_adapted (*this, 4, ttype, rinstar2, rstar2, center_plus, res) ;
    domains[5] = new Domain_shell_inner_adapted (*this, 5, ttype, rstar2, routstar2, center_plus, res) ;

    // Bispheric part
    domains[6] = new Domain_bispheric_chi_first(6, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[7] = new Domain_bispheric_rect(7, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[8] = new Domain_bispheric_eta_first(8, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[9] = new Domain_bispheric_rect(9, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[10] = new Domain_bispheric_chi_first(10, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;

    // Shell in 1/R
    Point center(3) ;
    domains[11] = new Domain_shell (11, ttype, rext, rshell,  center, res) ;

    // Compactified domain
    domains[12] = new Domain_compact (12, ttype, rshell,  center, res) ;

    const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[2]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[4]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[5]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}

Space_bin_ns::Space_bin_ns (Space_bin_ns const & sp) {

	nbr_domains = sp.nbr_domains;
  ndim = sp.ndim;
  type_base = sp.type_base;

  n_shells1 = sp.n_shells1;
  n_shells2 = sp.n_shells2;

  n_inner_shells1 = sp.n_inner_shells1;
  n_inner_shells2 = sp.n_inner_shells2;

  // We calculate domain indicies and ensure they match
  // the imported space - should be consistent!
  this->NS1 = 0;
  assert(NS1 == sp.NS1);

  this->NS2 = 3 + n_inner_shells1 + n_shells1;
  assert(NS2 == sp.NS2);

  this->ADAPTED1 = NS1 + n_inner_shells1 + 1;
  assert(ADAPTED1 == sp.ADAPTED1);

  this->ADAPTED2 = NS2 + n_inner_shells2 + 1;
  assert(ADAPTED2 == sp.ADAPTED2);

  this->OUTER = 6 + n_shells1 + n_inner_shells1 + n_shells2 + n_inner_shells2;
  assert(OUTER == sp.OUTER);

  this->n_shells_outer = nbr_domains - 12 - n_inner_shells1 - n_shells1 - n_inner_shells2 - n_shells2;
  assert(n_shells_outer >= 0 && n_shells_outer == sp.n_shells_outer);

	domains = new Domain* [nbr_domains] ;

  auto add_spherical_shells = [&](auto start_idx, auto nshells) {
    for(int i = 0; i < nshells; ++i) {
      const Domain_shell* d_shell = dynamic_cast<const Domain_shell*> (sp.get_domain(start_idx+i)) ;
      domains[start_idx+i] = new Domain_shell(*d_shell, true);
    }
  };

	//First NS :
  const Domain_nucleus* d_nuc1 = dynamic_cast<const Domain_nucleus*> (sp.get_domain(NS1)) ;
  domains[NS1] = new Domain_nucleus(*d_nuc1, true);

  add_spherical_shells(NS1+1, n_inner_shells1);

  const Domain_shell_outer_adapted* sp_pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (sp.get_domain(ADAPTED1)) ;
  domains[ADAPTED1] = new Domain_shell_outer_adapted(*this, *sp_pouter_1) ;

  const Domain_shell_inner_adapted* sp_pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (sp.get_domain(ADAPTED1+1)) ;
  domains[ADAPTED1+1] = new Domain_shell_inner_adapted(*this, *sp_pinner_1) ;

  add_spherical_shells(ADAPTED1+2, n_shells1);

	//second NS :
  const Domain_nucleus* d_nuc2 = dynamic_cast<const Domain_nucleus*> (sp.get_domain(NS2)) ;
  domains[NS2] = new Domain_nucleus(*d_nuc2, true);

  add_spherical_shells(NS2+1, n_inner_shells2);

  const Domain_shell_outer_adapted* sp_pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (sp.get_domain(ADAPTED2)) ;
  domains[ADAPTED2] = new Domain_shell_outer_adapted(*this, *sp_pouter_2) ;

  const Domain_shell_inner_adapted* sp_pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (sp.get_domain(ADAPTED2+1)) ;
  domains[ADAPTED2+1] = new Domain_shell_inner_adapted(*this, *sp_pinner_2) ;

  add_spherical_shells(ADAPTED2+2, n_shells2);

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

  add_spherical_shells(OUTER+5, n_shells_outer);

	// Compactified
  const Domain_compact* d_compact = dynamic_cast<const Domain_compact*> (sp.get_domain(nbr_domains-1)) ;
	domains[nbr_domains-1] = new Domain_compact(*d_compact, true) ;

  const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED1]) ;
  pouter_1->vars_to_terms() ;
  pouter_1->update() ;
  const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED1+1]) ;
  pinner_1->vars_to_terms() ;
  pinner_1->update() ;
  const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED2]) ;
  pouter_2->vars_to_terms() ;
  pouter_2->update() ;
  const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED2+1]) ;
  pinner_2->vars_to_terms() ;
  pinner_2->update() ;
}

Space_bin_ns::Space_bin_ns (FILE* fd, bool old) {
	fread_be (&nbr_domains, sizeof(int), 1, fd) ;
  fread_be (&n_inner_shells1, sizeof(int), 1, fd) ;
  fread_be (&n_inner_shells2, sizeof(int), 1, fd) ;

  // in the original BNS, n_shells referred to interior shells
  // whereas it now refers to exterior shells like for the BBH and BHNS
  if(!old) {
    fread_be (&n_shells1, sizeof(int), 1, fd) ;
    fread_be (&n_shells2, sizeof(int), 1, fd) ;
  }
	fread_be (&ndim, sizeof(int), 1, fd) ;
	fread_be (&type_base, sizeof(int), 1, fd) ;

  n_shells_outer = nbr_domains - 12 - n_inner_shells1 - n_shells1 - n_shells2 - n_inner_shells2;
  assert(n_shells_outer >= 0);

  NS1 = 0;
  NS2 = 3 + n_inner_shells1 + n_shells1;
  ADAPTED1 = NS1 + n_inner_shells1 + 1;
  ADAPTED2 = NS2 + n_inner_shells2 + 1;
  OUTER = 6 + n_shells1 + n_inner_shells1 + n_shells2 + n_inner_shells2;

	domains = new Domain* [nbr_domains] ;

	//First NS :
	domains[NS1] = new Domain_nucleus (NS1, fd) ;

	for(int i = 0; i < n_inner_shells1; ++i)
	  domains[NS1+1+i]  = new Domain_shell(NS1+1+i, fd);

	domains[ADAPTED1]   = new Domain_shell_outer_adapted (*this, ADAPTED1, fd) ;
	domains[ADAPTED1+1] = new Domain_shell_inner_adapted (*this, ADAPTED1+1, fd) ;

  for(int i = 0; i < n_shells1; ++i)
	  domains[ADAPTED1+2+i]  = new Domain_shell(ADAPTED1+2+i, fd);

	//second NS :
	domains[NS2] = new Domain_nucleus (NS2, fd) ;

	for(int i = 0; i < n_inner_shells2; ++i)
	  domains[NS2+1+i]  = new Domain_shell(NS2+1+i, fd);

	domains[ADAPTED2]   = new Domain_shell_outer_adapted (*this, ADAPTED2, fd) ;
	domains[ADAPTED2+1] = new Domain_shell_inner_adapted (*this, ADAPTED2+1, fd) ;

  for(int i = 0; i < n_shells2; ++i)
	  domains[ADAPTED2+2+i]  = new Domain_shell(ADAPTED2+2+i, fd);

	// Bispheric
  domains[OUTER]   = new Domain_bispheric_chi_first(OUTER, fd) ;
  domains[OUTER+1] = new Domain_bispheric_rect(OUTER+1, fd) ;
  domains[OUTER+2] = new Domain_bispheric_eta_first(OUTER+2, fd) ;
  domains[OUTER+3] = new Domain_bispheric_rect(OUTER+3, fd) ;
  domains[OUTER+4] = new Domain_bispheric_chi_first(OUTER+4, fd) ;

	for (int i=0 ; i<n_shells_outer ; i++)
		domains[OUTER+5+i] = new Domain_shell (OUTER+5+i, fd) ;

	// Compactified
	domains[OUTER+5+n_shells_outer] = new Domain_compact(OUTER+5+n_shells_outer, fd) ;

  const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED1]) ;
  pouter_1->vars_to_terms() ;
  pouter_1->update() ;
  const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED1+1]) ;
  pinner_1->vars_to_terms() ;
  pinner_1->update() ;
  const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED2]) ;
  pouter_2->vars_to_terms() ;
  pouter_2->update() ;
  const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED2+1]) ;
  pinner_2->vars_to_terms() ;
  pinner_2->update() ;

}

Space_bin_ns::~Space_bin_ns() {
  const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED1]) ;
  pouter_1->del_deriv() ;
  const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED1+1]) ;
  pinner_1->del_deriv() ;
  const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED2]) ;
  pouter_2->del_deriv() ;
  const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED2+1]) ;
  pinner_2->del_deriv() ;
}

void Space_bin_ns::save (FILE* fd) const  {
	fwrite_be (&nbr_domains, sizeof(int), 1, fd) ;
  fwrite_be (&n_inner_shells1, sizeof(int), 1, fd) ;
  fwrite_be (&n_inner_shells2, sizeof(int), 1, fd) ;
	fwrite_be (&n_shells1, sizeof(int), 1, fd) ;
	fwrite_be (&n_shells2, sizeof(int), 1, fd) ;
	fwrite_be (&ndim, sizeof(int), 1, fd) ;
	fwrite_be (&type_base, sizeof(int), 1, fd) ;
	for (int i=0 ; i<nbr_domains ; i++)
		domains[i]->save(fd) ;
}

int Space_bin_ns::nbr_unknowns_from_variable_domains () const {

    return (domains[ADAPTED1]->nbr_unknowns_from_adapted() + domains[ADAPTED2]->nbr_unknowns_from_adapted()) ;
}

void Space_bin_ns::affecte_coef_to_variable_domains (int& pos, int cc, Array<int>& zedoms) const {

  // In star 1 ?
  bool found_outer_1 = false ;
  int old_pos = pos ;
  domains[ADAPTED1]->affecte_coef(pos, cc, found_outer_1) ;
  pos = old_pos ;
  bool found_inner_1 = false ;
  domains[ADAPTED1+1]->affecte_coef(pos, cc, found_inner_1) ;
  assert (found_outer_1 == found_inner_1) ;
  if (found_outer_1) {
    zedoms.set(0) = ADAPTED1 ;
    zedoms.set(1) = ADAPTED1+1 ;
  }
  else {
    zedoms.set(0) = -1 ;
    zedoms.set(1) = -1 ;
  }

  // In star 2 ?
    bool found_outer_2 = false ;
    old_pos = pos ;
    domains[ADAPTED2]->affecte_coef(pos, cc, found_outer_2) ;
    pos = old_pos ;
    bool found_inner_2 = false ;
    domains[ADAPTED2+1]->affecte_coef(pos, cc, found_inner_2) ;
    assert (found_outer_2 == found_inner_2) ;
    if (found_outer_2) {
      zedoms.set(0) = ADAPTED2 ;
      zedoms.set(1) = ADAPTED2+1 ;
    }
}

void Space_bin_ns::xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {

  // Star 1
  int old_conte = conte ;
  domains[ADAPTED1]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[ADAPTED1+1]->xx_to_ders_from_adapted(xx, conte) ;
  // Star 2
  old_conte = conte ;
  domains[ADAPTED2]->xx_to_ders_from_adapted(xx, conte) ;
  conte = old_conte ;
  domains[ADAPTED2+1]->xx_to_ders_from_adapted(xx, conte) ;
}

void Space_bin_ns::xx_to_vars_variable_domains (System_of_eqs* sys, const Array<double>& xx, int& pos) const  {


  // First get the corrections :
  // Star 1
  int old_pos = pos ;
  Val_domain cor_outer_1 (domains[ADAPTED1]) ;
  domains[ADAPTED1]->xx_to_vars_from_adapted(cor_outer_1, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner_1 (domains[ADAPTED1+1]) ;
  domains[ADAPTED1+1]->xx_to_vars_from_adapted(cor_inner_1, xx, pos) ;

  // Star 2
  old_pos = pos ;
  Val_domain cor_outer_2 (domains[ADAPTED2]) ;
  domains[ADAPTED2]->xx_to_vars_from_adapted(cor_outer_2, xx, pos) ;
  pos = old_pos ;
  Val_domain cor_inner_2 (domains[ADAPTED2+1]) ;
  domains[ADAPTED2+1]->xx_to_vars_from_adapted(cor_inner_2, xx, pos) ;

  // Now update the variables :
  for (int i=0 ; i<sys->nvar ; i++) {
    for (int n=0 ; n<sys->var[i]->get_n_comp() ; n++) {
	Scalar res(*this) ;
	domains[ADAPTED1]->update_variable(cor_outer_1, *sys->var[i]->cmp[n], res) ;
	domains[ADAPTED1+1]->update_variable(cor_inner_1, *sys->var[i]->cmp[n], res) ;
	domains[ADAPTED2]->update_variable(cor_outer_2, *sys->var[i]->cmp[n], res) ;
	domains[ADAPTED2+1]->update_variable(cor_inner_2, *sys->var[i]->cmp[n], res) ;


	sys->var[i]->cmp[n]->set_domain(ADAPTED1) = res(ADAPTED1) ;
	sys->var[i]->cmp[n]->set_domain(ADAPTED1+1) = res(ADAPTED1+1) ;
	sys->var[i]->cmp[n]->set_domain(ADAPTED2) = res(ADAPTED2) ;
	sys->var[i]->cmp[n]->set_domain(ADAPTED2+1) = res(ADAPTED2+1) ;
    }
  }

  // Now the constants :
  for (int i=0 ; i<sys->ncst ; i++)
    if (sys->cst[i*(sys->dom_max-sys->dom_min+1)]->get_type_data()==TERM_T)
      for (int n=0 ; n<sys->cst[i*(sys->dom_max-sys->dom_min+1)]->val_t->get_n_comp() ; n++) {

	Scalar so (*this) ;
	so = 0 ;
  so.set_domain(ADAPTED1  ) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + ADAPTED1  ]->val_t->cmp[n])(ADAPTED1  ) ;
  so.set_domain(ADAPTED1+1) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + ADAPTED1+1]->val_t->cmp[n])(ADAPTED1+1) ;
  so.set_domain(ADAPTED2  ) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + ADAPTED2  ]->val_t->cmp[n])(ADAPTED2  ) ;
  so.set_domain(ADAPTED2+1) = (*sys->cst[i*(sys->dom_max-sys->dom_min+1) + ADAPTED2+1]->val_t->cmp[n])(ADAPTED2+1) ;

	Scalar res(*this) ;
	res = 0 ;
	domains[ADAPTED1  ]->update_constante(cor_outer_1, so, res) ;
	domains[ADAPTED1+1]->update_constante(cor_inner_1, so, res) ;
	domains[ADAPTED2  ]->update_constante(cor_outer_2, so, res) ;
	domains[ADAPTED2+1]->update_constante(cor_inner_2, so, res) ;

	sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED1  ]->val_t->cmp[n]->set_domain(ADAPTED1  ) = res(ADAPTED1  ) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED1+1]->val_t->cmp[n]->set_domain(ADAPTED1+1) = res(ADAPTED1+1) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED2  ]->val_t->cmp[n]->set_domain(ADAPTED2  ) = res(ADAPTED2  ) ;
	sys->cst[i*(sys->dom_max-sys->dom_min+1)+ADAPTED2+1]->val_t->cmp[n]->set_domain(ADAPTED2+1) = res(ADAPTED2+1) ;
      }

      // Update the mapping :
      domains[ADAPTED1  ]->update_mapping(cor_outer_1) ;
      domains[ADAPTED1+1]->update_mapping(cor_inner_1) ;
      domains[ADAPTED2  ]->update_mapping(cor_outer_2) ;
      domains[ADAPTED2+1]->update_mapping(cor_inner_2) ;
}


Array<int> Space_bin_ns::get_indices_matching_non_std(int dom, int bound) const {
  // FIXME not entirely sure if this is correct with outer shells

	if (dom == ADAPTED1 + 1 + n_shells1) {
	  // First star ;
	  Array<int> res (2,2) ;
	  switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = OUTER ; // Matching with chi first
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = OUTER+1 ; // Matching with rect
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;
		  abort() ;
	  }
	return res ;
	}

	if (dom == ADAPTED2 + 1 + n_shells1) {
		// second star ;
		Array<int> res(2, 2) ;
		switch (bound) {
		case OUTER_BC :
		  res.set(0,0) = OUTER+3 ; // Matching with rect
		  res.set(1,0) = INNER_BC ;
		  res.set(0,1) = OUTER+4 ; // Matching with chi_first
		  res.set(1,1) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
		}
	return res ;
	}

	if (dom == OUTER)	{
	  // first chi first :
	  Array<int> res(2,1) ;
	  switch (bound) {
	    case INNER_BC :
	      res.set(0,0) = ADAPTED1+1+n_shells1; // First star
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	if (dom == OUTER+1)	{
	  // first rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC :
	      res.set(0,0) = ADAPTED1+1+n_shells1; // First star
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
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
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	if (dom == OUTER+3) {
	  // second rect :
	  Array<int> res(2, 1) ;
	  switch (bound) {
		case INNER_BC :
		  res.set(0, 0) = ADAPTED2+1+n_shells2; // Second star
		  res.set(1, 0) = OUTER_BC ;
		  break ;
		case OUTER_BC :
		  res.set(0, 0) = OUTER+5 ; // Compactified domain or first shell
		  res.set(1, 0) = INNER_BC ;
		  break ;
		default :
		  cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;			abort() ;
	}
	return res ;
	}

	if (dom == OUTER+4) {
	  // second chi first :
	  Array<int> res(2, 1) ;
	  switch (bound) {
	    case INNER_BC :
	      res.set(0,0) = ADAPTED2+1+n_shells2; // second star
	      res.set(1,0) = OUTER_BC ;
	      break ;
	    case OUTER_BC :
	      res.set(0,0) = OUTER+5 ; // Compactified domain or first shell
	      res.set(1,0) = INNER_BC ;
	      break ;
	    default :
		cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;
		abort() ;
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
	      cerr << "Bad bound in Space_bin_ns::get_indices_matching_non_std" << endl ;
	      abort() ;
	}
	return res ;
	}

	cerr << "Bad domain in Space_bin_ns::get_indices_matching_non_std" << endl ;
	abort() ;
}

Space_bin_ns::Space_bin_ns (int ttype, double dist, const std::vector<double>& NS1_bounds, const std::vector<double>& NS2_bounds,
                            const std::vector<double>& outer_bounds, int nr, const int n_inner_shellsNS1, const int n_inner_shellsNS2)
                            : n_inner_shells1(n_inner_shellsNS1), n_inner_shells2(n_inner_shellsNS2) {

    ndim = 3 ;

    double rext = outer_bounds[0];
    n_shells1 = int(NS1_bounds.size())-3-n_inner_shells1;
    n_shells2 = int(NS2_bounds.size())-3-n_inner_shells2;

    NS1 = 0;
    NS2 = 3 + n_inner_shells1 + n_shells1;
    ADAPTED1 = NS1 + n_inner_shells1 + 1;
    ADAPTED2 = NS2 + n_inner_shells2 + 1;
    OUTER = 6 + n_shells1 + n_inner_shells1 + n_shells2 + n_inner_shells2;

    n_shells_outer = outer_bounds.size()-1;
    nbr_domains = 12 + n_inner_shells1 + n_shells1 + n_inner_shells2 + n_shells2 + n_shells_outer;
    type_base = ttype ;
    domains = new Domain* [nbr_domains] ;

    Dim_array res (ndim) ;
    res.set(0) = nr ; res.set(1) = nr ; res.set(2) = nr-1 ;
    Dim_array res_bi(ndim) ;
    res_bi.set(0) = nr ; res_bi.set(1) = nr ; res_bi.set(2) = nr ;

    // Bispheric :
    // Computation of aa
    Param par_a ;
    par_a.add_double(NS1_bounds.back(), 0) ;
    par_a.add_double(NS2_bounds.back(), 1) ;
    par_a.add_double(dist,2) ;
    double a_min = 0 ;
    double a_max = dist/2. ;
    double precis = PRECISION ;
    int nitermax = 500 ;
    int niter ;
    double aa = zerosec(func_abns, par_a, a_min, a_max, precis, nitermax, niter) ;
    double eta_plus  =  asinh(aa/NS2_bounds.back()) ;
    double eta_minus = -asinh(aa/NS1_bounds.back()) ;

    double chi_c = 2*atan(aa/rext) ;
    double eta_c = log((1+rext/aa)/(rext/aa-1)) ;
    double eta_lim = eta_c/2. ;
    double chi_lim = chi_lim_eta (eta_lim, rext, aa, chi_c) ;

    auto gen_ns_domains = [&](const int nuc_i, const int adapt_i, auto& bounds, const double& eta, const int& inner_shells, const int& outer_shells) {
        const int RIN  = 0;
        const int RMID = 1 + inner_shells;

        Point center (ndim) ;
        center.set(1) = aa*cosh(eta)/sinh(eta) ;
        domains[nuc_i] = new Domain_nucleus (nuc_i, ttype, bounds[RIN], center, res) ;
        for(int i = RIN; i < inner_shells; ++i) {
            domains[nuc_i+1+i] = new Domain_shell(nuc_i+1+i, ttype, bounds[i], bounds[i+1], center, res);
        }
        domains[adapt_i]   = new Domain_shell_outer_adapted (*this, adapt_i  , ttype, bounds[inner_shells] , bounds[RMID], center, res) ;
        domains[adapt_i+1] = new Domain_shell_inner_adapted (*this, adapt_i+1, ttype, bounds[RMID]   , bounds[RMID+1], center, res) ;

        for(int i = 0; i < outer_shells; ++i)
            domains[adapt_i+2 + i] = new Domain_shell(adapt_i+2+i, ttype, bounds[RMID+1 + i], bounds[RMID+1 + (i+1)], center, res);
    };
    gen_ns_domains(NS1, ADAPTED1, NS1_bounds, eta_minus, n_inner_shells1, n_shells1);
    gen_ns_domains(NS2, ADAPTED2, NS2_bounds, eta_plus, n_inner_shells2, n_shells2);

    // Bispheric part
    domains[OUTER] = new Domain_bispheric_chi_first(OUTER, ttype, aa, eta_minus, rext, chi_lim, res_bi) ;
    domains[OUTER+1] = new Domain_bispheric_rect(OUTER+1, ttype, aa, rext, eta_minus, -eta_lim, chi_lim, res_bi) ;
    domains[OUTER+2] = new Domain_bispheric_eta_first(OUTER+2, ttype, aa, rext, -eta_lim, eta_lim, res_bi) ;
    domains[OUTER+3] = new Domain_bispheric_rect(OUTER+3, ttype, aa, rext, eta_plus, eta_lim, chi_lim, res_bi) ;
    domains[OUTER+4] = new Domain_bispheric_chi_first(OUTER+4, ttype, aa, eta_plus, rext, chi_lim, res_bi) ;

    Point center(3) ;
  	for (int i=0 ; i<n_shells_outer ; i++)
	  	domains[OUTER+5+i] = new Domain_shell (OUTER+5+i, ttype, outer_bounds[i], outer_bounds[i+1], center, res) ;

	  // Compactified
	  domains[OUTER+5+n_shells_outer] = new Domain_compact(OUTER+5+n_shells_outer, ttype, outer_bounds[n_shells_outer], center, res) ;

    const Domain_shell_outer_adapted* pouter_1 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED1]) ;
    pouter_1->vars_to_terms() ;
    pouter_1->update() ;
    const Domain_shell_inner_adapted* pinner_1 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED1+1]) ;
    pinner_1->vars_to_terms() ;
    pinner_1->update() ;
    const Domain_shell_outer_adapted* pouter_2 = dynamic_cast<const Domain_shell_outer_adapted*> (domains[ADAPTED2]) ;
    pouter_2->vars_to_terms() ;
    pouter_2->update() ;
    const Domain_shell_inner_adapted* pinner_2 = dynamic_cast<const Domain_shell_inner_adapted*> (domains[ADAPTED2+1]) ;
    pinner_2->vars_to_terms() ;
    pinner_2->update() ;
}
}
