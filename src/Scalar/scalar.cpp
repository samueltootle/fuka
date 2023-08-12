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
#include "scalar.hpp"

namespace Kadath {
Scalar::Scalar (const Space& sp) : Tensor(sp) {
	val_zones = MemoryMapper::get_memory<Val_domain*>(ndom);

	for (int l=0 ; l<ndom ; l++)
	    val_zones[l] = new Val_domain(sp.get_domain(l)) ;
	cmp[0] = this ;
}

Scalar::Scalar (const Scalar& so, bool copie) : Tensor(so.espace) {

	val_zones = MemoryMapper::get_memory<Val_domain*>(ndom);
	for (int l=0 ; l<ndom ; l++)
	    val_zones[l] = new Val_domain(*so.val_zones[l], copie) ;
	cmp[0] = this ;
}

Scalar::Scalar (const Space& sp, const Scalar& so) : Tensor(sp) {

	val_zones = MemoryMapper::get_memory<Val_domain*>(ndom);
	for (int l=0 ; l<ndom ; l++)
	    val_zones[l] = new Val_domain(sp.get_domain(l), *so.val_zones[l]) ;
	cmp[0] = this ;
}

Scalar::Scalar (const Tensor& so, bool copie) : Tensor(so.espace) {

	assert (so.valence==0) ;

	val_zones = MemoryMapper::get_memory<Val_domain*>(ndom);
	for (int l=0 ; l<ndom ; l++)
	    val_zones[l] = new Val_domain(*so.cmp[0]->val_zones[l], copie) ;
	cmp[0] = this ;
}

Scalar::Scalar (const Space& sp, FILE* fd) : Tensor(sp) {
	val_zones = MemoryMapper::get_memory<Val_domain*>(ndom);
	for (int l=0 ; l<ndom ; l++)
	    val_zones[l] = new Val_domain(sp.get_domain(l), fd) ;
	cmp[0] = this ;
}

Scalar::~Scalar () {
	for (int i=0 ; i<ndom ; i++)
		delete val_zones[i] ;
  MemoryMapper::release_memory<Val_domain*>(val_zones,ndom);
	cmp[0] = 0x0 ;
}

void Scalar::save (FILE* fd) const {
	for (int i=0 ; i<get_nbr_domains() ; i++)
		val_zones[i]->save(fd) ;
}

void Scalar::operator= (const Scalar& so)  {
	assert (&espace==&so.espace) ;
	for (int i=0 ; i<ndom ; i++) {
		if (val_zones[i] != 0x0)
			delete val_zones[i] ;
		val_zones[i] = new Val_domain(*so.val_zones[i]) ;
	}
}

void Scalar::operator= (const Tensor& so)  {
	assert (&espace==&so.espace) ;
	assert (so.valence==0) ;
	for (int i=0 ; i<ndom ; i++) {
		if (val_zones[i] != 0x0)
			delete val_zones[i] ;
		val_zones[i] = new Val_domain(*so.cmp[0]->val_zones[i]) ;
	}
}

void Scalar::operator= (double xx) {
	for (int i=0 ; i<ndom ; i++)
		set_domain(i) = xx ;
}

void Scalar::annule_hard() {
	for (int i=0 ; i<ndom ; i++)
		set_domain(i).annule_hard() ;
}

void Scalar::annule_hard_coef() {
	for (int i=0 ; i<ndom ; i++)
		set_domain(i).annule_hard_coef() ;
}

void Scalar::set_in_conf() {
	for (int l=0 ; l<ndom ; l++)
	    val_zones[l]->set_in_conf() ;
}

void Scalar::set_in_coef() {
	for (int l=0 ; l<ndom ; l++)
	    val_zones[l]->set_in_coef() ;
}

void Scalar::std_base() {
    if (!is_m_quant_affected()) {
      for (int l=0 ; l<ndom ; l++)
	    val_zones[l]->std_base() ;
    }
    else {
      for (int l=0 ; l<ndom ; l++)
	    val_zones[l]->std_base(parameters->get_m_quant()) ;
    }
}

void Scalar::std_anti_base() {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_anti_base() ;
}

void Scalar::std_base(int m) {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_base(m) ;
}

void Scalar::std_anti_base(int m) {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_anti_base(m) ;
}

void Scalar::std_base_r_spher() {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_base_r_spher() ;
}

void Scalar::std_base_t_spher() {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_base_t_spher() ;
}

void Scalar::std_base_p_spher() {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_base_p_spher() ;
}

void Scalar::std_base_domain(int d) {
    val_zones[d]->std_base() ;
}

void Scalar::std_anti_base_domain(int d) {
    val_zones[d]->std_anti_base() ;
}

void Scalar::std_base_domain(int d, int m) {
    val_zones[d]->std_base(m) ;
}

void Scalar::std_base_x_cart_domain(int d) {
    val_zones[d]->std_base_x_cart() ;
}

void Scalar::std_base_y_cart_domain(int d) {
    val_zones[d]->std_base_y_cart() ;
}

void Scalar::std_base_z_cart_domain(int d) {
    val_zones[d]->std_base_z_cart() ;
}

void Scalar::std_base_r_spher_domain(int d) {
    val_zones[d]->std_base_r_spher() ;
}

void Scalar::std_base_t_spher_domain(int d) {
    val_zones[d]->std_base_t_spher() ;
}

void Scalar::std_base_p_spher_domain(int d) {
    val_zones[d]->std_base_p_spher() ;
}

void Scalar::std_base_xy_cart_domain(int d) {
    val_zones[d]->std_base_xy_cart() ;
}

void Scalar::std_base_xz_cart_domain(int d) {
    val_zones[d]->std_base_xz_cart() ;
}

void Scalar::std_base_yz_cart_domain(int d) {
    val_zones[d]->std_base_yz_cart() ;
}

void Scalar::std_base_rt_spher_domain(int d) {
    val_zones[d]->std_base_rt_spher() ;
}

void Scalar::std_base_rp_spher_domain(int d) {
    val_zones[d]->std_base_rp_spher() ;
}

void Scalar::std_base_tp_spher_domain(int d) {
    val_zones[d]->std_base_tp_spher() ;
}

void Scalar::std_xodd_base() {
    val_zones[0]->std_xodd_base() ;
    for (int l=1 ; l<ndom ; l++)
	   val_zones[l]->std_base() ;
}

void Scalar::std_todd_base() {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_todd_base() ;
}

void Scalar::std_xodd_todd_base() {
    val_zones[0]->std_xodd_todd_base() ;
    for (int l=1 ; l<ndom ; l++)
	   val_zones[l]->std_todd_base() ;
}

void Scalar::std_base_odd () {
    for (int l=0 ; l<ndom ; l++)
	   val_zones[l]->std_base_odd() ;
}

void Scalar::std_base_r_mtz_domain(int d) {
    val_zones[d]->std_base_r_mtz() ;
}

void Scalar::std_base_t_mtz_domain(int d) {
    val_zones[d]->std_base_t_mtz() ;
}

void Scalar::std_base_p_mtz_domain(int d) {
    val_zones[d]->std_base_p_mtz() ;
}


void Scalar::allocate_conf() {
 	for (int l=0 ; l<ndom ; l++)
	    val_zones[l]->allocate_conf() ;
}

void Scalar::allocate_coef() {
 	for (int l=0 ; l<ndom ; l++)
	    val_zones[l]->allocate_coef() ;
}


const Val_domain& Scalar::operator() (int i) const {
	assert ((i>=0) && (i<ndom)) ;
	return (*val_zones[i]) ;
}
const Val_domain& Scalar::at(int i) const {
	return operator()(i);
}


Val_domain& Scalar::set_domain (int l) {
	assert ((l>=0) && (l<ndom)) ;
	return *val_zones[l] ;
}

void Scalar::set_val_inf(double x) {
      val_zones[ndom-1]->get_domain()->set_val_inf(*val_zones[ndom-1], x) ;
}

void Scalar::set_val_inf(double x, int l) {
      val_zones[l]->get_domain()->set_val_inf(*val_zones[l], x) ;
}

double Scalar::val_point(const Point& xx, int sens) const {

	assert ((sens==+1) || (sens==-1)) ;

	// bool* inside = MemoryMapper::get_memory<bool>(ndom);

	// for (int l=ndom-1 ; l>=0 ; l--)
	    // inside[l] = get_domain(l)->is_in(xx) ;
	// First domain in which the point is :
	int ld = -1 ;
	// if (sens == -1) {
	//   for (int l=ndom-1 ; l>=0 ; l--)
	    //   if ((ld==-1) && (inside[l]))
		//   ld = l ;
	// }
	// else {
	//    for (int l=0 ; l<ndom ; l++)
	//       if ((ld==-1) && (inside[l]))
	// 	  ld = l ;
	// }
	for (int l=0 ; l<ndom ; l++) {
		bool isin = get_domain(l)->is_in(xx);
		if (isin) {
			ld = l ;
			break;
		}
	}
		  
	if (ld==-1) {
	     cout << "Point " << xx << "not found in the computational space..." << endl ;
	     abort() ;
	}
	else {
	if (val_zones[ld]->check_if_zero())
		return 0. ;
	else {
	Point num(get_domain(ld)->absol_to_num(xx)) ;

	coef() ;
//   MemoryMapper::release_memory<bool>(inside, ndom);
	return val_zones[ld]->base.summation(num, *val_zones[ld]->cf) ;
	}
	}
}

double Scalar::val_point_zeronotdef(const Point& xx, int sens) const {

	assert ((sens==+1) || (sens==-1)) ;

	bool* inside = MemoryMapper::get_memory<bool>(ndom);
	for (int l=ndom-1 ; l>=0 ; l--)
	    inside[l] = get_domain(l)->is_in(xx) ;
	// First domain in which the point is :
	int ld = -1 ;
	if (sens == -1) {
	  for (int l=ndom-1 ; l>=0 ; l--)
	      if ((ld==-1) && (inside[l]))
		  ld = l ;
	}
	else {
	   for (int l=0 ; l<ndom ; l++)
	      if ((ld==-1) && (inside[l]))
		  ld = l ;
	}

	if (ld==-1) {
	    return 0 ;
	}
	else {
	if (val_zones[ld]->check_if_zero())
		return 0. ;
	else {
	Point num(get_domain(ld)->absol_to_num(xx)) ;

	coef() ;
  MemoryMapper::release_memory<bool>(inside, ndom);
	return val_zones[ld]->base.summation(num, *val_zones[ld]->cf) ;
	}
	}
}
void Scalar::coef() const {
    // Boucles sur les domaines
    for (int l=0; l<ndom; l++)
	val_zones[l]->coef() ;
}

void Scalar::coef_i() const {
    // Boucles sur les domaines
    for (int l=0; l<ndom; l++)
	val_zones[l]->coef_i() ;
}


void Scalar::filter_phi (int dom, int ncf) {
	coef() ;
	Index pcf (get_space().get_domain(dom)->get_nbr_coefs()) ;
	int np = get_space().get_domain(dom)->get_nbr_coefs()(2) ;
	do {
		if (pcf(2)+ncf > np-1)
			set_domain(dom).set_coef(pcf) = 0 ;
	}
	while (pcf.inc()) ;
}


Scalar Scalar::der_var(int var) const {
	Scalar res (*this, false) ;
	for (int dom=0 ; dom<ndom ; dom++)
		res.set_domain(dom) = val_zones[dom]->der_var(var) ;
	return res ;
}

Scalar Scalar::der_abs(int var) const {
	Scalar res (*this, false) ;
	for (int dom=0 ; dom<ndom ; dom++)
		res.set_domain(dom) = val_zones[dom]->der_abs(var) ;

	return res ;
}

Scalar Scalar::der_spher(int var) const
{
   Scalar res(*this, false) ;
   for (int dom(0) ; dom < ndom ; ++dom)
      res.set_domain(dom) = val_zones[dom]->der_spher(var) ;
   return res ;
}

Scalar Scalar::der_r() const {
	Scalar res (*this, false) ;
	for (int dom=0 ; dom<ndom ; dom++)
		res.set_domain(dom) = val_zones[dom]->der_r() ;

	return res ;
}

double Scalar::integrale() const {
  double res = 0 ;
  for (int dom=0 ; dom<ndom ; dom++)
    res += val_zones[dom]->integrale() ;

	return res ;
}

unique_ptr<Scalar> Scalar::clone() const
{
   unique_ptr<Scalar> res(new Scalar(*this));   // USE make_unique when g++ 4.9 and std=c++14 available
   return res;
}

ostream& operator<< (ostream& o, const Scalar& so) {

    o << "Scalar" << endl ;
    for (int l=0 ; l<so.ndom ; l++) {
    	o << "Domain : " << l << endl ;
	o << *so.val_zones[l] << endl ;
    }
    return o ;
}

double Scalar::integ_volume() const {
  double res = 0 ;
  for (int dom=0 ; dom<ndom ; dom++)
    res += val_zones[dom]->integ_volume() ;

	return res ;
}

Scalar Scalar::zero(Space const& espace)
{
   Scalar res(espace);
   res = 0.0;
   return res;
}
double Scalar::val_point_bound(const Point& xx, int bound) const {

	bool* inside = MemoryMapper::get_memory<bool>(ndom);

	for (int l=ndom-1 ; l>=0 ; l--)
	    inside[l] = get_domain(l)->is_in(xx) ;
	// First domain in which the point is :
	int ld = -1 ;
  for (int l=0 ; l<ndom ; l++)
    if ((ld==-1) && (inside[l]))
      ld = l ;
	
  if (ld==-1) {
	     cout << "Point " << xx << "not found in the computational space..." << endl ;
	     abort() ;
	}
	else {
    if (val_zones[ld]->check_if_zero())
      return 0. ;
    else {
      Point num(get_domain(ld)->absol_to_num_bound(xx, bound)) ;

      coef() ;
      MemoryMapper::release_memory<bool>(inside, ndom);
      return val_zones[ld]->base.summation(num, *val_zones[ld]->cf) ;
	  }
	}
}


}
