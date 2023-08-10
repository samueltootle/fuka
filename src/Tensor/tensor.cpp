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

/*
 *  Methods of class Tensor
 *

 */


#include "headcpp.hpp"

#include "scalar.hpp"
#include "tensor.hpp"
#include "vector.hpp"
#include "metric_tensor.hpp"
namespace Kadath {
// Functions standard for storage :
int std_position_array (const Array<int>& idx, int ndim) {
    assert (idx.get_ndim() == 1) ;
    int valence =idx.get_size(0) ;

    for (int i=0 ; i<valence ; i++)
	assert ((idx(i)>=1) && (idx(i)<=ndim)) ;
    int res = 0 ;
    for (int i=0 ; i<valence ; i++)
        res = ndim*res+(idx(i)-1) ;

    return res;
}

int std_position_index (const Index& idx, int ndim) {
   int valence =idx.get_ndim() ;
   for (int i=0 ; i<valence ; i++)
	assert ((idx(i)>=0) && (idx(i)<ndim)) ;
    int res = 0 ;
    for (int i=0 ; i<valence ; i++)
        res = ndim*res+(idx(i)) ;

    return res;
}

Array<int> std_indices (int place, int valence, int ndim) {
    Array<int> res(valence) ;
    for (int i=valence-1 ; i>=0 ; i--) {
		res.set(i) = div(place, ndim).rem ;
		place = int((place-res(i))/ndim) ;
		res.set(i)++ ;
	}
    return res ;
}


			//--------------//
			// Constructors //
			//--------------//

// Standard constructor
// --------------------
Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, const Base_tensor& bb)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
		   valence(val), basis(bb), type_indice(tipe),
		   n_comp(int(pow(espace.get_ndim(), val))), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_size(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;

    // Storage methods :
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}

Tensor::Tensor(const Space& sp, int val, int tipe, const Base_tensor& bb)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
		   valence(val), basis(bb), type_indice(val),
		   n_comp(int(pow(espace.get_ndim(), val))), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    type_indice = tipe ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence);
    name_affected = false ;

    // Storage methods :
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}

Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, const Base_tensor& bb, int dim)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
		   valence(val), basis(bb), type_indice(tipe),
		   n_comp(int(pow(ndim, val))), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_size(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;

    // Storage methods :
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}

Tensor::Tensor(const Space& sp, int val, int tipe, const Base_tensor& bb, int dim)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
		   valence(val), basis(bb), type_indice(val),
		   n_comp(int(pow(ndim, val))), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    type_indice = tipe ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;

    // Storage methods :
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}



// Copy constructor
// ----------------
Tensor::Tensor (const Tensor& source, bool copy) :
    espace(source.espace), ndom(source.ndom), ndim(source.ndim), valence(source.valence), basis(source.basis),
    type_indice(source.type_indice), n_comp (source.n_comp) {

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);
    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(*source.cmp[i], copy) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;
    if ((copy) && (source.name_affected)) {
	name_affected = true ;
	for (int i=0 ; i<valence ; i++)
		name_indice[i] = source.name_indice[i] ;
	}

    // Storage methods :
    give_place_array = source.give_place_array ;
    give_place_index = source.give_place_index ;
    give_indices = source.give_indices ;

    if (source.parameters!=0x0)
      parameters = new Param_tensor(*source.parameters) ;
    else
      parameters = 0x0 ;
}

// THIS IS VERY DANGEROUS
// This assumes the spaces are identical, but different objects
// This also assumes the bases are identical in all domains
Tensor::Tensor (const Space& sp, const Tensor& source) :
    espace(sp), ndom(source.ndom), ndim(source.ndim), valence(source.valence), basis(sp, source.get_basis().get_basis(0)),
    type_indice(source.type_indice), n_comp (source.n_comp) {

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);
    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(*source.cmp[i]) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;
    if ((source.name_affected)) {
        name_affected = true ;
        for (int i=0 ; i<valence ; i++)
            name_indice[i] = source.name_indice[i] ;
	}

    // Storage methods :
    give_place_array = source.give_place_array ;
    give_place_index = source.give_place_index ;
    give_indices = source.give_indices ;

    if (source.parameters!=0x0)
      parameters = new Param_tensor(*source.parameters) ;
    else
      parameters = 0x0 ;
}

//  Constructor for a scalar field: to be used by the derived
//  class {\tt Scalar}
//-----------------------------------------------------------
Tensor::Tensor(const Space& sp) : espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
		valence(0), basis(sp),
		type_indice(0), n_comp(1), parameters(0x0) {

  cmp = MemoryMapper::get_memory<Scalar*>(1);
  cmp[0] = 0x0 ;
  name_indice = 0x0 ;
  name_affected = false ;

   // Storage methods :
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}

Tensor::Tensor(const Space& sp, int val, int tipe, int ncompi, const Base_tensor& bb)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
		   valence(val), basis(bb), type_indice(val),
		   n_comp(ncompi), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    type_indice = tipe ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;
}

Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, int ncompi, const Base_tensor& bb)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
		   valence(val), basis(bb), type_indice(tipe),
		   n_comp(ncompi), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_size(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;
}

Tensor::Tensor(const Space& sp, int val, int tipe, int ncompi, const Base_tensor& bb, int dim)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
		   valence(val), basis(bb), type_indice(val),
		   n_comp(ncompi), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert ((tipe == COV) || (tipe == CON)) ;
    type_indice = tipe ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;
}

Tensor::Tensor(const Space& sp, int val, const Array<int>& tipe, int ncompi, const Base_tensor& bb, int dim)
		: espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
		   valence(val), basis(bb), type_indice(tipe),
		   n_comp(ncompi), parameters(0x0) {

    // Des verifs :
    assert (valence >= 0) ;
    assert (tipe.get_ndim() == 1) ;
    assert (valence == tipe.get_size(0)) ;
    for (int i=0 ; i<valence ; i++)
	assert ((tipe(i) == COV) || (tipe(i) == CON)) ;

    cmp = MemoryMapper::get_memory<Scalar*>(n_comp);

    for (int i=0 ; i<n_comp ; i++)
	cmp[i] = new Scalar(espace) ;

    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
    name_affected = false ;
}

Tensor::Tensor (const Space& sp, FILE* fd) :
	espace(sp), ndom(espace.get_nbr_domains()), ndim(espace.get_ndim()),
	basis(sp, fd), type_indice(fd), parameters(0x0) {

	fread_be (&valence, sizeof(int), 1, fd) ;
	assert (type_indice.get_size(0) == valence) ;
	fread_be (&n_comp, sizeof(int), 1, fd) ;
  cmp = MemoryMapper::get_memory<Scalar*>(n_comp);
	for (int i=0 ; i<n_comp ; i++)
		cmp[i] = new Scalar (espace, fd) ;
	// name of indices are not saved :
    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
        name_affected = false ;

    // Storage methods (overwritten in case of non std things)
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}

Tensor::Tensor (const Space& sp, int dim, FILE* fd) :
	espace(sp), ndom(espace.get_nbr_domains()), ndim(dim),
	basis(sp, fd), type_indice(fd), parameters(0x0) {

	fread_be (&valence, sizeof(int), 1, fd) ;
	assert (type_indice.get_size(0) == valence) ;
	fread_be (&n_comp, sizeof(int), 1, fd) ;
  cmp = MemoryMapper::get_memory<Scalar*>(n_comp);
	for (int i=0 ; i<n_comp ; i++)
		cmp[i] = new Scalar (espace, fd) ;
	// name of indices are not saved :
    name_indice = (valence==0) ? 0x0 : MemoryMapper::get_memory<char>(valence) ;
        name_affected = false ;

    // Storage methods (overwritten in case of non std things)
    give_place_array = std_position_array ;
    give_place_index = std_position_index ;
    give_indices = std_indices ;
}

			//--------------//
			//  Destructor  //
			//--------------//


Tensor::~Tensor () {
    for (int i=0 ; i<n_comp ; i++)
	if (cmp[i]!=0x0)
		delete cmp[i] ;
    MemoryMapper::release_memory<Scalar*>(cmp, n_comp);
    if (name_indice!=0x0)
    MemoryMapper::release_memory<char>(name_indice, valence);
    if (parameters!=0x0)
      delete parameters ;
}

void Tensor::save (FILE* fd) const {
	basis.save(fd) ;
	type_indice.save(fd) ;
	fwrite_be (&valence, sizeof(int), 1, fd) ;
	assert (type_indice.get_size(0) == valence) ;
	fwrite_be (&n_comp, sizeof(int), 1, fd) ;
	for (int i=0 ; i<n_comp ; i++)
		cmp[i]->save(fd) ;
}

void Tensor::annule_hard() {
	for (int i=0 ; i<n_comp ; i++)
		cmp[i]->annule_hard() ;
}

		//-------------
		// Accessors
		//-------------


int Tensor::position (const Array<int>& idx) const {
    return (give_place_array(idx, ndim));
}

int Tensor::position (const Index& idx) const {

    return (give_place_index(idx, ndim)) ;
}

Array<int> Tensor::indices (int place) const {
    return (give_indices(place, valence, ndim)) ;
}

Scalar& Tensor::set() {

    assert (valence == 0) ;
    return *cmp[0] ;
}


// Affectation d'un tenseur d'ordre 1 :
Scalar& Tensor::set(int i) {

    assert (valence == 1) ;

    Array<int> ind (valence) ;
    ind.set(0) = i ;

    int place = position(ind) ;

    return *cmp[place] ;
}

// Affectation d'un tenseur d'ordre 2 :
Scalar& Tensor::set(int ind1, int ind2) {

    assert (valence == 2) ;

    Array<int> ind (valence) ;
    ind.set(0) = ind1 ;
    ind.set(1) = ind2 ;

    int place = position(ind) ;

    return *cmp[place] ;
}

// Affectation d'un tenseur d'ordre 3 :
Scalar& Tensor::set(int ind1, int ind2, int ind3) {

    assert (valence == 3) ;

    Array<int> idx(valence) ;
    idx.set(0) = ind1 ;
    idx.set(1) = ind2 ;
    idx.set(2) = ind3 ;
    int place = position(idx) ;

    return *cmp[place] ;
}


// Affectation d'un tenseur d'ordre 4 :
Scalar& Tensor::set(int ind1, int ind2, int ind3, int ind4) {

    assert (valence == 4) ;

    Array<int> idx(valence) ;
    idx.set(0) = ind1 ;
    idx.set(1) = ind2 ;
    idx.set(2) = ind3 ;
    idx.set(3) = ind4 ;
    int place = position(idx) ;

    return *cmp[place] ;
}


// Affectation cas general
Scalar& Tensor::set(const Array<int>& idx) {

    assert (idx.get_ndim() == 1) ;
    assert (idx.get_size(0) == valence) ;

    int place = position(idx) ;
    return *cmp[place] ;
}

// Affectation cas general from an Index
Scalar& Tensor::set(const Index& idx) {

    Array<int> ind (valence) ;
    for (int i=0 ; i<valence ; i++)
	ind.set(i) = idx(i)+1 ;

    return set(ind) ;
}

const Scalar& Tensor::operator()() const {

    assert(valence == 0) ;

    return *cmp[0] ;

}

const Scalar& Tensor::operator()(int indice) const {

    assert(valence == 1) ;

    Array<int> idx(1) ;
    idx.set(0) = indice ;
    return *cmp[position(idx)] ;

}

const Scalar& Tensor::operator()(int indice1, int indice2) const {

    assert(valence == 2) ;

    Array<int> idx(2) ;
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    return *cmp[position(idx)] ;

}

const Scalar& Tensor::operator()(int indice1, int indice2, int indice3) const {

    assert(valence == 3) ;

    Array<int> idx(3) ;
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    idx.set(2) = indice3 ;
    return *cmp[position(idx)] ;
}


const Scalar& Tensor::operator()(int indice1, int indice2, int indice3,
                                 int indice4) const {

    assert(valence == 4) ;

    Array<int> idx(4) ;
    idx.set(0) = indice1 ;
    idx.set(1) = indice2 ;
    idx.set(2) = indice3 ;
    idx.set(3) = indice4 ;
    return *cmp[position(idx)] ;
}

const Scalar& Tensor::at(int indice1, int indice2) const {
    return operator()(indice1,indice2);

}

const Scalar& Tensor::operator()(const Array<int>& ind) const {

    assert (ind.get_ndim() == 1) ;
    assert (ind.get_size(0) == valence) ;
    return *cmp[position(ind)] ;

}

const Scalar& Tensor::operator()(const Index& idx) const {
    Array<int> ind (valence) ;
    for (int i=0 ; i<valence ; i++)
	ind.set(i) = idx(i)+1 ;

    return operator()(ind) ;

}

void Tensor::set_name_ind (int pos, char name) {
	assert ((pos>=0) && (pos<valence)) ;
	if (!name_affected)
		name_affected = true ;
	name_indice[pos] = name ;
}

// Le cout :
ostream& operator<<(ostream& flux, const Tensor &source ) {

  	flux << '\n' ;
  	flux << "Class : " << typeid(source).name()
        << "           Valence : " << source.valence << '\n' ;

    if (source.valence!=0) {
		flux << source.get_basis() << endl ;
    }

    if (source.valence != 0) {
		flux <<    "Type of the indices : " ;
		for (int i=0 ; i<source.valence ; i++) {
			flux << "index " << i << " : " ;
			if (source.type_indice(i) == CON)
	    		flux << " contravariant." << '\n' ;
			else
	    		flux << " covariant." << '\n' ;
			if ( i < source.valence-1 ) flux << "                      " ;
		}
		flux << '\n' ;
	}

    for (int i=0 ; i<source.n_comp ; i++) {

		if (source.valence == 0) {
			flux <<
			"===================== Scalar field ========================= \n" ;
		}
		else {
      		flux << "================ Component " ;
			Array<int> num_indices (source.indices(i)) ;
			for (int j=0 ; j<source.valence ; j++) {
	  			flux << " " << num_indices(j) ;
      		}
			flux << " ================ \n" ;
		}
		flux << '\n' ;

      	flux << *source.cmp[i] << '\n' ;
    }

	return flux ;
}

// Sets the standard spectal bases of decomposition for each component
void Tensor::std_base() {
	// Loop on the domains
	for (int d=0 ; d<ndom ; d++) {

	switch (valence) {

		case 0 : {
			if (!is_m_quant_affected())
			  cmp[0]->std_base_domain(d) ;
			else
			  cmp[0]->std_base_domain(d, parameters->get_m_quant()) ;
			break ;
		}
		case 1 : {
			bool done = false ;
			if (basis.get_basis(d) ==CARTESIAN_BASIS) {
			  cmp[0]->std_base_x_cart_domain(d) ;
			  cmp[1]->std_base_y_cart_domain(d) ;
			  cmp[2]->std_base_z_cart_domain(d) ;
			  done = true ;
			}
			if (basis.get_basis(d) == SPHERICAL_BASIS) {
			  cmp[0]->std_base_r_spher_domain(d) ;
			  cmp[1]->std_base_t_spher_domain(d) ;
			  cmp[2]->std_base_p_spher_domain(d) ;
			  done = true ;
			  }
			if (basis.get_basis(d) == MTZ_BASIS) {
			  cmp[0]->std_base_r_mtz_domain(d) ;
			  cmp[1]->std_base_t_mtz_domain(d) ;
			  cmp[2]->std_base_p_mtz_domain(d) ;
			  done = true ;
			  }
			if (!done) {
			    cerr << "Tensor::std_base not yet implemented for " << basis << endl ;
			    abort() ;
			}
			break ;
			}
	    default : {
			Vector auxi (espace, CON, basis) ;
			auxi.std_base() ;

			for (int cc=0 ; cc<n_comp ; cc++) {
				  Array<int> ind (indices(cc)) ;
				  Base_spectral result (auxi(ind(0))(d).get_base()) ;
				  for (int i=1 ; i<valence ; i++)
					result =  espace.get_domain(d)->mult(result, auxi(ind(i))(d).get_base()) ;
				  cmp[cc]->set_domain(d).set_base() = result ;

			  }
			}
			break ;
		}
	}
}


bool Tensor::find_indices (const Tensor& tt, Array<int>& perm) const {

	assert (name_affected) ;
	bool res = true ;
	Array<bool> found (valence) ;
	found = false ;
	for (int ncmp=0 ; ncmp<valence ; ncmp++) {
		int pos=0 ;
		bool finloop = false ;
		do {
			if ((!found(pos)) && (tt.name_indice[pos]==name_indice[ncmp])) {
				found.set(pos) = true ;
				perm.set(ncmp) = pos ;
				finloop = true ;
				// Check valence :
				if (type_indice(ncmp) != tt.type_indice(pos))
					res = false ;
			}
			pos ++ ;
			if ((pos == valence) && (!finloop)) {
				finloop = true ;
				res = false ;
			}
		}
		while ((!finloop) && (res));
	}
	return res ;
}

int& Tensor::set_basis (int d) {
    return basis.set_basis(d) ;
}

const Param_tensor* Tensor::get_parameters() const {
   return parameters ;
}

Param_tensor*& Tensor::set_parameters() {
   return parameters ;
}

void Tensor::affect_parameters() {
    if (parameters==0x0)
      parameters = new Param_tensor() ;
}

bool Tensor::is_m_order_affected() const {
      if (parameters==0x0)
	return false ;
      else
	return parameters->m_order_affected ;
}

bool Tensor::is_m_quant_affected() const {
      if (parameters==0x0)
	return false ;
      else
	return parameters->m_quant_affected ;
}

void Tensor::coef() const
{
   Index pos(*this);
   do
   {
      (*this)(pos).coef();
   }while(pos.inc());
}

void Tensor::coef_i() const
{
   Index pos(*this);
   do
   {
      (*this)(pos).coef_i();
   }while(pos.inc());
}
void Tensor::filter_phi(int dom, int ncf)
{
   Index pos(*this);
   do
   {
      set(pos).filter_phi(dom, ncf);
   }while(pos.inc());
}

void Tensor::change_basis_spher_to_cart()
{
   for (int d(0) ; d <= espace.get_nbr_domains() - 1 ; ++d)
   {
      Tensor auxi(espace, get_valence(), get_index_type(), get_basis());
      auxi = espace.get_domain(d)->change_basis_spher_to_cart(d, *this);
      Index pos(auxi);
      do
      {
         set(pos).set_domain(d) = auxi(pos)(d);
      }while(pos.inc());
      set_basis(d) = CARTESIAN_BASIS;
   }
}

void Tensor::change_basis_cart_to_spher()
{
   for (int d(0) ; d <= espace.get_nbr_domains() - 1 ; ++d)
   {
      Tensor auxi(espace, get_valence(), get_index_type(), get_basis());
      auxi = espace.get_domain(d)->change_basis_cart_to_spher(d, *this);
      Index pos(auxi);
      do
      {
         set(pos).set_domain(d) = auxi(pos)(d);
      }while(pos.inc());
      set_basis(d) = SPHERICAL_BASIS;
   }
}


void Tensor::filter (double threshold)  {

   for (int d=0 ; d<espace.get_nbr_domains() ; d++) {
	espace.get_domain(d)->filter(*this, d, threshold) ;
   }

}


}
