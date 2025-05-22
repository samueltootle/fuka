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
 *  Definition of class Vector
 *
 */

#ifndef __VECTOR_HPP_
#define __VECTOR_HPP_

#include "tensor.hpp"

namespace Kadath {

			//-------------------------//
			//       class Vector      //
			//-------------------------//
			

/**
 * A class derived from \c Tensor to deal specificaly with objects of valence 1 (and so also 1-forms).
 * \ingroup fields
 */
class Vector : public Tensor {

    public:

	/** Standard constructor 
	 * @param sp : the \c Space.
	 * @param tipe : the type tensor (COV vs CON).
	 * @param ba : the tensorial basis used.
	 */
	Vector(const Space& sp, int tipe, const Base_tensor& ba) ;

	Vector(const Vector& a) ;       ///< Copy constructor
        Vector (const Space& sp, FILE*) ; ///< Constructor from a file.

	/** Constructor from a \c Tensor .
	 *  The \c Tensor  must be of valence one.
	 */
	Vector(const Tensor& a) ;

	virtual ~Vector() ;			///< Destructor

	/**
	 * @brief Use with Caution! Construct a new Vector object by copying
	 * the coefficients/values/etc, but to a new Tensor with a different Space
	 * pointer.  THIS ASSUMES BOTH SPACES ARE THE SAME, but with different
	 * memory addresses
	 * 
	 * @param sp New space
	 * @param a Tensor to copy
	 */
	Vector(const Space& sp, const Tensor& a) ;
	Vector(const Space& sp, const Vector& a) ;       ///< Copy constructor

    // Mutators / assignment
    // ---------------------
    public:

	void operator=(const Vector&) ; ///< Assignment to another \c Vector.	
	virtual void operator=(const Tensor&) ;	
	virtual void operator=(double) ;
	virtual void annule_hard() ;

    // Accessors
    // ---------
    public:
        using Tensor::set;
	Scalar& set(int ) ; ///< Read/write access to a component
	const Scalar& operator()(int ) const; ///<Readonly access to a component
	const Scalar& at(int) const; ///<Readonly access to a component

	
	virtual int position(const Array<int>& idx) const {
	  assert (idx.get_ndim() == 1) ;
	  assert (idx.get_size(0) == 1) ;
	  assert ((idx(0) >= 1) && (idx(0) <= espace.get_ndim())) ;

	  return (idx(0) - 1) ;
	} ;

	virtual int position(const Index& idx) const {
		assert ((idx(0)>=0) && (idx(0)<ndim)) ;
		return (idx(0)) ;
	}

	/**
	 * Returns the type of the objects (CON or COV)
	 */
	int get_index_type() const {return type_indice(0) ;};

	virtual Array<int> indices(int place) const {
	  assert((place>=0) && (place<espace.get_ndim())) ;
	  
	  Array<int> res(1) ;
	  res = place + 1;
	  return res ;
	};
};
}
#endif
