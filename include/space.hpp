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

#ifndef __SPACE_HPP_
#define __SPACE_HPP_

#include "headcpp.hpp"
#include "base_spectral.hpp"
#include "array.hpp"
#include "param.hpp"
#include "utilities.hpp"
#include "memory.hpp"

#define CHEB_TYPE 1
#define LEG_TYPE 2



#define OUTER_BC 1
#define INNER_BC 2
#define CHI_ONE_BC 3
#define ETA_PLUS_BC 4
#define ETA_MINUS_BC 5
#define TIME_INIT 6

#include "point.hpp"

namespace Kadath  {
class Space ;
class Scalar ;
class System_of_eqs ;
class Val_domain ;
class Term_eq ;
class Base_tensor ;
class Metric ;
class Tensor ;

/**
* Abstract class that implements the fonctionnalities common to all the type of domains.
*
* In particular, it relates the numerical coordinates (for which the spectral expansion is done), to the
* physical ones.
* \ingroup domain
*/

class Domain : public MemoryMappable {

 protected:
  int num_dom ; ///< Number of the current domain (used by the \c Space)
  int ndim ; ///< Number of dimensions.
  Dim_array nbr_points ; ///< Number of colocation points.
  mutable Dim_array nbr_coefs ; ///< Number of coefficients.

  /**
  * Type of colocation point :
  * \li \c CHEB_TYPE : Gauss-Lobato of Chebyshev polynomials.
  * \li \c LEG_TYPE : Gauss-Lobato of Legendre polynomials.
  */
  int type_base ;

  /**
   * Colocation points in each dimension (stored in \c ndim 1d- arrays)
   */
  Array<double>** coloc ;
  mutable Val_domain** absol ; ///< Asbolute coordinates (if defined ; usually Cartesian-like)
  mutable Val_domain** cart ; ///< Cartesian coordinates
  mutable Val_domain* radius ; ///< The generalized radius.
  mutable Val_domain** cart_surr ; ///< Cartesian coordinates divided by the radius

  explicit Domain (int num, int ttype, const Dim_array& res) ; ///< Constructor from a number of points and a type of base
  explicit Domain (int, FILE*) ; ///< Constructor from a file
  Domain (const Domain& so) ; ///< Copy constructor.
  Domain (const Domain& so, bool import) ; ///< Copy constructor.

 public:
  virtual ~Domain() ; ///< Destructor.
  virtual void save (FILE*) const ; ///< Saving function
   int get_num() const  /// Returns the index of the curent domain.
	{return num_dom ;} ;
   Dim_array get_nbr_points() const /// Returns the number of points.
	{return nbr_points ;} ;
   Dim_array get_nbr_coefs() const  /// Returns the number of coefficients.
	{return nbr_coefs ;} ;
   int get_ndim() const /// Returns the number of dimensions.
	{return ndim; };
   int get_type_base() const /// Returns the type of the basis.
	{return type_base ;};
   Array<double> get_coloc (int) const ; ///< Returns the colocation points for a given variable.
   virtual Point get_center() const ; ///< Returns the center.
   virtual Val_domain get_chi() const ; ///< Returns the variable \f$ \chi \f$.
   virtual Val_domain get_eta() const ; ///< Returns the variable \f$ \eta \f$.
   virtual Val_domain get_X() const ; ///< Returns the variable \f$ X \f$.
   virtual Val_domain get_T() const ; ///< Returns the variable \f$ T \f$.

 public:
    Val_domain get_absol(int i) const ; ///<Returns the absolute coordinates
    Val_domain get_radius() const ; ///< Returns the generalized radius.
    /**
    * Returns a Cartesian coordinates.
    * @param i [input] : composant (\f$ 0 \leq i < \f$ \c ndim).
    */
    Val_domain get_cart(int i) const ;
/**
    * Returns a Cartesian coordinates divided by the radius
    * @param i [input] : composant (\f$ 0 \leq i < \f$ \c ndim).
    */
    Val_domain get_cart_surr(int i) const ;

 protected:
    /**
    * Destroys the derivated members (like \c coloc, \c cart and \c radius),
    * when changing the type of colocation points.
    */
    virtual void del_deriv() const ;
    virtual void do_radius () const ; ///< Computes the generalized radius.
    virtual void do_cart () const ; ///< Computes the Cartesian coordinates.
    virtual void do_cart_surr () const ; ///< Computes the Cartesian coordinates over the radius
    virtual void do_absol() const ; ///< Computes the absolute coordinates

 public:
     void operator= (const Domain&) ; ///< Assignement operator.

  private:
 /**
    * Gives the standard base for Chebyshev polynomials.
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base(Base_spectral& so) const ;
    /**
    * Gives the standard base for Legendre polynomials.
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base(Base_spectral& so) const ;
    /**
    * Gives the base using Chebyshev polynomials, for functions antisymetric with respect to \f$z\f$
    * @param so [intput] : the returned base.
    */
     virtual void set_anti_cheb_base(Base_spectral& so) const ;
   /**
    * Gives the base using Legendre polynomials, for functions antisymetric with respect to \f$z\f$
    * @param so [input] : the returned base.
    */
     virtual void set_anti_legendre_base(Base_spectral& so) const ;
    /**
    * Gives the standard base using Chebyshev polynomials. Used for a spherical harmonic \f$ m\f$.
    * @param so [input] : the returned base.
    * @param m [input] : index of the spherical harmonic.
    */
     virtual void set_cheb_base_with_m(Base_spectral& so, int m) const ;
    /**
    * Gives the stnadard base using Legendre polynomials. Used for a spherical harmonic \f$ m\f$.
    * @param so [input] : the returned base.
     * @param m [input] : index of the spherical harmonic.
    */
     virtual void set_legendre_base_with_m(Base_spectral& so, int m) const ;
    /**
    * Gives the base using Chebyshev polynomials, for functions antisymetric with respect to \f$z\f$ . Used for a spherical harmonic \f$ m\f$.
    * @param so [input] : the returned base.
    * @param m [input] : index of the spherical harmonic.
    */
     virtual void set_anti_cheb_base_with_m(Base_spectral& so, int m) const ;
   /**
    * Gives the base using Legendre polynomials, for functions antisymetric with respect to \f$z\f$. Used for a spherical harmonic \f$ m\f$.
    * @param so [input] : the returned base.
    * @param m [input] : index of the spherical harmonic.
    */
     virtual void set_anti_legendre_base_with_m(Base_spectral& so, int m) const ;


 /**
    * Gives the base using Chebyshev polynomials, for the radial component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_r_spher(Base_spectral& so) const ;

/**
    * Gives the base using Chebyshev polynomials, for the \f$ \theta\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_t_spher(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ \varphi\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_p_spher(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the radial component of a vector in the MTZ setting
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_r_mtz(Base_spectral& so) const ;

/**
    * Gives the base using Chebyshev polynomials, for the \f$ \theta\f$ component of a vector in the MTZ setting
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_t_mtz(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ \varphi\f$ component of a vector in the MTZ setting
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_p_mtz(Base_spectral& so) const ;

/**
    * Gives the base using Chebyshev polynomials, for the \f$(r, \theta)\f$ component of a 2-tensor
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_rt_spher(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$(r, \varphi)\f$ component of a 2-tensor
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_rp_spher(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$(\theta, \varphi)\f$ component of a 2-tensor
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_tp_spher(Base_spectral& so) const ;

	 /**
    * Gives the base using Legendre polynomials, for the radial component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_r_spher(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ \theta\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_t_spher(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ \varphi\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_p_spher(Base_spectral& so) const ;
  	 /**
    * Gives the base using Legendre polynomials, for the radial component of a vector in the MTZ context
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_r_mtz(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ \theta\f$ component of a vector in the MTZ context
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_t_mtz(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ \varphi\f$ component of a vector in the MTZ context
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_p_mtz(Base_spectral& so) const ;

	/**
    * Gives the base using Chebyshev polynomials, for the \f$ x\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_x_cart(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ y\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_y_cart(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ z\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_z_cart(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ (x,y) \f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_xy_cart(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ (x,z) \f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_xz_cart(Base_spectral& so) const ;
/**
    * Gives the base using Chebyshev polynomials, for the \f$ (y, z) \f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_base_yz_cart(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ x\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_x_cart(Base_spectral& so ) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ y\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_y_cart(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for the \f$ z\f$ component of a vector
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_base_z_cart(Base_spectral& so) const ;

	/**
    * Gives the base using Chebyshev polynomials, for odd functions in\f$ X\f$ (critic space case)
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_xodd_base(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for odd functions in\f$ X\f$ (critic space case)
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_xodd_base(Base_spectral&) const ;
    /**
    * Gives the base using Chebyshev polynomials, for odd functions in\f$ T\f$ (critic space case)
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_todd_base(Base_spectral& so) const ;
/**
    * Gives the base using Legendre polynomials, for odd functions in\f$ T\f$ (critic space case)
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_todd_base(Base_spectral&) const ;
 /**
    * Gives the base using Chebyshev polynomials, for odd functions in\f$ X\f$ and in \f$ T\f$ (critic space case)
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_xodd_todd_base(Base_spectral& so) const ;
 /**
    * Gives the base using Chebyshev polynomials, for odd functions in\f$ X\f$ and in \f$ T\f$(critic space case)
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_xodd_todd_base(Base_spectral& so) const ;

   /**
    * Gives the base using odd Chebyshev polynomials$
    * @param so [input] : the returned base.
    */
    virtual void set_cheb_base_odd(Base_spectral& so) const ;
 /**
    * Gives the base using odd Legendre polynomials$
    * @param so [input] : the returned base.
    */
    virtual void set_legendre_base_odd(Base_spectral&) const ;

/**
    * Gives the base using odd Chebyshev polynomials$ for the radius
    * @param so [input] : the returned base.
    */
     virtual void set_cheb_r_base(Base_spectral& so) const ;
/**
    * Gives the base using odd Legendre polynomials$ for the radius
    * @param so [input] : the returned base.
    */
     virtual void set_legendre_r_base(Base_spectral& so) const ;

     virtual void do_coloc () ; ///< Computes the colocation points.

   public:
   /**
     * Check whether a point lies inside \c Domain.
     * @param xx [input] : the point.
     * @param prec [input] : precision of the computation (used when comparing doubles).
     * @returns a \c true if the point is in the domain and \c false otherwise.
     */
     virtual bool is_in(const Point& xx, double prec=1e-12) const ;
     /**
     * Computes the numerical coordinates from the physical ones.
     * @param xxx [input] : the absolute Cartesian coordinates of the point.
     * @returns the numerical coordinates.
     */
     virtual const Point absol_to_num(const Point& xxx) const ;
 /**
     * Computes the numerical coordinates from the physical ones for a point lying on a boundary
     * @param xxx [input] : the absolute Cartesian coordinates of the point.
     * @param bound [input] : the boundary.
     * @returns the numerical coordinates.
     */
     virtual const Point absol_to_num_bound(const Point& xxx, int bound) const ;
     /**
     * Computes the derivative with respect to the absolute Cartesian coordinates from the
     * derivative with respect to the numerical coordinates.
     * @param der_var [input] : the \c ndim derivatives with respect to the numerical coordinates.
     * @param der_abs [output] : the \c ndim derivatives with respect to the absolute Cartesian coordinates.
     */
     virtual void do_der_abs_from_der_var(Val_domain** der_var, Val_domain** der_abs) const ;
     /**
     * Method for the multiplication of two \c Base_spectral.
     * @returns the output base is undefined if the result is not implemented (i.e. if one tries to multiply cosines with Chebyshev polynomials for instance).
     */
     virtual Base_spectral mult (const Base_spectral&, const Base_spectral&) const ;

  public:

    virtual double get_rmin() const ; ///< Returns the minimum radius.
    virtual double get_rmax() const ; ///< Returns the maximum radius.

    /**
     * Multiplication by \f$ \cos \varphi\f$.
     */
     virtual Val_domain mult_cos_phi (const Val_domain&) const ;
    /**
     * Multiplication by \f$ \sin \varphi\f$.
     */
     virtual Val_domain mult_sin_phi (const Val_domain&) const ;
    /**
     * Multiplication by \f$ \cos \theta\f$.
     */
     virtual Val_domain mult_cos_theta (const Val_domain&) const ;
    /**
     * Multiplication by \f$ \sin \theta\f$.
     */
     virtual Val_domain mult_sin_theta (const Val_domain&) const ;
    /**
     * Division by \f$ \sin \theta\f$.
     */
     virtual Val_domain div_sin_theta (const Val_domain&) const ;
     /**
     * Division by \f$ \cos \theta\f$.
     */
     virtual Val_domain div_cos_theta (const Val_domain&) const ;
    /**
     * Division by \f$ x\f$.
     */
     virtual Val_domain div_x (const Val_domain&) const ;
    /**
     * Division by \f$ \chi\f$.
     */
     virtual Val_domain div_chi (const Val_domain&) const ;
    /**
     * Division by \f$ (x-1)\f$.
     */
     virtual Val_domain div_xm1 (const Val_domain&) const ;
 /**
     * Division by \f$ (1-x^2)\f$.
     */
     virtual Val_domain div_1mx2 (const Val_domain&) const ;
      /**
     * Division by \f$ (x+1)\f$.
     */
     virtual Val_domain div_xp1 (const Val_domain&) const ;
    /**
     * Multiplication by \f$ (x-1)\f$.
     */
     virtual Val_domain mult_xm1 (const Val_domain&) const ;
    /**
     * Division by \f$ \sin \chi\f$ .
     */
     virtual Val_domain div_sin_chi (const Val_domain&) const ;
      /**
     * Multiplication by \f$ \cos \omega t\f$.
     */
     virtual Val_domain mult_cos_time (const Val_domain&) const ;
    /**
     * Multiplication by \f$ \sin \omega t\f$.
     */
     virtual Val_domain mult_sin_time (const Val_domain&) const ;

	/**
	* Changes the tensorial basis from Cartsian to spherical in a given domain.
	* @param dd [input] : the domain. Should be consistent with *this.
	* @param so [input] : the input tensor.
	* @returns the tensor in spherical tensorial basis in the current domain.
	*/
     virtual Tensor change_basis_cart_to_spher (int dd, const Tensor& so) const ;

/**
	* Changes the tensorial basis from spherical to Cartesian in a given domain.
	* @param dd [input] : the domain. Should be consistent with *this.
	* @param so [input] : the input tensor.
	* @returns the tensor in Cartesian tensorial basis in the current domain.
	*/
     virtual Tensor change_basis_spher_to_cart (int dd, const Tensor&) const ;

	/**
	* Computes the ordinary flat Laplacian for a scalar field with an harmonic index \c m.
	* @param so [input] : the input scalar field.
	* @param m [input] : harmonic index.
	* @returns the Laplacian.
	*/
     virtual Val_domain laplacian (const Val_domain& so, int m) const ;
/**
	* Computes the ordinary flat 2dè- Laplacian for a scalar field with an harmonic index \c m.
	* @param so [input] : the input scalar field.
	* @param m [input] : harmonic index.
	* @returns the 2d-Laplacian.
	*/
     virtual Val_domain laplacian2 (const Val_domain& so, int m) const ;

	/**
	* Compute the radial derivative of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the derivative.
	*/
     virtual Val_domain der_r (const Val_domain&) const ;

/**
	* Compute the derivative with respect to \f$ \theta\f$ of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the derivative.
	*/
     virtual Val_domain der_t (const Val_domain&) const ;

/**
	* Compute the derivative with respect to \f$ \varphi\f$ of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the derivative.
	*/
     virtual Val_domain der_p (const Val_domain&) const ;

/**
	* Compute the radial derivative multiplied by \f$ r^2\f$ of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the result.
	*/
     virtual Val_domain der_r_rtwo (const Val_domain& so) const ;
/**
	* Compute the \f$ f' / r\f$ of a scalar field \f$ f\f$.
	* @param so [input] : the input scalar field.
	* @returns the result.
	*/
     virtual Val_domain srdr (const Val_domain& so) const ;
/**
	* Compute the second radial derivative w of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the result.
	*/
     virtual Val_domain ddr (const Val_domain&) const ;
/**
	* Compute the second derivative with respect to \f$\varphi\f$ of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the result.
	*/
     virtual Val_domain ddp (const Val_domain&) const ;
/**
	* Compute the second derivative with respect to \f$\theta\f$ of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the result.
	*/
     virtual Val_domain ddt (const Val_domain&) const ;
/**
	* Compute the derivative with respect to \f$\theta\f$ of a scalar field.
	* @param so [input] : the input scalar field.
	* @returns the result.
	*/
     virtual Val_domain dt (const Val_domain&) const ;


	/**
	* Computes the time derivative of a field.
	* @param so [input] : the input field.
	* @returns the result.
	*/
     virtual Val_domain dtime (const Val_domain&) const ;

	/**
	* Computes the second time derivative of a field
	* @param so [input] : the input field.
	* @returns the result
	*/
	virtual Val_domain ddtime (const Val_domain&) const ;


	/**
	* Returns the vector normal to a surface. Must be a \c Term_eq because the dmain can be a variable one.
	* @param bound [input] : the boundary where the normal is computed.
	* @param tipe [input] : tensorial coordinates used. (Cartesian or spherical basis)
	* @returns  the normal derivative.
	*/
     virtual const Term_eq* give_normal (int bound, int tipe) const ;

     // Multipoles extraction
	/**
	* Extraction of a given multipole, at some boundary, for a  symmetric scalar function.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param bound [input] : the boundary at which the computation is done.
	* @param so [input] : input scalar field.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the multipolar coefficient.
	*/
     virtual double multipoles_sym (int k, int j, int, const Val_domain& so, const Array<double>& passage) const ;

	/**
	* Extraction of a given multipole, at some boundary, for a anti-symmetric scalar function.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param bound [input] : the boundary at which the computation is done.
	* @param so [input] : input scalar field.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the multipolar coefficient.
	*/
     virtual double multipoles_asym (int, int, int, const Val_domain&, const Array<double>&) const ;

  	/**
	* Extraction of a given multipole, at some boundary, for a symmetric scalar function.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param bound [input] : the boundary at which the computation is done.
	* @param so [input] : input scalar field.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the multipolar coefficient.
	*/
     virtual Term_eq multipoles_sym (int k, int j, int bound, const Term_eq& so, const Array<double>& passage) const ;
   /**
	* Extraction of a given multipole, at some boundary, for an anti-symmetric scalar function.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param bound [input] : the boundary at which the computation is done.
	* @param so [input] : input scalar field.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the multipolar coefficient.
	*/
     virtual Term_eq multipoles_asym (int k, int j, int bound, const Term_eq& so, const Array<double>& passage) const ;

	 /**
	* Gives some radial fit for a given multipole, intended for symmetric scalar function.
	* @param space [input] : the concerned \c Space.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param omega [input] : angular velocity.
	* @param f [input] : pointer on the radial function describing the fit.
	* @param param [input] : parameters of the radial fit.
	* @returns the radial fit.
	*/
     virtual Term_eq radial_part_sym (const Space& space, int k, int j, const Term_eq& omega, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param& param) const ;
 	/**
	* Gives some radial fit for a given multipole, intended for anti-symmetric scalar function.
	* @param space [input] : the concerned \c Space.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param omega [input] : angular velocity.
	* @param f [input] : pointer on the radial function describing the fit.
	* @param param [input] : parameters of the radial fit.
	* @returns the radial fit.
	*/
     virtual Term_eq radial_part_asym (const Space& space, int k, int j, const Term_eq& omega, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param& param) const ;

	/**
	* Fit, spherical harmonic by spherical harmonic, for a symmetric function.
	* @param so [input] : input scalar field.
	* @param omega [input] : angular velocity.
	* @param bound [input] : the boundary at which the computation is done.
	* @param f [input] : pointer on the radial function describing the radial fit.
	* @param param [input] : parameters of the radial fit.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the fit.
	*/
     virtual Term_eq harmonics_sym (const Term_eq& so, const Term_eq& omega, int bound , Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param& param, const Array<double>& passage) const ;

	/**
	* Fit, spherical harmonic by spherical harmonic, for an anti-symmetric function.
	* @param so [input] : input scalar field.
	* @param omega [input] : angular velocity.
	* @param bound [input] : the boundary at which the computation is done.
	* @param f [input] : pointer on the radial function describing the radial fit.
	* @param param [input] : parameters of the radial fit.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the fit.
	*/
     virtual Term_eq harmonics_asym (const Term_eq&, const Term_eq&, int, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param&, const Array<double>&) const ;

	/**
	* Extraction of a given multipole, at some boundary, for the radial derivative a  symmetric scalar function.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param bound [input] : the boundary at which the computation is done.
	* @param so [input] : input scalar field.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the multipolar coefficient.
	*/
     virtual Term_eq der_multipoles_sym (int k, int j, int bound, const Term_eq& so, const Array<double>& passage) const ;

	/**
	* Extraction of a given multipole, at some boundary, for the radial derivative of an anti-symmetric scalar function.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param bound [input] : the boundary at which the computation is done.
	* @param so [input] : input scalar field.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the multipolar coefficient.
	*/
     virtual Term_eq der_multipoles_asym (int k, int j, int bound, const Term_eq& so, const Array<double>& passage) const ;



 	/**
	* Gives some radial fit for a given multipole, intended for the radial derivative of an anti-symmetric scalar function.
	* @param space [input] : the concerned \c Space.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param omega [input] : angular velocity.
	* @param f [input] : pointer on the radial function describing the fit.
	* @param param [input] : parameters of the radial fit.
	* @returns the radial fit.
	*/
     virtual Term_eq der_radial_part_asym (const Space& space, int k, int j, const Term_eq& omega, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param& param) const ;

	/**
	* Gives some radial fit for a given multipole, intended for the radial derivative of a symmetric scalar function.
	* @param space [input] : the concerned \c Space.
	* @param k [input] : index of the spherical harmonic for \f$\varphi\f$.
	* @param j [input] : index of the spherical harmonic for \f$\theta\f$.
	* @param omega [input] : angular velocity.
	* @param f [input] : pointer on the radial function describing the fit.
	* @param param [input] : parameters of the radial fit.
	* @returns the radial fit.
	*/
     virtual Term_eq der_radial_part_sym (const Space& space, int k, int j, const Term_eq& omega, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param& param), const Param& param) const ;


	/**
	* Fit, spherical harmonic by spherical harmonic, for the radial derivative of a symmetric function.
	* @param so [input] : input scalar field.
	* @param omega [input] : angular velocity.
	* @param bound [input] : the boundary at which the computation is done.
	* @param f [input] : pointer on the radial function describing the radial fit.
	* @param param [input] : parameters of the radial fit.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the fit.
	*/
     virtual Term_eq der_harmonics_sym (const Term_eq& so, const Term_eq& omega, int bound, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param& param, const Array<double>& passage) const ;

	/**
	* Fit, spherical harmonic by spherical harmonic, for the radial derivative of an anti-symmetric function.
	* @param so [input] : input scalar field.
	* @param omega [input] : angular velocity.
	* @param bound [input] : the boundary at which the computation is done.
	* @param f [input] : pointer on the radial function describing the radial fit.
	* @param param [input] : parameters of the radial fit.
	* @param passage [input] : passage matrix describing the spherical harmonics.
	* @returns the fit.
	*/
     virtual Term_eq der_harmonics_asym (const Term_eq& so, const Term_eq& omega, int bound, Term_eq (*f) (const Space&, int, int, const Term_eq&, const Param&), const Param& param, const Array<double>& passage) const ;

     // Abstract stuff
  protected:
     /**
      * Function used to apply the same operation to all the components of a tensor, in the current domain.
      * It works at the \c Term_eq level and the same operation is applied to both the value and the variation of the tensor.
      * @param so : the source tensor.
      * @param pfunc : pointer onf the function to be applied.
      * @returns the resulting \c Term_eq.
      */
     Term_eq do_comp_by_comp (const Term_eq& so, Val_domain (Domain::*pfunc) (const Val_domain&) const) const ;
 	/**
      * Function used to apply the same operation to all the components of a tensor, in the current domain.
      * The operation applied depends on an integer parameter.
      * It works at the \c Term_eq level and the same operation is applied to both the value and the variation of the tensor.
      * @param so : the source tensor.
      * @param val : the integer parameter.
      * @param pfunc : pointer onf the function to be applied.
      * @returns the resulting \c Term_eq.
      */
     Term_eq do_comp_by_comp_with_int (const Term_eq& so, int val, Val_domain (Domain::*pfunc) (const Val_domain&, int) const) const ;

  public:
	/**
	* Returns the normal derivative of a \c Term_eq
	* @param so : input field.
	* @param bound : the boundary.
	* @returns  the normal derivative.
	*/
     virtual Term_eq der_normal_term_eq (const Term_eq& so, int bound) const ;

	/**
	* Returns the division by \f$1-x^2\f$ of a \c Term_eq
	* @param so : input field.
	* @returns  the result of the division.
	*/
	virtual Term_eq div_1mx2_term_eq (const Term_eq&) const ;

	/**
	* Returns the flat Laplacian of \c Term_eq, for a given harmonic.
	* @param so : input field.
	* @param m : the index of the harmonic.
	* @returns  the Laplacian.
	*/
     virtual Term_eq lap_term_eq (const Term_eq& so, int m) const ;

     	/**
	* Returns the flat 2d-Laplacian of \c Term_eq, for a given harmonic.
	* @param so : input field.
	* @param m : the index of the harmonic.
	* @returns  the 2d-Laplacian.
	*/
     virtual Term_eq lap2_term_eq (const Term_eq& so, int m) const ;

     	/**
	* Multiplication by \f$r\f$ of a \c Term_eq.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq mult_r_term_eq (const Term_eq& so) const ;

	 /**
	* Volume integral of a \c Term_eq. The volume term must be provided by the actual function.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq integ_volume_term_eq (const Term_eq& so) const ;

	 /**
	* Gradient of \c Term_eq.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq grad_term_eq (const Term_eq& so) const ;

     	/**
	* Division by \f$r\f$ of a \c Term_eq.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq div_r_term_eq (const Term_eq&) const ;


	 /**
	* Surface integral of a \c Term_eq. The surface term must be provided by the actual function.
	* @param so : input field.
	* @param bound : the surface on which the integral is performed.
	* @returns  result.
	*/
     virtual Term_eq integ_term_eq (const Term_eq& so, int bound) const ;

	/**
	* Radial derivative of a \c Term_eq.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq dr_term_eq (const Term_eq& so) const ;

	/**
	* Time derivative of a \c Term_eq.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq dtime_term_eq (const Term_eq& so) const ;

	/**
	* Second time derivative of a \c Term_eq.
	* @param so : input field.
	* @returns  result.
	*/
     virtual Term_eq ddtime_term_eq (const Term_eq& so) const ;

	 /**
	* Computes the flat derivative of a \c Term_eq, in spherical orthonormal coordinates.
	* If the index of the derivative is present in the source, appropriate contraction is performed.
	* If the contravariant version is called for, the index is raised using an arbitrary metric.
	* @param tipe : type of derivative (\t COV or \t CON)
	* @param ind : name of the index corresponding to the derivative.
	* @param so : input field.
	* @param manip : pointer on the metric used to manipulate the derivative index, if need be.
	* @returns  result.
	*/
     virtual Term_eq derive_flat_spher (int tipe, char ind, const Term_eq& so, const Metric* manip) const ;


	 /**
	* Computes the flat derivative of a \c Term_eq, in spherical coordinates
	* where the constant radii sections have a negative curvature.
	* If the index of the derivative is present in the source, appropriate contraction is performed.
	* If the contravariant version is called for, the index is raised using an arbitrary metric.
	* @param tipe : type of derivative (\t COV or \t CON)
	* @param ind : name of the index corresponding to the derivative.
	* @param so : input field.
	* @param manip : pointer on the metric used to manipulate the derivative index, if need be.
	* @returns  result.
	*/
     virtual Term_eq derive_flat_mtz (int tipe, char ind, const Term_eq& so, const Metric* manip) const ;


 	/**
	* Computes the flat derivative of a \c Term_eq, in Cartesian coordinates.
	* If the index of the derivative is present in the source, appropriate contraction is performed.
	* If the contravariant version is called for, the index is raised using an arbitrary metric.
	* @param tipe : type of derivative (\t COV or \t CON)
	* @param ind : name of the index corresponding to the derivative.
	* @param so : input field.
	* @param manip : pointer on the metric used to manipulate the derivative index, if need be.
	* @returns  result.
	*/
     virtual Term_eq derive_flat_cart (int tipe, char ind, const Term_eq& so, const Metric* manip) const ;

	/**
	* Gives the number of unknowns coming from the variable shape of the domain.
	*/
     virtual int nbr_unknowns_from_adapted () const {return 0;} ;
	/**
	* The \c Term_eq describing the variable shape of the \c Domain are updated.
	*/
     virtual void vars_to_terms() const {} ;
	/**
	* The variation of the functions describing the shape of the \c Domain are affected from the unknowns of the system.
	* @param conte : current position if the unknowns.
	* @param cc : position of the unknown to be set.
	* @param found : \c true if the index was indeed corresponding to one of the coefficients of the variable domain.
	*/
     virtual void affecte_coef(int& conte, int cc, bool& found) const {} ;
	/**
	* Computes the new boundary of a \c Domain from a set of values.
	* @param shape : \c Val_domain describing the variable boundary of the \c Domain.
	* @param xx : set of values used by the affectation
	* @param conte : current position in the values vector.
	*/
     virtual void xx_to_vars_from_adapted (Val_domain& shape, const Array<double>& xx, int& conte) const {} ;

	/**
	* Computes the new boundary of a \c Domain from a set of values.
	* @param bound : Variable boundary of the \c Domain.
	* @param xx : set of values used by the affectation
	* @param conte : current position in the values vector.
	*/
     virtual void xx_to_vars_from_adapted (double bound, const Array<double>& xx, int& conte) const {} ;

	/**
	* Affects the derivative part of variable a \c Domain from a set of values.
	* @param xx : set of values used by the affectation
	* @param conte : current position in the values vector.
	*/
     virtual void xx_to_ders_from_adapted (const Array<double>& xx, int& conte) const {} ;

	/**
	* Update the value of a field, after the shape  of the \c Domain has been changed by the system.
	* @param so : pointer on the \c Term_eq to be updated.
	*/
     virtual void update_term_eq (Term_eq* so) const  ;
	/**
	* Update the value of a scalar, after the shape  of the \c Domain has been changed by the system.
	* This is intended for variable fields and the new valued is computed using a first order Taylor expansion.
	* @param shape : modification to the shape of the \c Domain.
	* @param oldval : old value of the scalar field.
	* @param newval : new value of the scalar field.
	*/
     virtual void update_variable (const Val_domain& shape, const Scalar& oldval, Scalar& newval) const {} ;


	/**
	* Update the value of a scalar, after the shape  of the \c Domain has been changed by the system.
	* This is intended for variable fields and the new valued is computed using a first order Taylor expansion.
	* @param bound : modification of the bound of the \c Domain.
	* @param oldval : old value of the scalar field.
	* @param newval : new value of the scalar field.
	*/
     virtual void update_variable (double bound, const Scalar& oldval, Scalar& newval) const {} ;

	/**
	* Update the value of a scalar, after the shape  of the \c Domain has been changed by the system.
	* This is intended for constant field and the new valued is computed using a true spectral summation.
	* @param shape : modification to the shape of the \c Domain.
	* @param oldval : old value of the scalar field.
	* @param newval : new value of the scalar field.
	*/
     virtual void update_constante (const Val_domain& shape, const Scalar& oldval, Scalar& newval) const {} ;

	/**
	* Update the value of a scalar, after the shape  of the \c Domain has been changed by the system.
	* This is intended for constant field and the new valued is computed using a true spectral summation.
	* @param bound : modification to the bound of the \c Domain.
	* @param oldval : old value of the scalar field.
	* @param newval : new value of the scalar field.
	*/
     virtual void update_constante (double bound, const Scalar& oldval, Scalar& newval) const {} ;


	/**
	* Updates the variables parts of the \c Domain
	* @param shape : correction to the variable boundary.
	*/
     virtual void update_mapping(const Val_domain& shape) {} ;
     /**
	* Updates the variables parts of the \c Domain
	* @param  double : correction to the variable boundary.
	*/
     virtual void update_mapping(double bound) {} ;

     /**
      * Sets the value at infinity of a \c Val_domain : not implemented for this type of \c Domain.
      * @param so [input/output] : the \c Val_domain.
      * @param xx [input] : value at infinity.
      */
     virtual void set_val_inf (Val_domain& so, double xx) const ;
	/**
	* Multiplication by \f$ r\f$
	* @param so : the input
	* @returns  the output
	*/
     virtual Val_domain mult_r (const Val_domain& so) const ;
	/**
	* Multiplication by \f$ x\f$
	* @param so : the input
	* @returns  the output
	*/
     virtual Val_domain mult_x (const Val_domain& so) const ;
	/**
	* Division by \f$ r\f$
	* @param so : the input
	* @returns  the output
	*/
     virtual Val_domain div_r (const Val_domain& so) const ;
	/**
	* Division by \f$ 1 - r/L\f$
	* @param so : the input
	* @returns  the output
	*/
     virtual Val_domain div_1mrsL (const Val_domain& so) const ;
	/**
	* Multiplication by \f$ 1 - r/L\f$
	* @param so : the input
	* @returns  the output
	*/
     virtual Val_domain mult_1mrsL (const Val_domain& so) const ;

	/**
	* Volume integral. The volume element is provided by the function.
	* @param so : the input scalar field.
	* @returns  the integral.
	*/
     virtual double integ_volume (const Val_domain&) const ;

	/**
	* Gives the informations corresponding the a touching neighboring domain.
	* @param dom : index of the currect domain (deprecated, should be the same as num_dom)
	* @param bound : the boundary between the two domains.
	* @param otherdom : index of the other domain.
	* @param otherbound : name of the boundary, as seen by the other domain.
	*/
     virtual void find_other_dom (int dom, int bound, int& otherdom, int& otherbound) const ;
	/**
	* Normal derivative with respect to a given surface.
	* @param so : the input scalar field.
	* @param bound : boundary at which the normal derivative is computed.
	* @returns  the normal derivative.
	*/
     virtual Val_domain der_normal (const Val_domain& so, int bound) const ;
	/**
	* Partial derivative with respect to a coordinate.
	* @param so : the input scalar field.
	* @param ind : index of the variable used by the derivative.
	* @returns  the partial derivative.
	*/
     virtual Val_domain der_partial_var (const Val_domain& so, int ind) const ;

	/**
	* Surface integral on a given boundary.
	* @param so : the input scalar field.
	* @param bound : boundary at which the integral is computed.
	* @returns the surface integral.
	*/
     virtual double integ (const Val_domain&, int) const ;
	/**
	* Volume integral. The volume element is provided by the function (need some cleaning : same as integ_volume)
	* @param so : the input scalar field.
	* @returns the integral.
	*/
     virtual double integrale (const Val_domain&) const ;

	/**
	* Computes the part of the gradient containing the partial derivative of the field, in spherical orthonormal coordinates.
	* @param so : the input \c Term_eq
	* @returns the result, being \f$(\partial_r, \partial_\theta /r , \partial_\varphi/r/sin\theta)\f$
	*/
     virtual Term_eq partial_spher (const Term_eq& so) const ;
	  /**
	* Computes the part of the gradient containing the partial derivative of the field, in Cartesian coordinates.
	* @param so : the input \c Term_eq
	* @returns the result.
	*/
     virtual Term_eq partial_cart (const Term_eq& so) const ;
	/**
	* Computes the part of the gradient containing the partial derivative of the field, in  orthonormal coordinates where
	* the constant radius sections have negative curvature.
	* @param so : the input \c Term_eq
	* @returns the result, being \f$(\partial_r, \frac{\cos \theta}{r} \partial_\theta , \frac{\cos \theta}{r \sin \theta} \partial_\varphi)\f$
	*/
     virtual Term_eq partial_mtz (const Term_eq& so) const ;
	/**
	* Computes the part of the gradient involving the connections, in spherical orthonormal coordinates.
	* @param so : the input \c Term_eq
	* @returns the result.
	*/
     virtual Term_eq connection_spher (const Term_eq& so) const ;
     /**
	* Computes the part of the gradient involving the connections, in spherical coordinates
	* where the constant radius sections have negative curvature.
	* @param so : the input \c Term_eq
	* @returns the result.
	*/
     virtual Term_eq connection_mtz (const Term_eq& so) const ;

	/**
	* Computes the value of a field at a boundary. The result correspond to one particular coefficient.
	* @param bound : name of the boundary at which the result is computed.
	* @param so : input scalar field.
	* @param ind : indexes describing which coefficient is computed. The index corresponding the surface is irrelevant.
	* @returns one of the coefficients of the field, on the surface (given by ind).
	*/
     virtual double val_boundary (int bound, const Val_domain& so, const Index& ind) const ;
	/**
	* Computes the number of relevant collocation points on a boundary.
	* @param bound : the boundary.
	* @param base : the spectral basis. Its symmetries are used to get the right result.
	* @returns the true number of degrees of freedom on the boundary.
	*/
     virtual int nbr_points_boundary (int bound, const Base_spectral& base) const ;
	/**
	* Lists all the indices corresponding to true collocation points on a boundary.
	* @param bound: the boundary.
	* @param base : the spect basis. Its symmetries are used to get the right result.
	* @param ind : the list of indices.
	* @param start : states where the indices are stored in ind. The first is in ind[start]
	*/
     virtual void do_which_points_boundary (int bound, const Base_spectral& base, Index** ind, int start) const ;

	/**
	* Computes the number of true unknowns of a \c Tensor, in a given domain.
	* It takes into account the various symmetries and regularity conditions to determine the precise number of degrees of freedom.
	* @param so : the tensorial field.
	* @param dom : the domain considered (should be the same as num_dom).
	* @returns the number of true unknowns.
	*/
     virtual int nbr_unknowns (const Tensor& so, int dom) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	* @returns the number of true conditions, component by component.
	*/
     virtual Array<int> nbr_conditions (const Tensor& eq, int dom, int order, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	* @returns the number of true conditions, component by component.
	*/
     virtual Array<int> nbr_conditions_boundary (const Tensor& eq, int dom, int bound, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports all the residual equations corresponding to a tensorial one in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param order : describes the order of the equation (2 for a Laplacian for instance).
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : \c Array containing the number of equations corresponding to each component. It is used when some of the components are null.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	*/
     virtual void export_tau (const Tensor& eq, int dom, int order, Array<double>& res, int& pos_res, const Array<int>& ncond,int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports all the residual equations corresponding to a tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : \c Array containing the number of equations corresponding to each component. It is used when some of the components are null.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	*/
     virtual void export_tau_boundary (const Tensor& eq , int dom, int bound, Array<double>& res, int& pos_res, const Array<int>& ncond, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports all the residual equations corresponding to one tensorial one on a given boundary, excepted for some coefficients where another equation is used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : \c Array containing the number of equations corresponding to each component. It is used when some of the components are null.
	* @param param : parameters describing the coefficients where the alternative condition is enforced.
	* @param type_exception : states which type of exception (value or derivative ; current domain or the other one). Highly specialized...
	* @param exception : the equation used for the alternative condition.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	* @returns the number of true conditions.
	*/
     virtual void export_tau_boundary_exception (const Tensor& eq, int dom, int bound, Array<double>& res, int& pos_res, const Array<int>& ncond,
			const Param& param, int type_exception, const Tensor& exception, int n_cmp=-1,  Array<int>** p_cmp=0x0) const ;

	/**
	* Affects some coefficients to a \c Tensor.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Tensor to be affected.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param cf : \c Array of the coefficients used.
	* @param pos_cf : current position in the array of coefficients.
	*/
     virtual void affecte_tau (Tensor& so, int dom, const Array<double>& cf, int& pos_cf) const ;

	/**
	* Sets at most one coefficient of a \c Tensor to 1.
	* It takes into account the various symmetries and regularity conditions (by means of Garlekin basis).
	* @param so : the \c Tensor to be affected. It is set to zero if cc does not corresponds to another field.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param cc : location, in the overall system, of the coefficient to be set to 1.
	* @param pos_cf : current position.
	*/
     virtual void affecte_tau_one_coef (Tensor& so, int dom, int cc, int& pos_cf) const ;

	/**
	* Computes number of discretized equations associated with a given tensorial equation in the bulk.
	* It takes into account the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param order : describes the order of the equation, with respect to each variable.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	* @returns the number of true conditions, component by component.
	*/
     virtual Array<int> nbr_conditions_array (const Tensor& eq, int dom, const Array<int>& order, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Computes number of discretized equations associated with a given tensorial equation on a boundary.
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param order : describes the order of the equation, with respect to each variable. The one normal to the surface is irrelevant.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	* @returns the number of true conditions, component by component.
	*/
     virtual Array<int> nbr_conditions_boundary_array (const Tensor& eq, int dom, int bound, const Array<int>& order, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports all the residual equations corresponding to one tensorial one in the bulk.
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param order : describes the order of the equation, with respect to each variable.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : \c Array containing the number of equations corresponding to each component. It is used when some of the components are null.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	*/
     virtual void export_tau_array (const Tensor& eq, int dom, const Array<int>& order, Array<double>& res, int& pos_res, const Array<int>& ncond,
														int n_cmp=-1, Array<int>** p_cmp=0x0) const ;
	/**
	* Exports all the residual equations corresponding to one tensorial one on a given boundary
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param order : describes the order of the equation, with respect to each variable. The one normal to the surface is irrelevant.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : \c Array containing the number of equations corresponding to each component. It is used when some of the components are null.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	*/
     virtual void export_tau_boundary_array (const Tensor& eq, int dom, int bound, const Array<int>& order, Array<double>& res, int& pos_res,
				const Array<int>& ncond,  int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Computes number of discretized equations associated with a given tensorial equation on a boundary.
	* The boundary is assumed to also have boundaries. (Used for bispherical coordinates).
	* It takes into account the various Galerkin basis used.
	* It is used for implementing boundary conditions and matching ones.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	* @returns the number of true conditions, component by component.
	*/
    virtual Array<int> nbr_conditions_boundary_one_side (const Tensor& eq, int dom, int bound, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

	/**
	* Exports all the residual equations corresponding to one tensorial one on a given boundary.
	* The boundary is assumed to also have boundaries. (Used for bispherical coordinates).
	* It makes use of the various Galerkin basis used.
	* @param eq : the residual of the equation.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary at which the equation is enforced.
	* @param res : The \c Array where the discretized equations are stored.
	* @param pos_res : current position in res.
	* @param ncond : \c Array containing the number of equations corresponding to each component. It is used when some of the components are null.
	* @param n_cmp : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param p_cmp : pointer on the indexes of the components to be considered. Not used of n_cmp = -1 .
	*/
     virtual void export_tau_boundary_one_side (const Tensor& eq, int dom, int bound, Array<double>& res, int& pos_res,
				const Array<int>& ncond, int n_cmp=-1, Array<int>** p_cmp=0x0) const ;

      /**
	* Translates a name of a coordinate into its corresponding numerical name
	* @param name : name of the variable (like 'R', for the radius).
	* @returns the corresponding number (1 for 'R' for instance). The result is -1 if the name is not recognized.
	*/
     virtual int give_place_var (char* name) const ;

	/**
	* Gets the value of a \c Term_eq by importing data from neighboring domains, on a boundary.
	* This makes use of the spectral summation.
	* @param numdom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary where the values from the other domains are computed.
	* @param n_ope : number of other domains touching the current one at the given boundary.
	* @param parts : pointers on the various \c Term_eq that have to be imported.
	* @returns the imported field. Assumed to be constant in the direction normal to the boundary.
	*/

     Term_eq import (int numdom, int bound, int n_ope, Term_eq** parts) const ;

	/**
	* Gets the value of a \c Tensor by importing data from neighboring domains, on a boundary.
	* This makes use of the spectral summation.
	* @param numdom : the domain considered (should be the same as num_dom).
	* @param bound : the boundary where the values from the other domains are computed.
	* @param n_ope : number of other domains touching the current one at the given boundary.
	* @param other_doms : numbers of the other pertinent domains.
	* @param parts : pointers on the various \c Tensors that have to be imported.
	* @returns the imported field. Assumed to be constant in the direction normal to the boundary.
	*/
     virtual Tensor import (int numdom, int bound, int n_ope, const Array<int>& other_doms, Tensor** parts) const ;

	/**
	* Puts to zero all the coefficients below a given treshold.
	* Regularity is maintained by dealing properly with the various Galerkin basis
	* @param tt : the \c Tensor to be filtered.
	* @param dom : the domain considered (should be the same as num_dom).
	* @param treshold : all coefficients below this are set to zero.
	*/
     virtual void filter (Tensor& tt, int dom, double treshold) const ;

public:
	/**
	* Print function, to allow override in operator<<.
	* @param out : the output stream.
	*/
     virtual ostream& print (ostream& o) const = 0 ;

     friend ostream& operator<< (ostream& o, const Domain& so) {return so.print(o);}; ///< Display
     friend class Val_domain ;
     friend class Metric_ADS  ;
     friend class Metric_AADS ;
} ;

/**
 * The \c Space class is an ensemble of domains describing the whole space of the computation.
 * It is a purely abstract class and so space must be constructed via the derived classes.
 * \ingroup domain
 */
class Space : public MemoryMappable {

  protected:
    int nbr_domains ; ///< Number od \c Domains.
    int ndim ; ///< Number of dimensions (should be the same for all the \c Domains).
    int type_base ; ///< Type of basis used (i.e. using either Chebyshev or Legendre polynomials).
    Domain** domains ; ///< Pointers on the various \c Domains.
    Space () ; ///< Standard constructor
    virtual ~Space () ; ///< Destructor

  public:
     int get_ndim() const /// Returns the number of dimensions
	{return ndim ;} ;
     int get_nbr_domains() const /// Returns the number of \c Domains
	{return nbr_domains ;} ;
     int get_type_base() const /// Returns the type of basis
	{return type_base ;} ;
     virtual void save (FILE*) const ; ///< Saving function

    /**
     * returns a pointer on the domain.
     * @param i [input] : the index of the domain.
     */
    const Domain* get_domain(int i) const {
    	assert ((i>=0) && (i<nbr_domains)) ;
	return domains[i] ;
    }

		Scalar get_cart_field(int const cart) const;
		Scalar get_cart_field_bound(int const cart, int const bound) const;

    // Things for adapted domains
	/**
	* Gives the number of unknowns coming from the variable shape of the domain.
	*/
    virtual int nbr_unknowns_from_variable_domains() const {return 0;} ;
	/**
	* The variation of the functions describing the shape of the \c Domain are affected from the unknowns of the system.
	* @param conte : current position if the unknowns.
	* @param cc : position of the unknown to be set.
	* @param doms : Array containing the index of the variable domains present in this space.
	*/
    virtual void affecte_coef_to_variable_domains (int&conte, int cc, Array<int>& doms) const {};
	/**
	* Update the vairable domains from a set of values.
	* @param xx : set of values used by the affectation
	* @param conte : current position in the values vector.
	*/
    virtual void xx_to_ders_variable_domains (const Array<double>& xx, int& conte) const {};
	/**
	* Update the variables of a system, from the variation of the shape of the domains.
	* @param syst : the \c System_of_eqs considered.
	* @param xx : set of values used by the affectation
	* @param conte : current position in the values vector.
	*/
    virtual void xx_to_vars_variable_domains (System_of_eqs* syst, const Array<double>& xx, int& pos) const {};

	/**
	* Gives the number of the other domains, touching a given boundary.
	* It also gives the name of the boundary, as seen by the other domain.
	* @param dom : the domain considered.
	* @param bound : the boundary.
	* @returns a 2d-array containing the numbers of the other domains (stored in (0,i)) and the names of the boundary (stored in (1,i)).
	*/
    virtual Array<int> get_indices_matching_non_std(int dom, int bound) const ;

  friend ostream& operator << (ostream& o, const Space& so) ; ///< Display
} ;
}
#endif
