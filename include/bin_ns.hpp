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

#ifndef __BIN_NS_HPP_
#define __BIN_NS_HPP_

#include "space.hpp"
#include "spheric.hpp"
#include "adapted.hpp"
#include "list_comp.hpp"
#include "bispheric.hpp"

namespace Kadath {
/**
 * Spacetime intended for binary neutron stars configurations (see constructor for details about the domains used)
 * \ingroup domain
 */

class Space_bin_ns : public Space {

     protected :
  	int n_shells_outer{0} ; ///< Number of outer shells.
    int n_shells1{0}; ///< Number of shells around the first nucleus.
    int n_shells2{0}; ///< Number of shells around the second nucleus.

     public:
    int NS1{0}; ///< Starting index of the first spheres.
    int NS2{3}; ///< Starting index of the second spheres.
    int ADAPTED1{1}; ///< Starting index of the first adapted domains.
    int ADAPTED2{4}; ///< Starting index of the second adapted domains.
    int OUTER{6}; ///< Starting index of the outer domains.

    /**
     * Constructor with an outer shell surrounding the bispheric region; stars are initial spherical.
     * @param ttype [input] : the type of basis.
     * @param dist [input] : distance \f$ d \f$ between the centers of the two spheres.
     * @param rinstar1 [input] : radius  the first nucleus (inside the star).
     * @param rstar1 [input] : radius  the first shell with variable outer surface (inside the star).
     * @param routstar1 [input] : radius  the first shell with variable inner surface (outide the star).
     * @param rinstar2 [input] : radius  the second nucleus (inside the star).
     * @param rstar2 [input] : radius  the second shell with variable outer surface (inside the star).
     * @param routstar2 [input] : radius  the second shell with variable inner surface (outide the star).
     * @param rext [input] : outer radius the bispherical part.
     * @param rshell [input] : outer radius the outer shell.
     * @param nr [input] : number of points in each dimension
     *
     * The various domains are then :
     * \li One \c Domain_nucleus, of radius  rinstar1, centered at \f$ x_1 \f$.
     * \li One \c Domain_shell_outer_adapted  centered at \f$ x_1 \f$.
     * \li One \c Domain_shell_inner_homothetic  centered at \f$ x_1 \f$.
     * \li One \c Domain_nucleus, of radius  rinstar2, centered at \f$ x_2 \f$.
     * \li One \c Domain_shell_outer_adapted  centered at \f$ x_2 \f$.
     * \li One \c Domain_shell_inner_homothetic  centered at \f$ x_2 \f$.
     * \li One \c Domain_bispheric_chi_first near the first sphere.
     * \li One \c Domain_bispheric_rect near the first sphere.
     * \li One \c Domain_bispheric_eta_first inbetween the two spheres.
     * \li One \c Domain_bispheric_rect near the second sphere.
     * \li One \c Domain_bispheric_chi_first near the second sphere.
     * \li One \c Domain_shell centered on the origin
     * \li One \c Domain_compact centered on the origin, with inner radius \f$ R \f$.
     */
	Space_bin_ns (int ttype, double dist, double rinstar1, double rstar1, double routstar1,
			double rinstar2, double rstar2, double routstar2, double rext, double rshell, int nr) ;

    /**
     * Standard constructor ; stars are initial spherical.
     * @param ttype [input] : the type of basis.
     * @param dist [input] : distance \f$ d \f$ between the centers of the two spheres.
     * @param rinstar1 [input] : radius  the first nucleus (inside the star).
     * @param rstar1 [input] : radius  the first shell with variable outer surface (inside the star).
     * @param routstar1 [input] : radius  the first shell with variable inner surface (outide the star).
     * @param rinstar2 [input] : radius  the second nucleus (inside the star).
     * @param rstar2 [input] : radius  the second shell with variable outer surface (inside the star).
     * @param routstar2 [input] : radius  the second shell with variable inner surface (outide the star).
     * @param rext [input] : radius \f$ R \f$ of the outer boundary of the bispherical part.
     * @param nr [input] : number of points in each dimension
     * @param nshells1 [input] : number of additional shells around the first nucleus
     * @param nshells2 [input] : number of additional shells around the second nucleus
     *
     * The various domains are then :
     * \li One \c Domain_nucleus, of radius  rinstar1, centered at \f$ x_1 \f$.
     * \li \f$ nshells1 \f$ \c Domain_shell centered at \f$ x_1 \f$.
     * \li One \c Domain_shell_outer_adapted  centered at \f$ x_1 \f$.
     * \li One \c Domain_shell_inner_homothetic  centered at \f$ x_1 \f$.
     * \li One \c Domain_nucleus, of radius  rinstar2, centered at \f$ x_2 \f$.
     * \li \f$ nshells2 \f$ \c Domain_shell centered at \f$ x_2 \f$.
     * \li One \c Domain_shell_outer_adapted  centered at \f$ x_2 \f$.
     * \li One \c Domain_shell_inner_homothetic  centered at \f$ x_2 \f$.
     * \li One \c Domain_bispheric_chi_first near the first sphere.
     * \li One \c Domain_bispheric_rect near the first sphere.
     * \li One \c Domain_bispheric_eta_first inbetween the two spheres.
     * \li One \c Domain_bispheric_rect near the second sphere.
     * \li One \c Domain_bispheric_chi_first near the second sphere.
     * \li One \c Domain_compact centered on the origin, with inner radius \f$ R \f$.
     */
	Space_bin_ns (int ttype, double dist, double rinstar1, double rstar1, double routstar1,
			double rinstar2, double rstar2, double routstar2, double rext, int nr, int nshells1 = 0, int nshells2 = 0) ;

  Space_bin_ns (int ttype, double dist, const std::vector<double>& NS1_bounds, const std::vector<double>& NS2_bounds,
                            const std::vector<double>& outer_bounds, int nr);

	Space_bin_ns (const Space_bin_ns &) ; ///< Constructor from existing space.
	Space_bin_ns (FILE*, bool old = false) ; ///< Constructor from a file for old data.

 	virtual ~Space_bin_ns() ; ///< Destructor
	virtual void save(FILE*) const ;

	/**
	* Adds a bulk equation and two matching conditions.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ;

	/**
	* Adds a bulk equation and two matching conditions between nuclii and adapted domains.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq_matter (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ;

	/**
	* Adds a bulk equation and two matching conditions. The nuclei and inner adapted domains are excluded from the computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq_excised (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ;

	/**
	* Sets a boundary condition at the inner radius of the first outer shell.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_bc_sphere_one (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
	/**
	* Sets a boundary condition at the inner radius of the second outer shell.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_bc_sphere_two (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;
	/**
	* Sets a boundary condition at the outer boundary.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the boundary condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_bc_outer (System_of_eqs& syst, const char* eq, int nused=-1, Array<int>** pused=0x0)  ;

	/**
	* Adds a bulk equation and two matching conditions. The compactified domain is excluded from the computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq_nozec (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, int nused=-1, Array<int>** pused=0x0)  ;

/**
	* Adds a bulk equation and two matching conditions. The compactified domain is excluded from the computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param used : list of components used
	*/
	void add_eq_nozec (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, const List_comp& used)  ;

	/**
	* Adds a bulk equation and two matching conditions. The outer shells are excluded from the computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param nused : number of components of \c eq to be considered. All the components are used of it is -1.
	* @param pused : pointer on the indexes of the components to be considered. Not used of nused = -1 .
	*/
	void add_eq_noshell (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der,int nused=-1, Array<int>** pused=0x0)  ;


	/**
	* Adds a bulk equation and two matching conditions. The outer shells are excluded from the computational space.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the bulk equation.
	* @param rac : the string describing the first matching condition.
	* @param rac_der : the string describing the second matching condition.
	* @param used : list of components used
	*/
	void add_eq_noshell (System_of_eqs& syst, const char* eq, const char* rac, const char* rac_der, const List_comp& used)  ;


	virtual int nbr_unknowns_from_variable_domains() const ;
	virtual void affecte_coef_to_variable_domains(int& , int, Array<int>&) const ;
	virtual void xx_to_ders_variable_domains(const Array<double>&, int&) const ;
	virtual void xx_to_vars_variable_domains(System_of_eqs*, const Array<double>&, int&) const ;
	virtual Array<int> get_indices_matching_non_std(int dom, int bound) const ;

	/**
	* Adds an equation being a surface integral at infinity.
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the equation (should contain something like integ(f)=b)
	*/
	void add_eq_int_inf (System_of_eqs& syst, const char* eq) ;
	/**
	* Adds an equation being the value of some field at the origin of the first nucleus
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the quantity that must be zero at the origin
	*/
	void add_eq_ori_one (System_of_eqs& syst, const char* eq) ;
	/**
	* Adds an equation being the value of some field at the origin of the first nucleus
	* @param syst : the \c System_of_eqs.
	* @param eq : the string describing the quantity that must be zero at the origin
	*/
	void add_eq_ori_two (System_of_eqs& syst, const char* eq) ;
	/**
	* Adds an equation being a volume integral over contiguous domains.
	* @param syst : the \c System_of_eqs.
	* @param dmin : the integral is taken beginning at this domain.
	* @param dmax : the integral is taken up to this domain.
	* @param eq : the string describing the equation (should contain something like integvolume(f)=b)
	*/
	void add_eq_int_volume (System_of_eqs& sys, int dmin, int dmax, const char* nom) ;
	/**
	* Adds an equation being a surface integral over the outer boundary of the first inner adapted domain
	* @param syst : the \c System_of_eqs.
	* @param dmin : the integral is taken beginning at this domain.
	* @param dmax : the integral is taken up to this domain.
	* @param eq : the string describing the equation (should contain something like integvolume(f)=b)
	*/
	void add_eq_int_outer_sphere_one (System_of_eqs& sys, const char* nom) ;
	/**
	* Adds an equation being a surface integral over the outer boundary of the second inner adapted domain
	* @param syst : the \c System_of_eqs.
	* @param dmin : the integral is taken beginning at this domain.
	* @param dmax : the integral is taken up to this domain.
	* @param eq : the string describing the equation (should contain something like integvolume(f)=b)
	*/
	void add_eq_int_outer_sphere_two (System_of_eqs& sys, const char* nom) ;
  
  int const & get_n_shells_outer() const { return n_shells_outer; }
} ;
}
#endif
