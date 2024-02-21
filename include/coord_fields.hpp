/*
 * This file is part of the KADATH library.
 * Copyright (C) 2019, Samuel Tootle
 *                     <tootle@th.physik.uni-frankfurt.de>
 * Copyright (C) 2019, Ludwig Jens Papenfort
 *                      <papenfort@th.physik.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once
#include "kadath.hpp"
#include <optional>

/**
 * \addtogroup domain
 * @{
 */
namespace Kadath {
enum coord_vector {GLOBAL_ROT, BCO1_ROT, BCO2_ROT, EX, EY, EZ, S_BCO1, S_BCO2, S_INF, NUM_VECTORS};
enum coord_scalar {R_BCO1=0, R_BCO2, NUM_SCALARS};
                   
using vec_ary_t = std::array<std::optional<Vector>, NUM_VECTORS>;
using scalar_ary_t = std::array<std::optional<Scalar>, NUM_SCALARS>;

/**
 * gen_cv_names()
 * in order to make sure fields are updated during each solver iteration,
 * we need to maintain a catalogue of the field names.  They must be in
 * the system of equations with the same name!
 * We generate this using a function such that the array used in the update functions
 * is const
 *
 * @return array of pre-defined names that correspond to enum coord_vector
 */
inline std::array<std::string, NUM_VECTORS> gen_cv_names() {
  std::array<std::string, NUM_VECTORS> cv_names;
  cv_names[GLOBAL_ROT] = "mg^i";
  cv_names[BCO1_ROT]   = "mm^i";
  cv_names[BCO2_ROT]   = "mp^i";
  cv_names[EX]         = "ex^i";
  cv_names[EY]         = "ey^i";
  cv_names[EZ]         = "ez^i";
  cv_names[S_BCO1]     = "sm^i";
  cv_names[S_BCO2]     = "sp^i";
  cv_names[S_INF]      = "einf^i";
  return cv_names;
}
const std::array<std::string, NUM_VECTORS> cv_names = gen_cv_names();

/**
 * gen_cs_names()
 * same as gen_cv_names only for Scalar fields
 *
 * return cs_names array of predefined names that correspond to enum coord_scalar
 */
inline std::array<std::string, NUM_SCALARS> gen_cs_names() {
  std::array<std::string, NUM_SCALARS> cs_names;
  cs_names[R_BCO1]     = "rm";
  cs_names[R_BCO2]     = "rp";
  return cs_names;
}
const std::array<std::string, NUM_SCALARS> cs_names = gen_cs_names();


// Forward declarations
template<typename space_t>
class CoordFields;

/**
 * update_fields
 *
 * this function updates the constant fields used in the generation of
 * various ID solvers.  Since the fields are a function of the cartesian basis,
 * they must be updated every iteration as the coordinates change.  However, since
 * they are a "constant" in the solver, they must be updated manually.  This includes
 * the fields themselves and their value within the system of equation. See
 * update_field()
 *
 * @tparam space_t type of computation space
 * @param [input] cf_generator: CoordFields object for a given space
 * @param [input] coord_vectors: ref to array of pointers to vector fields
 * @param [input] coord_scalars: ref to array of pointers to scalar fields
 * @param [input] xo: origin of the global fields
 * @param [input] xc1: origin of BCO1
 * @param [input] xc2: origin of BCO2
 * @param [input] syst optional pointer to the system of equations.
 */
template<typename space_t>
void update_fields (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t & coord_scalars,
                   const double xo, const double xc1, const double xc2, 
                   std::unique_ptr<System_of_eqs>& syst);

template<typename space_t>
void update_fields (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t & coord_scalars,
                   const double xo, const double xc1, const double xc2, 
                   System_of_eqs* syst=nullptr);

// Helper function to avoid declaring unnecessary coord_scalar arrays
template<typename space_t>
void update_fields (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t&& coord_scalars,
                   const double xo, const double xc1, const double xc2, 
                   System_of_eqs* syst=nullptr);
/**
 * update_fields_co
 *
 * this function is an abbreviated version of updates_fields
 * in order to update only the fields required for an isolated compact object
 *
 * @tparam space_t type of computation space
 * @param [input] cf_generator: CoordFields object for a given space
 * @param [input] coord_vectors: ref to array of coordinate vector field pointers
 * @param [input] coord_scalars: pointer to the needed radius scalar field
 * @param [input] xo: origin of the global fields
 * @param [input] syst: optional pointer to the system of equations.
 */
template<typename space_t>
void update_fields_co(CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t & coord_scalars, 
                   const double xo, System_of_eqs* syst=nullptr);

// Helper function to avoid declaring unnecessary coord_scalar arrays
template<typename space_t>
void update_fields_co(CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t && coord_scalars, 
                   const double xo, System_of_eqs* syst=nullptr);

/**
 * update_field
 *
 * this function updates a field in the given system of equations
 * based on the statically definded names: see cv_names and cs_names
 *
 * @param [input] syst: optional pointer to the system of equations.
 * @param [input] dom: dom of the field to be updated
 * @param [input] so the value to set 'n' to in the system of equations
 * @param [input] so: value to apply
 */
bool update_field(System_of_eqs& syst, int dom, const char* n, Tensor& so);

/**
 * default_binary_vector_ary
 *
 * generates the default array of Vector fields needed for XCTS
 * binary initial data.
 *
 * @tparam space_t type of computation space
 * @param [input] space: numerical spoace
 * @return array of Vector fields
 */
template<class space_t>
vec_ary_t default_binary_vector_ary(space_t& space);

/**
 * default_co_vector_ary
 *
 * generates the default array of Vector fields needed for XCTS
 * isolated compact object initial data.
 *
 * @tparam space_t type of computation space
 * @param [input] space: numerical spoace
 * @return array of Vector fields
 */
template<class space_t>
vec_ary_t default_co_vector_ary(space_t& space);

/**
 * @class CoordFields
 *
 * generates and updates necessary coordinate fields
 * @tparam space_t Space type (e.g. Space_bin_bh)
 */
template<typename space_t>
class CoordFields {
	private:
		space_t const & space; ///< store reference to computational space
		Kadath::Base_tensor basis; ///< store basis

	public:
		CoordFields(space_t const & space) : space(space), basis (space, CARTESIAN_BASIS) 
    { }

    /**
     * CoordFields::cart
     *
     * generate and return vector field of the cartesian coordinates
     *
     * @param [input] shift_x: coordinate shift in the x direction
     * @param [input] shift_y: coordinate shift in the y direction
     * @param [input] shift_z: coordinate shift in the z direction
     * @return Vector field containing the shifted cartesian coordinates
     */
    template<int ind_t = CON>
		Kadath::Vector cart(double shift_x = 0., double shift_y = 0., double shift_z= 0.) const;
    
    /**
     * CoordFields::radius
     *
     * generate and return a shifted radius field
     *
     * @param [input] shift_x: coordinate shift in the x direction
     * @param [input] shift_y: coordinate shift in the y direction
     * @param [input] shift_z: coordinate shift in the z direction
     * @return Scalar field containing the shifted radius field
     */
		Kadath::Scalar radius(double shift_x = 0., double shift_y = 0., double shift_z= 0.) const;
    
    /**
     * CoordFields::rot_z
     *
     * generate and return a shifted rotation field about the z axis
     *
     * @param [input] shift_x: coordinate shift in the x direction
     * @param [input] shift_y: coordinate shift in the y direction
     * @param [input] shift_z: coordinate shift in the z direction
     * @return Vector field containing the shifted rotation field
     */
    template<int ind_t = CON>
		Kadath::Vector rot_z(double shift_x = 0., double shift_y = 0., double shift_z= 0.) const;
    
    /**
     * CoordFields::e_rad
     *
     * generate and return a shifted, radial pointing, unit vector field
     *
     * @param [input] shift_x: coordinate shift in the x direction
     * @param [input] shift_y: coordinate shift in the y direction
     * @param [input] shift_z: coordinate shift in the z direction
     * @return Vector field containing the shifted unit vector field
     */
    template<int ind_t = CON>
		Kadath::Vector e_rad(double shift_x = 0., double shift_y = 0., double shift_z= 0.) const;
    
    /**
     * CoordFields::e_cart
     *
     * generate and return a unit vector field for a given cartesian coordinate
     *
     * @param [input] dir: coordinate direction (e.g. x = 1, y = 2, z = 3)
     * @return Vector field containing the vector field of the specified coordinate direction
     */
    template<int ind_t = CON>
		Kadath::Vector e_cart(int dir) const;
};
/**
 * @}
 */

template<typename space_t>
template<int ind_t>
Kadath::Vector CoordFields<space_t>::cart(double shift_x, double shift_y, double shift_z) const {
	Kadath::Vector cart(space, ind_t, basis);

  int ndom = space.get_nbr_domains();

  for (int d = 0; d < ndom - 1; d++) {
    cart.set(1).set_domain(d) = space.get_domain(d)->get_cart(1) - shift_x;
    cart.set(2).set_domain(d) = space.get_domain(d)->get_cart(2) - shift_y;
    cart.set(3).set_domain(d) = space.get_domain(d)->get_cart(3) - shift_z;
  }

  // set outter boundary to surface coordinates (1/r) for compact domain
  // FIXME: these are not coordinate shifted
  cart.set(1).set_domain(ndom - 1) = space.get_domain(ndom - 1)->get_cart_surr(1);
  cart.set(2).set_domain(ndom - 1) = space.get_domain(ndom - 1)->get_cart_surr(2);
  cart.set(3).set_domain(ndom - 1) = space.get_domain(ndom - 1)->get_cart_surr(3);

  cart.set(1).std_base();
  cart.set(2).std_base();
  cart.set(3).std_anti_base();

	return cart;
}

template<typename space_t>
Kadath::Scalar CoordFields<space_t>::radius(double shift_x, double shift_y, double shift_z) const {
	auto coords = this->cart(shift_x, shift_y, shift_z);
  int ndom = space.get_nbr_domains();
  
	Kadath::Scalar r_sq = coords(1) * coords(1) + coords(2) * coords(2) + coords(3) * coords(3);
  r_sq.std_base();

  Kadath::Scalar r = sqrt(r_sq);
  r.std_base();
  
  // set outter boundary to constant large radius 
  Index pos(space.get_domain(ndom-1)->get_nbr_points());
  const int npts_r = space.get_domain(ndom-1)->get_nbr_points()(0);
  for(int i = 1; i <= 3; ++i) {
    const int dom = ndom - 1;
    do{
      if(pos(0) == npts_r - 1) {
        r.set_domain(dom).set(pos) = 1e10;
      }
    }while(pos.inc());
  }

  return r;
}

template<typename space_t>
template<int ind_t>
Kadath::Vector CoordFields<space_t>::rot_z(double shift_x, double shift_y, double shift_z) const {
	auto coords = this->cart(shift_x, shift_y, shift_z);

	Kadath::Vector rot_z(space, ind_t, basis);
  rot_z.set(1) = -coords(2);
  rot_z.set(2) = coords(1);
  rot_z.set(3) = 0.;

  rot_z.std_base();
	return rot_z;
}

template<typename space_t>
template<int ind_t>
Kadath::Vector CoordFields<space_t>::e_rad(double shift_x, double shift_y, double shift_z) const {
	auto coords = this->cart(shift_x, shift_y, shift_z);

  int ndom = space.get_nbr_domains();

  // this can be problematic around the given origin
	Kadath::Vector e_rad(space, ind_t, basis);
	for(int i : {1,2,3}) {
	  e_rad.set(i) = coords(i) / this->radius(shift_x, shift_y, shift_z);
    
    // fix non-finite values likely at (0,0,0)
    // Note: this fix is sufficient for the FUKA codes as
    // we are only concerned with e_rad on the surfaces of compact
    // objects or at infinity.  This fix is to remove NaNs and their
    // impact on the system of equations
    for(int dom=0; dom < ndom - 1; ++dom) {
      Index pos(space.get_domain(dom)->get_nbr_points());
      __attribute__((unused)) const int npts_r = space.get_domain(dom)->get_nbr_points()(0);
      do{
        if(!std::isfinite(e_rad(i)(dom)(pos)) && pos(0) == 0) {
          Index tempos(pos);
          tempos.set(0) = pos(0) + 1;
          e_rad.set(i).set_domain(dom).set(pos) = e_rad(i)(dom)(tempos);
        }
      }while(pos.inc());
    }

		// FIXME this is not correct for shifted cartesian coords
    e_rad.set(i).set_domain(ndom - 1) = space.get_domain(ndom - 1)->get_cart_surr(i);
	}

  e_rad.set(1).std_base();
  e_rad.set(2).std_base();
  e_rad.set(3).std_anti_base();

	return e_rad;
}

template<typename space_t>
template<int ind_t>
Kadath::Vector CoordFields<space_t>::e_cart(int dir) const {
	Kadath::Vector e_cart(space, ind_t, basis);
	e_cart = 0;
	e_cart.set(dir) = 1.;

  // all components should act like simple scalar fields,
  // so decompose them separately
  e_cart.set(1).std_base();
  e_cart.set(2).std_base();
  e_cart.set(3).std_base();

	return e_cart;
}

template<typename space_t>
void update_fields (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t & coord_scalars,
                   const double xo, const double xc1, const double xc2, 
                   std::unique_ptr<System_of_eqs>& syst) {

  if(coord_vectors[GLOBAL_ROT]) 
    *coord_vectors[GLOBAL_ROT] = cf_generator.template rot_z(xo);
	if(coord_vectors[BCO1_ROT])   
    *coord_vectors[BCO1_ROT]   = cf_generator.template rot_z(xc1) ;
	if(coord_vectors[BCO2_ROT])   
    *coord_vectors[BCO2_ROT]   = cf_generator.template rot_z(xc2);
	if(coord_vectors[EX])
    *coord_vectors[EX]         = cf_generator.template e_cart<COV>(1);
	if(coord_vectors[EY])         
    *coord_vectors[EY]         = cf_generator.template e_cart<COV>(2);
	if(coord_vectors[EZ])         
    *coord_vectors[EZ]         = cf_generator.template e_cart<COV>(3);
  if(coord_vectors[S_BCO1])     
    *coord_vectors[S_BCO1]     = cf_generator.template e_rad<COV>(xc1);
  if(coord_vectors[S_BCO2])     
    *coord_vectors[S_BCO2]     = cf_generator.template e_rad<COV>(xc2);
  if(coord_vectors[S_INF])      
    *coord_vectors[S_INF]      = cf_generator.template e_rad<COV>(xo);
  if(coord_scalars[R_BCO1]) 
    *coord_scalars[R_BCO1]     = cf_generator.radius(xc1);
  if(coord_scalars[R_BCO2])     
    *coord_scalars[R_BCO2]     = cf_generator.radius(xc2);
  
  int ndom = coord_vectors[GLOBAL_ROT]->get_space().get_nbr_domains();
  
  auto update = [&] (auto& name, auto& field) {
    for(int dom = 0; dom < ndom; ++dom){
      __attribute__((unused)) bool succ = update_field(*syst, dom, name.c_str(), field);
      #ifdef DEBUG
      if(!succ)
        std::cout << name << " failed\n";
      #endif
    }
  };
  if(syst) {
    for(int i = 0; i < NUM_VECTORS; ++i) {
      if(coord_vectors[i])
        update(cv_names[i], *coord_vectors[i]);
    }
    for(int i = 0; i < NUM_SCALARS; ++i) {
      if(coord_scalars[i])
        update(cs_names[i], *coord_scalars[i]);
    }
    syst->sec_member();
  }
}

template<typename space_t>
void update_fields (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t& coord_scalars,
                   const double xo, const double xc1, const double xc2, 
                   System_of_eqs* syst_) {
  auto syst = std::make_unique<System_of_eqs>(*syst_);
  update_fields(cf_generator, coord_vectors, coord_scalars, xo, xc1, xc2, syst);
}

template<typename space_t>
void update_fields (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t&& coord_scalars,
                   const double xo, const double xc1, const double xc2, 
                   System_of_eqs* syst_) {
  auto syst = std::make_unique<System_of_eqs>(*syst_);
  update_fields(cf_generator, coord_vectors, coord_scalars, xo, xc1, xc2, syst);
}

template<typename space_t>
void update_fields_co (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t && coord_scalars,
                   const double xo, System_of_eqs* syst_) {
  auto syst = std::make_unique<System_of_eqs>(*syst_);
  update_fields(cf_generator, coord_vectors, coord_scalars, xo, xo, 0., syst);
}

template<typename space_t>
void update_fields_co (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t & coord_scalars,
                   const double xo, System_of_eqs* syst_) {
  auto syst = std::make_unique<System_of_eqs>(*syst_);
  update_fields(cf_generator, coord_vectors, coord_scalars, xo, xo, 0., syst);
}

template<typename space_t>
void update_fields_co (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t && coord_scalars,
                   const double xo, std::unique_ptr<System_of_eqs>& syst) {
  update_fields(cf_generator, coord_vectors, coord_scalars, xo, xo, 0., syst);
}

template<typename space_t>
void update_fields_co (CoordFields<space_t> const & cf_generator,
                   vec_ary_t & coord_vectors,
                   scalar_ary_t & coord_scalars,
                   const double xo, std::unique_ptr<System_of_eqs>& syst) {
  update_fields(cf_generator, coord_vectors, coord_scalars, xo, xo, 0., syst);
}

inline bool update_field(System_of_eqs& syst, int dom, const char* n, Tensor& so) {
  char* name;
  name = new char[LMAX];
  trim_spaces(name, n);
  int which = -1 ;
  int valence;
  char* name_ind = 0x0  ;
  Array<int>* type_ind = 0x0 ;
  bool found = false ;

  found = syst.iscst (name, which, valence, name_ind, type_ind) ;
  if(found) {
    Term_eq* mm = syst.give_cst(which, dom);
    mm->set_val_t(so);
  }
  return found;
}

template<class space_t>
vec_ary_t default_binary_vector_ary(space_t& space) {
  vec_ary_t coord_vectors{};
  Base_tensor basis (space, CARTESIAN_BASIS) ;
  for(auto i = 0; i < NUM_VECTORS; ++i) {
    switch(i) {
      case GLOBAL_ROT:
      case BCO1_ROT:
      case BCO2_ROT:
        coord_vectors[i] = Vector(space, CON, basis);
        break;
      default:
        coord_vectors[i] = Vector(space, COV, basis);
    }
  }
  return coord_vectors;
}

template<class space_t>
vec_ary_t default_co_vector_ary(space_t& space) {
  vec_ary_t coord_vectors{};
  Base_tensor basis (space, CARTESIAN_BASIS) ;
  coord_vectors[GLOBAL_ROT] = Vector(space,CON,basis);
  coord_vectors[BCO1_ROT]   = Vector(space,CON,basis);
  coord_vectors[EX]         = Vector(space,COV,basis);
  coord_vectors[EY]         = Vector(space,COV,basis);
  coord_vectors[EZ]         = Vector(space,COV,basis);
  coord_vectors[S_BCO1]     = Vector(space,COV,basis);
  coord_vectors[S_INF]      = Vector(space,COV,basis);
  return coord_vectors;
}
}
