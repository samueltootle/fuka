#pragma once
#include "point.hpp"
#include <vector>
#include "kadath.hpp"

//class Point ;
/** 
 * @namespace export_utils
 * Utilities for exporting Kadath initial data to evolution kits such as ETK
 */
namespace export_utils {

std::vector<int> const R2TensorSymmetricIndices = {0, 1, 2, 4, 5, 8};
/// id quantities
enum {
  PSI,
  ALP,
  BETX,
  BETY,
  BETZ,
  AXX,
  AXY,
  AXZ,
  AYY,
  AYZ,
  AZZ,
  NUM_VQUANTS
};

/// chained enumeration to include matter quantities
enum id_matter_quants {
  H = NUM_VQUANTS,
  UX,
  UY,
  UZ,
  NUM_QUANTS
};

/// enumeration for quantities to export to evolution codes
enum sim_vac_quants {
  ALPHA,
  BETAX,
  BETAY,
  BETAZ,
  GXX,
  GXY,
  GXZ,
  GYY,
  GYZ,
  GZZ,
  KXX,
  KXY,
  KXZ,
  KYY,
  KYZ,
  KZZ,
  NUM_VOUT
};

/// chained enums together to include matter quantities
enum sim_matter_quants {
  RHO=NUM_VOUT,
  EPS,
  PRESS,
  VELX,
  VELY,
  VELZ,
  NUM_OUT
};

/**
 * generate the kth lagrange polynomial
 *
 * [input] interp_order: polynomial interpolation order
 * [input] x: current x
 * [input] xp: pointer to array of x points
 * [input] yp: pointer to array of y points
 */
double lagrange_gen_k(int interp_order, double x, const double *xp, const double *yp);

/**
 * Create a Kadath point in Cartesian coordinates from the given input spherical coordinates
 *
 * [input] r: radius
 * [input] theta: angle theta
 * [input] phi: angle phi
 * [input] shift_x: offset of x
 */
Kadath::Point point_spherical(double r, double theta, double phi, double shift_x);

/**
 * Need to refine the excision radius to ensure that it is outside
 * the excision surface otherwise we get NaNs if we
 * interpolate inside or worse - wrong data.
 *
 * [input] r_: radius guess
 * [input] theta_: fixed theta
 * [input] phi_: fixed angle phi
 * returns: corrected radius approximately equal to the excision radius
 */
template<class T, class space_t>
T get_excision_r (space_t const & space, 
  T r, T const & theta_, T const & phi_, int const dom_, T const xshift_) {
  // Generate cartesian point based on (r,t,p)
  auto p = point_spherical(r, theta_, phi_, xshift_);

  // bool to ascertain if the point is inside the excision region
  auto is_in_excision = [&] (auto start, auto stop) -> bool  {
    bool result = false;
    for(auto d = start; d < stop; ++d)
      result = result || space.get_domain(d)->is_in(p);
    return result;
  };  
  
  size_t cnt{0};
  constexpr size_t max_iter = 1000;
  // Excision region consists of two domains, hence, dom_ - 2
  while(is_in_excision(dom_ - 2, dom_)) {
    r *= 1.01;
    p = point_spherical(r, theta_, phi_, xshift_);
    cnt++;

    // FIXME in some cases this is not very efficient - fails for max_iter = 100
    // usually at the pole
    if(cnt > max_iter){
      cerr<< p << "not found by radius increase. inc failed at " << r << ".\n";
      std::_Exit(EXIT_FAILURE);
    }
  };
  
  // the point should lay in the domain just outside the excision surface
  // if it's not, there must be a problem.
  if(!space.get_domain(dom_)->is_in(p)){
    cerr << p << " not in dom " << dom_ << endl;
    std::_Exit(EXIT_FAILURE);
  }

  return r;
}

/**
 * interpolate_radial
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface for a given field before extrapolating inside the
 * excision surface
 * 
 * [input] space: numerical space
 * [input] field_in: field to interpolate on
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class T = double, class space_t>
T interpolate_radial(space_t const & space, Kadath::Scalar const & field_in,
  int const order_, T const dr_, T const offset_, T const ah_r, 
  T const r_, T const theta_, T const phi_, const int dom_, 
  T const xshift_) {
    
  // build vector of radial points
  std::vector<T> r_points(order_);
  for (int j = 0; j < order_; j++) {
    r_points[j] = (1. + offset_) * (1. + j * dr_) * ah_r;
  }

  // build vector of points to interpolate
  std::vector<T> vals(order_);
  for (int j = 0; j < order_; j++) {
    auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
    vals[j] = field_in.val_point(p);
  }

  return lagrange_gen_k(order_, r_, r_points.data(), vals.data());
}

template <class qarray_t>
void add_tensor_refs(qarray_t& quants, std::vector<int>&& ary_indicies, Kadath::Tensor& field) {
  auto c = 0;
  for(auto i : ary_indicies) {
    auto tidx = export_utils::R2TensorSymmetricIndices[c];
    Kadath::Array<int> ind (field.indices(tidx));
    quants[i] = std::cref(field(ind));
    c++;
  }
}

/**
 * spherical_turduck
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface before filling excision
 * 
 * [input] quant_vals: vector to be manipulated for alp, psi, bet, and kij
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class T, class quant_ary_t, class fields_ary_t, size_t N = NUM_VQUANTS>
void spherical_turduck(fields_ary_t& quants, quant_ary_t& quant_vals,
  int const order_, T const dr_, T const offset_, T const r_bound_, 
  T const r_, T const theta_, T const phi_, const int dom_, 
  T const xshift_) {
  
  //auto& space = quants[0].get().get_space();
  T const ah_r = r_bound_;
  // export_utils::get_excision_r(space,
  //   r_bound_, theta_, phi_, dom_, xshift_);
  
  std::vector<T> r_points(order_);
  for (int j = 0; j < order_; j++) {
    r_points[j] = (1. + offset_) * (1. + j * dr_) * ah_r;
  }

  for (size_t k = 0; k < N; ++k) {
    std::vector<T> vals(order_);

    // Avoid computations if fluid quantities are encountered
    // Necessary only for BHNS
    if constexpr(N == NUM_QUANTS) {
      if(k == H) {
        quant_vals[k] = 0;
        continue;
      }
      else if(k == UX || k == UY || k == UZ) {
        quant_vals[k] = 0;
        continue;
      }
    } 
    for (int j = 0; j < order_; j++) {
      auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
      vals[j] = quants[k].get().val_point(p);
    }
    
    quant_vals[k] =
      lagrange_gen_k(order_, r_, r_points.data(), vals.data());
  }
}

/**
 * spherical_turduck_fit_origin
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface before filling excision but matching to the provided
 * values (quant_vals_origin) for the grid functions at the origin.
 * 
 * [input] quant_vals: vector to be manipulated for alp, psi, bet, and kij
 * [input] quant_vals_origin: values of grid functions at the origin to be fit to
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class T, class quant_ary_t, class fields_ary_t, size_t N = NUM_VQUANTS>
void spherical_turduck_fit_origin(fields_ary_t& quants, quant_ary_t& quant_vals, quant_ary_t& quant_vals_origin,
  int const order_, T const dr_, T const offset_, T const r_bound_, 
  T const r_, T const theta_, T const phi_, const int dom_, 
  T const xshift_) {
  
  //auto& space = quants[0].get().get_space();
  T const ah_r = r_bound_;
  
  std::vector<T> r_points(order_);
  r_points[0] = 0.;
  for (int j = 1; j < order_; j++) {
    r_points[j] = (1. + offset_) * (1. + j * dr_) * ah_r;
  }

  for (size_t k = 0; k < N; ++k) {
    std::vector<T> vals(order_);
    vals[0] = quant_vals_origin[k];

    // Avoid computations if fluid quantities are encountered
    // Necessary only for BHNS
    if constexpr(N == NUM_QUANTS) {
      if(k == H) {
        quant_vals[k] = 0;
        continue;
      }
      else if(k == UX || k == UY || k == UZ) {
        quant_vals[k] = 0;
        continue;
      }
    } 
    for (int j = 1; j < order_; j++) {
      auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
      vals[j] = quants[k].get().val_point(p);
    }
    
    quant_vals[k] =
      lagrange_gen_k(order_, r_, r_points.data(), vals.data());
  }
}

/**
 * spherical_turduck
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface before filling excision
 * 
 * [input] quant_vals: vector to be manipulated for alp, psi, bet, and kij
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class quant_ary_t, class fields_ary_t>
void spherical_turduck__gf(fields_ary_t& quants, quant_ary_t& quant_vals,
  int const order_, std::vector<double> r_points,
  double const r_, double const theta_, double const phi_, int const gf, double const xshift_) {

  std::vector<double> vals(order_);

  for (int j = 0; j < order_; j++) {
    auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
    vals[j] = quants[gf].get().val_point(p);
  }
  
  quant_vals[gf] =
    lagrange_gen_k(order_, r_, r_points.data(), vals.data());
}

/**
 * spherical_turduck
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface before filling excision
 * 
 * [input] quant_vals: vector to be manipulated for alp, psi, bet, and kij
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class quant_ary_t, class fields_ary_t, class slice_t>
void spherical_turduck__all_gfs(fields_ary_t& quants, quant_ary_t& quant_vals, slice_t& slice,
  double const x, double const y, double const z, double const ah_r, double extrap_r,
  int const order_, double const dr_, double const offset_, double const bh_ori) {

  // Avoid division by "0"
  extrap_r = (extrap_r <= 1e-12) ? 1e-12 : extrap_r;
  double x_shifted = (x == 0.) ? 1e-14 : x;
  double theta = std::acos(z / extrap_r);
  
  // atan2 is needed here
  double phi = std::atan2(y, (x_shifted - bh_ori)); 

  std::vector<double> r_points(order_);
  for (int j = 0; j < order_; j++) {
    r_points[j] = (1. + offset_) * (1. + j * dr_) * ah_r;
  }

  for(auto & gf : slice) {
    // Where the filling takes places
    export_utils::spherical_turduck__gf(
      quants, quant_vals, order_, r_points, extrap_r, theta, phi, gf, bh_ori
    );
  }
}


/**
 * spherical_turduck_fit_origin
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface before filling excision but matching to the provided
 * values (quant_vals_origin) for the grid functions at the origin.
 * 
 * [input] quant_vals: vector to be manipulated for alp, psi, bet, and kij
 * [input] quant_vals_origin: values of grid functions at the origin to be fit to
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class quant_ary_t, class fields_ary_t>
void spherical_turduck_fit_origin__gf(fields_ary_t& quants, quant_ary_t& quant_vals, quant_ary_t& quant_vals_origin,
  int const order_, std::vector<double> r_points,
  double const r_, double const theta_, double const phi_, int const gf, double const xshift_) {

  std::vector<double> vals(order_);
  vals[0] = quant_vals_origin[gf];

  for (int j = 1; j < order_; j++) {
    auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
    vals[j] = quants[gf].get().val_point(p);
  }
  
  quant_vals[gf] =
    lagrange_gen_k(order_, r_, r_points.data(), vals.data());
}

/**
 * spherical_turduck
 * 
 * We interpolate radially in 3D - f(r, theta, phi) - for an arbitrary
 * excision surface before filling excision
 * 
 * [input] quant_vals: vector to be manipulated for alp, psi, bet, and kij
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] phi_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] xshift_: coordinate x shift (for binaries)
 */
template<class quant_ary_t, class fields_ary_t, class slice_t>
void spherical_turduck_fit_origin__all_gfs(fields_ary_t& quants, quant_ary_t& quant_vals, quant_ary_t& quant_vals_origin,
  slice_t& slice, double const x, double const y, double const z, double const ah_r, double extrap_r,
  int const order_, double const dr_, double const offset_, double const bh_ori) {

  // Avoid division by "0"
  extrap_r = (extrap_r <= 1e-12) ? 1e-12 : extrap_r;
  double x_shifted = (x == 0.) ? 1e-14 : x;
  double theta = std::acos(z / extrap_r);
  
  // atan2 is needed here
  double phi = std::atan2(y, (x_shifted - bh_ori)); 

  std::vector<double> r_points(order_);
  r_points[0] = 0.;
  for (int j = 1; j < order_; j++) {
    r_points[j] = (1. + offset_) * (1. + j * dr_) * ah_r;
  }

  for(auto & gf : slice) {
    export_utils::spherical_turduck_fit_origin__gf(
      quants, quant_vals, quant_vals_origin, order_, r_points, 
      extrap_r, theta, phi, gf, bh_ori
    );
  }
}


template<class bh_exporter_t>
void partial_fill_excision(bh_exporter_t& bh_exporter, const int dom_to_fill, const int src_dom) {
  using namespace Kadath;
  
  // We set this to false to ensure that bh_exporter.interpolate_pointwise will not fit to origin
  bh_exporter.set__export_ready(false);
  
  auto dom = bh_exporter.get_space()->get_domain(dom_to_fill);
  auto npts = dom->get_nbr_points();
  
  Index pos(npts);
  Val_domain xx = dom->get_cart(1);
  Val_domain yy = dom->get_cart(2);
  Val_domain zz = dom->get_cart(3);

  auto& lapse = bh_exporter.get_lapse();
  auto& conf  = bh_exporter.get_conformal_factor();
  auto& shift = bh_exporter.get_shift();

  do {
    if(pos(0) < npts(0)-1){
      double x = xx(pos);
      double y = yy(pos);
      double z = zz(pos);
      bh_exporter.interpolate_pointwise__solution_gfs(x, y, z);
      auto qv = bh_exporter.get__quant_vals();
      conf->set_domain(dom_to_fill).set(pos) = qv[bh_exporter_t::XCTS_VARS::XCTS_PSI];
      lapse->set_domain(dom_to_fill).set(pos) = qv[bh_exporter_t::XCTS_VARS::XCTS_ALPHA];
      
      shift->set(1).set_domain(dom_to_fill).set(pos) = qv[bh_exporter_t::XCTS_VARS::XCTS_BETA1];
      shift->set(2).set_domain(dom_to_fill).set(pos) = qv[bh_exporter_t::XCTS_VARS::XCTS_BETA2];
      shift->set(3).set_domain(dom_to_fill).set(pos) = qv[bh_exporter_t::XCTS_VARS::XCTS_BETA3];
    } else {
      Index bcpos(pos);
      bcpos.set(0) = 0;
      conf->set_domain(dom_to_fill).set(pos) = conf->set_domain(src_dom)(bcpos);
      lapse->set_domain(dom_to_fill).set(pos) = lapse->set_domain(src_dom)(bcpos);
      for(int i = 1; i <=3; ++i) {
        shift->set(i).set_domain(dom_to_fill).set(pos) = shift->set(i).set_domain(src_dom)(bcpos);
      }
    }
  }while(pos.inc());

  conf->std_base();
  lapse->std_base();
  shift->std_base();
}
/** @}*/
}
