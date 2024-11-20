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
 * [input] KMQ_vals: vector to be manipulated for gij
 * [input] order_: Interpolation order
 * [input] dr_: spacing of points to use in interpolation
 * [input] offset_: offset from excision to start interpolation
 * [input] r_: radius to evaluate
 * [input] r_bound_: initial guess of boundary radius
 * [input] theta_: angle from z to xy plane [0, pi]
 * [input] theta_: angle inside xy plane [0, 2pi]
 * [input] dom_: dom just outside excision surface
 * [input] extrapolate_metric: toggle whether gij is extrapolated inside
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

std::string throw_no_multithreaded_support_error(std::string not_implemented);
/** @}*/
}