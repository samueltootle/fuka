#pragma once
#include <vector>
#include "kadath.hpp"
#include "point.hpp"

//class Point ;
/**
 * @namespace export_utils
 * Utilities for exporting Kadath initial data to evolution kits such as ETK
 */
namespace export_utils {

std::vector<int> const R2TensorSymmetricIndices = {0, 1, 2, 4, 5, 8};

// clang-format off
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

// clang-format on

/**
 * generate the kth lagrange polynomial
 *
 * [input] interp_order: polynomial interpolation order
 * [input] x: current x
 * [input] xp: pointer to array of x points
 * [input] yp: pointer to array of y points
 */
double lagrange_gen_k(int interp_order,
                      double x,
                      const double* xp,
                      const double* yp);

/**
 * Create a Kadath point in Cartesian coordinates from the given input spherical coordinates
 *
 * [input] r: radius
 * [input] theta: angle theta
 * [input] phi: angle phi
 * [input] shift_x: offset of x
 */
Kadath::Point point_spherical(double r,
                              double theta,
                              double phi,
                              double shift_x);

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
template <class T, class space_t>
T get_excision_r(space_t const& space,
                 T r,
                 T const& theta_,
                 T const& phi_,
                 int const dom_,
                 T const xshift_) {
  // Generate cartesian point based on (r,t,p)
  auto p = point_spherical(r, theta_, phi_, xshift_);

  // bool to ascertain if the point is inside the excision region
  auto is_in_excision = [&](auto start, auto stop) -> bool {
    bool result = false;
    for (auto d = start; d < stop; ++d)
      result = result || space.get_domain(d)->is_in(p);
    return result;
  };

  size_t cnt{0};
  constexpr size_t max_iter = 1000;
  // Excision region consists of two domains, hence, dom_ - 2
  while (is_in_excision(dom_ - 2, dom_)) {
    r *= 1.01;
    p = point_spherical(r, theta_, phi_, xshift_);
    cnt++;

    // FIXME in some cases this is not very efficient - fails for max_iter = 100
    // usually at the pole
    if (cnt > max_iter) {
      cerr << p << "not found by radius increase. inc failed at " << r << ".\n";
      std::_Exit(EXIT_FAILURE);
    }
  };

  // the point should lay in the domain just outside the excision surface
  // if it's not, there must be a problem.
  if (!space.get_domain(dom_)->is_in(p)) {
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
template <class T = double, class space_t>
T interpolate_radial(space_t const& space,
                     Kadath::Scalar const& field_in,
                     int const order_,
                     T const dr_,
                     T const offset_,
                     T const ah_r,
                     T const r_,
                     T const theta_,
                     T const phi_,
                     const int dom_,
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
void add_tensor_refs(qarray_t& quants,
                     std::vector<int>&& ary_indicies,
                     Kadath::Tensor& field) {
  auto c = 0;
  for (auto i : ary_indicies) {
    auto tidx = export_utils::R2TensorSymmetricIndices[c];
    Kadath::Array<int> ind(field.indices(tidx));
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
template <class T,
          class quant_ary_t,
          class fields_ary_t,
          size_t N = NUM_VQUANTS>
void spherical_turduck(fields_ary_t& quants,
                       quant_ary_t& quant_vals,
                       int const order_,
                       T const dr_,
                       T const offset_,
                       T const r_bound_,
                       T const r_,
                       T const theta_,
                       T const phi_,
                       const int dom_,
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
    if constexpr (N == NUM_QUANTS) {
      if (k == H) {
        quant_vals[k] = 0;
        continue;
      } else if (k == UX || k == UY || k == UZ) {
        quant_vals[k] = 0;
        continue;
      }
    }
    for (int j = 0; j < order_; j++) {
      auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
      vals[j] = quants[k].get().val_point(p);
    }

    quant_vals[k] = lagrange_gen_k(order_, r_, r_points.data(), vals.data());
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
template <class T,
          class quant_ary_t,
          class fields_ary_t,
          size_t N = NUM_VQUANTS>
void spherical_turduck_fit_origin(fields_ary_t& quants,
                                  quant_ary_t& quant_vals,
                                  quant_ary_t& quant_vals_origin,
                                  int const order_,
                                  T const dr_,
                                  T const offset_,
                                  T const r_bound_,
                                  T const r_,
                                  T const theta_,
                                  T const phi_,
                                  const int dom_,
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
    if constexpr (N == NUM_QUANTS) {
      if (k == H) {
        quant_vals[k] = 0;
        continue;
      } else if (k == UX || k == UY || k == UZ) {
        quant_vals[k] = 0;
        continue;
      }
    }
    for (int j = 1; j < order_; j++) {
      auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
      vals[j] = quants[k].get().val_point(p);
    }

    quant_vals[k] = lagrange_gen_k(order_, r_, r_points.data(), vals.data());
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
template <class quant_ary_t, class fields_ary_t>
void spherical_turduck__gf(fields_ary_t& quants,
                           quant_ary_t& quant_vals,
                           int const order_,
                           std::vector<double> r_points,
                           double const r_,
                           double const theta_,
                           double const phi_,
                           int const gf,
                           double const xshift_) {

  std::vector<double> vals(order_);

  for (int j = 0; j < order_; j++) {
    auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
    vals[j] = quants[gf].get().val_point(p);
  }

  quant_vals[gf] = lagrange_gen_k(order_, r_, r_points.data(), vals.data());
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
template <class quant_ary_t, class fields_ary_t, class slice_t>
void spherical_turduck__all_gfs(fields_ary_t& quants,
                                quant_ary_t& quant_vals,
                                slice_t& slice,
                                double const x,
                                double const y,
                                double const z,
                                double const ah_r,
                                double extrap_r,
                                int const order_,
                                double const dr_,
                                double const offset_,
                                double const bh_ori) {

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

  for (auto& gf : slice) {
    // Where the filling takes places
    export_utils::spherical_turduck__gf(quants, quant_vals, order_, r_points,
                                        extrap_r, theta, phi, gf, bh_ori);
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
template <class quant_ary_t, class fields_ary_t>
void spherical_turduck_fit_origin__gf(fields_ary_t& quants,
                                      quant_ary_t& quant_vals,
                                      quant_ary_t& quant_vals_origin,
                                      int const order_,
                                      std::vector<double> r_points,
                                      double const r_,
                                      double const theta_,
                                      double const phi_,
                                      int const gf,
                                      double const xshift_) {

  std::vector<double> vals(order_);
  vals[0] = quant_vals_origin[gf];

  for (int j = 1; j < order_; j++) {
    auto p = point_spherical(r_points[j], theta_, phi_, xshift_);
    vals[j] = quants[gf].get().val_point(p);
  }

  quant_vals[gf] = lagrange_gen_k(order_, r_, r_points.data(), vals.data());
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
template <class quant_ary_t, class fields_ary_t, class slice_t>
void spherical_turduck_fit_origin__all_gfs(fields_ary_t& quants,
                                           quant_ary_t& quant_vals,
                                           quant_ary_t& quant_vals_origin,
                                           slice_t& slice,
                                           double const x,
                                           double const y,
                                           double const z,
                                           double const ah_r,
                                           double extrap_r,
                                           int const order_,
                                           double const dr_,
                                           double const offset_,
                                           double const bh_ori) {

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

  for (auto& gf : slice) {
    export_utils::spherical_turduck_fit_origin__gf(quants, quant_vals,
                                                   quant_vals_origin, order_,
                                                   r_points, extrap_r, theta,
                                                   phi, gf, bh_ori);
  }
}

template <class bh_exporter_t>
void partial_fill_excision(bh_exporter_t& bh_exporter,
                           const int dom_to_fill,
                           const int src_dom) {
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
  auto& conf = bh_exporter.get_conformal_factor();
  auto& shift = bh_exporter.get_shift();

  do {
    if (pos(0) < npts(0) - 1) {
      double x = xx(pos);
      double y = yy(pos);
      double z = zz(pos);
      bh_exporter.interpolate_pointwise__solution_gfs(x, y, z);
      auto qv = bh_exporter.get__quant_vals();
      conf->set_domain(dom_to_fill).set(pos) =
          qv[bh_exporter_t::XCTS_VARS::XCTS_PSI];
      lapse->set_domain(dom_to_fill).set(pos) =
          qv[bh_exporter_t::XCTS_VARS::XCTS_ALPHA];

      shift->set(1).set_domain(dom_to_fill).set(pos) =
          qv[bh_exporter_t::XCTS_VARS::XCTS_BETA1];
      shift->set(2).set_domain(dom_to_fill).set(pos) =
          qv[bh_exporter_t::XCTS_VARS::XCTS_BETA2];
      shift->set(3).set_domain(dom_to_fill).set(pos) =
          qv[bh_exporter_t::XCTS_VARS::XCTS_BETA3];
    } else {
      Index bcpos(pos);
      bcpos.set(0) = 0;
      conf->set_domain(dom_to_fill).set(pos) = conf->set_domain(src_dom)(bcpos);
      lapse->set_domain(dom_to_fill).set(pos) =
          lapse->set_domain(src_dom)(bcpos);
      for (int i = 1; i <= 3; ++i) {
        shift->set(i).set_domain(dom_to_fill).set(pos) =
            shift->set(i).set_domain(src_dom)(bcpos);
      }
    }
  } while (pos.inc());

  conf->std_base();
  lapse->std_base();
  shift->std_base();
}

std::string throw_no_multithreaded_support_error(std::string not_implemented);

struct basis_transform_spherical_tofrom_cart {
  static constexpr int dim = 3;
  using vector_t = std::array<double, dim>;
  using matrix_t = std::array<vector_t, dim>;

 private:
  matrix_t Jac_dSph_dCart;
  matrix_t Jac_dCart_dSph;

  double r_{0.};
  double theta_{0.};
  double phi_{0.};

  void set__Jacobians() {
    using std::cos;
    using std::sin;
    const double sint = sin(theta_);
    const double cost = cos(theta_);
    const double sinp = sin(phi_);
    const double cosp = cos(phi_);

    const double rsq = r_ * r_;
    const double rsint = r_ * sint;
    const double rsint_sq = rsint * rsint;

    const double div_r = 1. / r_;
    const double div_rsint = 1. / rsint;

    // dr^i/dx^i'
    Jac_dSph_dCart[0][0] = sint * cosp;
    Jac_dSph_dCart[0][1] = sint * sinp;
    Jac_dSph_dCart[0][2] = cost;

    // dtheta^i/dx^i'
    Jac_dSph_dCart[1][0] = cost * cosp * div_r;
    Jac_dSph_dCart[1][1] = cost * sinp * div_r;
    Jac_dSph_dCart[1][2] = -sint * div_r;

    // dphi^i/dx^i'
    Jac_dSph_dCart[2][0] = -sinp * div_rsint;
    Jac_dSph_dCart[2][1] = cosp * div_rsint;
    Jac_dSph_dCart[2][2] = 0.0;

    // dx^i/dx^i'
    Jac_dCart_dSph[0][0] = Jac_dSph_dCart[0][0];
    Jac_dCart_dSph[0][1] = Jac_dSph_dCart[1][0] * rsq;
    Jac_dCart_dSph[0][2] =
        Jac_dSph_dCart[2][0] * rsint * rsint;  //-sint * sinp * r_ ;

    // dy^i/dx^i'
    Jac_dCart_dSph[1][0] = Jac_dSph_dCart[0][1];
    Jac_dCart_dSph[1][1] = Jac_dSph_dCart[1][1] * rsq;
    Jac_dCart_dSph[1][2] = Jac_dSph_dCart[2][1] * rsint_sq;

    // dz^i/dx^i'
    Jac_dCart_dSph[2][0] = Jac_dSph_dCart[0][2];
    Jac_dCart_dSph[2][1] = Jac_dSph_dCart[1][2] * rsq;
    Jac_dCart_dSph[2][2] = 0.0;
  }

 public:
  void set__cart(const double x, const double y, const double z) {
    using std::sqrt;
    const double rsq = x * x + y * y + z * z;
    r_ = sqrt(rsq);
    if (std::fabs(r_) < 1e-12) {
      r_ = 1e-12;
      theta_ = 1e-12;
      phi_ = 1e-12;
    } else if (std::fabs(x) < 1e-12 && std::fabs(y) < 1e-12) {
      theta_ = 1e-12;
      phi_ = 1e-12;
    } else {
      theta_ = acos(z / r_);  //theta (angle from Z to xy plane)
      phi_ = atan2(y, x);     // phi (angle from x to y axis)
    }
    set__Jacobians();
  }

  void set__sph(const double r, const double theta, const double phi) {
    r_ = r;
    theta_ = theta;
    phi_ = phi;
    set__Jacobians();
  }

  basis_transform_spherical_tofrom_cart() = default;

  vector_t vectorU__sph_to_cart(const vector_t& V__sph) const {
    vector_t res_vectorU;
    for (int iCart = 0; iCart < dim; ++iCart) {
      res_vectorU[iCart] = 0.;
      for (int iSph = 0; iSph < dim; ++iSph) {
        res_vectorU[iCart] += Jac_dCart_dSph[iCart][iSph] * V__sph[iSph];
      }
    }
    return res_vectorU;
  }

  vector_t vectorD__sph_to_cart(const vector_t& V__sph) const {
    vector_t res_vectorD;
    for (int iCart = 0; iCart < dim; ++iCart) {
      res_vectorD[iCart] = 0.;
      for (int iSph = 0; iSph < dim; ++iSph) {
        res_vectorD[iCart] += Jac_dSph_dCart[iSph][iCart] * V__sph[iSph];
      }
    }
    return res_vectorD;
  }

  matrix_t matrixDD__sph_to_cart(const matrix_t& MDD__sph) const {
    matrix_t res_matrixDD;
    for (auto& col : res_matrixDD)
      col.fill(0);

    for (int iSph = 0; iSph < dim; ++iSph) {
      auto tmp_CartvectorD = vectorD__sph_to_cart(MDD__sph[iSph]);

      for (int iCart = 0; iCart < dim; ++iCart) {
        for (int jCart = 0; jCart < dim; ++jCart) {
          res_matrixDD[iCart][jCart] +=
              Jac_dSph_dCart[iSph][iCart] * tmp_CartvectorD[jCart];
        }
      }
    }
    return res_matrixDD;
  }

  vector_t vectorD__cart_to_sph(const vector_t& V__sph) const {
    vector_t res_vectorD;
    for (int iSph = 0; iSph < dim; ++iSph) {
      res_vectorD[iSph] = 0.;
      for (int iCart = 0; iCart < dim; ++iCart) {
        res_vectorD[iSph] += Jac_dCart_dSph[iCart][iSph] * V__sph[iCart];
      }
    }
    return res_vectorD;
  }

  matrix_t matrixDD__cart_to_sph(const matrix_t& MDD__cart) const {
    matrix_t res_matrixDD;
    for (auto& col : res_matrixDD)
      col.fill(0);

    for (int iCart = 0; iCart < dim; ++iCart) {
      auto tmp_SphvectorD = vectorD__cart_to_sph(MDD__cart[iCart]);

      for (int iSph = 0; iSph < dim; ++iSph) {
        for (int jSph = 0; jSph < dim; ++jSph) {
          res_matrixDD[iSph][jSph] +=
              Jac_dCart_dSph[iCart][iSph] * tmp_SphvectorD[jSph];
        }
      }
    }
    return res_matrixDD;
  }

  double get_r() const { return r_; }

  double get_theta() const { return theta_; }

  double get_phi() const { return phi_; }

  double get_r_sq() const { return r_ * r_; }
};

/** @}*/
}  // namespace export_utils