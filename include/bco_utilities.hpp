/*
 * This file is part of the KADATH library.
 * Author: Samuel Tootle
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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "kadath.hpp"
#include "kadath_adapted.hpp"

/**
 * \addtogroup FUKA_BCO_Utils
 * @ingroup FUKA
 * @{*/

namespace Kadath {
namespace bco_utils {
// using namespace Kadath;
using namespace ::Kadath::FUKA_Config;

// approximate value of psi on the horizon
// to make coordinate estimates
constexpr double psi = 1.55;
constexpr double psisq = psi * psi;
constexpr double invpsisq = 1. / psisq;

static constexpr double shell_factor =
    std::log(1.855) / (100.0 * std::log(2.0));

#define EQUI -11
#define INNER_EQUI -12

// this is in principle a joke and came from the need to have
// a scaling factor (R) s.t. 1.5 < R < 2.
// Since this works well, it embarrassingly persists...
const double gold_ratio = (1 + std::sqrt(5.)) / 2.;

/**
 * bco_utils::KadathPNOrbitalParams \n
 * 3.5 PN orbital parameters for irrotational compact object
 * This is a very rough, but useful guess for any CO
 *
 * @param[input]  bconfig:  binary Configurator file
 * @param[input]  Madm1: ADM mass (FIXED_MADM or MCH)
 * @param[input]  Madm2: ADM mass (FIXED_MADM or MCH)
 */
void KadathPNOrbitalParams(kadath_config_boost<BIN_INFO>& bconfig,
                           double Madm1,
                           double Madm2);

/**
 * bco_utils::save_to_file\n
 * Standard print config, space, and fields to their respective files
 * Config is saved to base_fname.info, space and fields in base_fname.dat
 * Output directory is stored in the config file. See config_utils_boost.hpp
 *
 * @tparam config_t Configurator type
 * @tparam space_t Space type
 * @tparam fields_t parameter pack of Kadath field types to save to file
 * @param[input]  base_fname: base filename (no extension)
 * @param[input]  bconfig:  Configurator file
 * @param[input]  space:  Numerical space
 * @param[input]  fields: Parameter pack of fields to save to file
 */
template <typename space_t, typename config_t, typename... fields_t>
void save_to_file(std::stringstream& base_fname,
                  space_t& space,
                  config_t& bconfig,
                  fields_t&... fields) {
  bconfig.set_filename(base_fname.str());
  bconfig.write_config();
  std::string kadath_filename = bconfig.space_filename();
  FILE* ff = fopen(kadath_filename.c_str(), "w");
  space.save(ff);

  (fields.save(ff), ...);
  fclose(ff);
}

/**
 * bco_utils::save_to_file
 *
 * Specialization that excludes the base_filename stringstream.\n
 * Standard print config, space, and fields to their respective files.\n
 *  - Config is saved to its set preset filename\n
 *  - Fields are save to the paired filename. See config_utils_boost.hpp.\n
 *
 * @tparam config_t Configurator type
 * @tparam space_t Space type
 * @tparam fields_t parameter pack of Kadath field types to save to file
 * @param[input]  space:  Numerical space
 * @param[input]  bconfig:  Configurator file
 * @param[input]  fields: Parameter pack of fields to save to file
 */
template <typename space_t, typename config_t, typename... fields_t>
void save_to_file(space_t& space, config_t& bconfig, fields_t&... fields) {
  bconfig.write_config();

  std::string kadath_filename = bconfig.space_filename();
  FILE* ff = fopen(kadath_filename.c_str(), "w");
  space.save(ff);

  (fields.save(ff), ...);
  fclose(ff);
}

/**
 * bco_utils::get_radius
 *
 * When all we want is a simply radius value (inner or outer radius)
 *
 * Note: For variable domains, this is just R at a single point on the surface
 *
 * @tparam dom_t Domain type
 * @param[input]  dom: pointer to a domain with type dom_t
 * @param[input]  bound: need to know if you want inner or outer radius
 * @return a radius of the domain.   */
template <typename dom_t>
double get_radius(const dom_t* dom, const int bound) {
  auto const npts = dom->get_nbr_points();
  Index pos(npts);
  auto dim = pos.get_ndim();
  switch (bound) {
    case INNER_BC:
      break;
    case OUTER_BC:
      pos.set(0) = npts(0) - 1;
      pos.set(1) = npts(1) - 1;
      if (dim == 3)
        pos.set(2) = npts(2) - 1;
      break;
    case EQUI:
      pos.set(0) = npts(0) - 1;
      pos.set(1) = npts(1) - 1;
      break;
    case INNER_EQUI:
      pos.set(1) = npts(1) - 1;
      break;
    default:
      std::cout << "Unknown bound sent to get_radius: " << bound << std::endl;
      break;
  }
  return dom->get_radius()(pos);
}

/**
 * bco_utils::get_min_max
 *
 * For variable domains, find the min and max of a Scalar field on a boundary
 *
 * @tparam space_t type of numerical space
 * @param[input]  field: Scalar field
 * @param[input]  dom: domain to find min and max of
 * @return {fmin, fmax} std::array containing fmin and fmax
 */
inline std::array<double, 2> get_field_min_max(const Scalar& field,
                                               const int dom,
                                               const int bound = OUTER_BC) {
  const int npts_r = field.get_domain(dom)->get_nbr_points()(0);
  const int dim = field.get_domain(dom)->get_ndim();

  // position index to loop over
  Index pos(field.get_domain(dom)->get_nbr_points());
  // position index to update with a fixed radial boundary
  Index bpos(field.get_domain(dom)->get_nbr_points());

  int r_bound = 0;
  switch (bound) {
    case INNER_BC:
      break;
    case OUTER_BC:
      r_bound = npts_r - 1;
      break;
    default:
      std::cout << "Unknown bound sent to get_field_min_max: " << bound
                << std::endl;
      break;
  }

  double fmax = field(dom)(pos);
  double fmin = field(dom)(pos);

  // FIXME this currently loops over theta and phi multiple times, but the syntax is straightforward.
  do {
    bpos.set(0) = r_bound;
    bpos.set(1) = pos(1);
    if (dim == 3)
      bpos.set(2) = pos(2);
    double f = field(dom)(bpos);

    if (f > fmax)
      fmax = f;
    if (f < fmin)
      fmin = f;
  } while (pos.inc());
  return std::array<double, 2>{fmin, fmax};
}

/**
 * bco_utils::get_rmin_rmax
 *
 * For variable domains, find the min and max radii
 *
 * @tparam space_t type of numerical space
 * @param[input]  space: Numerical space
 * @param[input]  dom: domain to find min and max radii of
 * @return {rmin, rmax} std::array containing rmin and rmax
 */
template <typename space_t>
std::array<double, 2> get_rmin_rmax(const space_t& space, const int dom) {
  const int npts_r = space.get_domain(dom)->get_nbr_points()(0);
  Index pos(space.get_domain(dom)->get_nbr_points());
  pos.set_start();
  pos.set(0) = npts_r - 1;

  double rmax = space.get_domain(dom)->get_radius()(pos);
  double rmin = space.get_domain(dom)->get_radius()(pos);
  do {
    pos.set(0) = npts_r - 1;
    double r = space.get_domain(dom)->get_radius()(pos);

    if (r > rmax)
      rmax = r;
    if (r < rmin)
      rmin = r;
  } while (pos.inc());
  return std::array<double, 2>{rmin, rmax};
}

/**
 * bco_utils::update_adapted_field
 *
 * extrapolate variable fields to domains that are/should be zero based on the coefficients on the boundary, e.g:
 * 	- extrapolate space-time fields into the horizon of the BH from outside the horizon
 * 	- extrapolate matter fields from inside the star to outside of the star
 *
 * Note 1:\n
 * This is important before importing fields from one numerical space to another as interpolating the fields to the new\n
 * space may have problems if the adjacent domain is zero.  After import, remember to set relevant domains back to zero\n
 *
 * Note 2:\n
 * could be generalized to Tensors.  Currently only done for Scalars
 *
 * @tparam dom_t type of domain.
 * @param[input]  res: reference to modifiable Scalar field that contains the current solution to the field.  This will be changed!
 * @param[input]  from_dom: index of domain to extrapolate FROM
 * @param[input]  to_dom: index of domain to extrapolate TO
 * @param[input]  dom: pointer to domain to extrapolate TO
 * @param[input]  bound: which boundary we extrapolating from (OUTER_BC, INNER_BC)
 * @param[output] res: Modified scalar field with values interpolated from from_dom to to_dom
 */
template <typename dom_t>
void update_adapted_field(Scalar& res,
                          const int from_dom,
                          const int to_dom,
                          const dom_t* dom,
                          const int bound) {
  Tensor* ptr = &res;
  Array<int> doms(2);
  doms.set(0) = from_dom;
  doms.set(1) = to_dom;
  Scalar import = dom->import(to_dom, bound, 1, doms, &ptr);
  ptr->set().set_domain(to_dom) = import(to_dom);
}

/**
 * bco_utils::interp_adapted_mapping\n
 * Update the mapping of an adapted domain based on an old numerical space before importing fields
 *
 * @tparam adapted_t type of adapted domain
 * @param[input]  new_shell: const pointer to adapted domain of the new shell that we're creating a mapping for
 * @param[input]  old_outer_adapted_dom: index of old outer adapted dom that we are creating a mapping from
 * @param[input]  old_radius_field: Scalar field describing the radii of all the fields in the numerical space to interpolate on
 */
template <typename adapted_t>
void interp_adapted_mapping(const adapted_t* new_shell,
                            const int old_outer_adapted_dom,
                            const Scalar& old_radius_field) {
  Val_domain new_mapping = new_shell->get_radius();
  int const ndim = new_shell->get_ndim();

  //Need to normalize by the constant radius at the INNER_BC of the outer_adapted shell in the old space
  auto old_shell =
      old_radius_field.get_space().get_domain(old_outer_adapted_dom);
  int const old_ndim = old_shell->get_ndim();
  double rinner = get_radius(old_shell, INNER_BC);

  Index new_pos(new_shell->get_nbr_points());
  double xc_new = new_shell->get_center()(1);
  double xc_old = old_shell->get_center()(1);
  do {

    auto interp_val = [&](auto const x, auto const y, auto const z) {
      Kadath::Point absol(old_ndim);
      switch (old_ndim) {
        case 2: {
          auto rsq_xy = x * x + y * y;
          auto r_xy = std::sqrt(rsq_xy);
          absol.set(1) = r_xy + xc_old;
          absol.set(2) = z;
          break;
        }
        case 3:
          absol.set(1) = x + xc_old;
          absol.set(2) = y;
          absol.set(3) = z;
          break;
      }
      return old_radius_field.val_point(absol);
    };
    double x, y, z;

    if (ndim == 3) {
      x = new_shell->get_cart(1)(new_pos) - xc_new;
      y = new_shell->get_cart(2)(new_pos);
      z = new_shell->get_cart(3)(new_pos);
    } else {
      x = new_shell->get_cart(1)(new_pos) - xc_new;
      y = 0;  // phi symmetry
      z = new_shell->get_cart(2)(new_pos);
    }
    double rsq = x * x + y * y + z * z;
    double r = std::sqrt(rsq);

    x /= r / rinner;
    y /= r / rinner;
    z /= r / rinner;

    new_mapping.set(new_pos) = interp_val(x, y, z);

  } while (new_pos.inc());

  new_mapping.std_base();
  new_shell->set_mapping(new_mapping);
}

/**
 * bco_utils::get_center
 * Quickly get the X coordinate of the  cartesian center of a domain\n
 * while using the correct Index and letting it be deleted when we're done.
 *
 * @tparam space_t type of numerical space
 * @param[input]  space: numerical space
 * @param[input]  dom: index of domain of interest
 * @return coordinate center on the X axis since we're always at (X,0,0)
 */
template <typename space_t>
double get_center(const space_t& space, const int dom) {
  Index pos(space.get_domain(dom)->get_nbr_points());
  return space.get_domain(dom)->get_cart(1)(pos);
}

/**
 * bco_utils::get_boundary_val
 * Quickly get the value at a specified boundary of a scalar field in a given domain\n
 * and letting the setup be deleted when we're done.
 *
 * @tparam space_t type of numerical space
 * @param[input]  space: numerical space
 * @param[input]  dom: index of domain of interest
 * @param[input]  field: scalar field of interest
 * @param[input]  bound: bondary of interest
 * @return value at the (0,0,0) collocation point
 */
inline double get_boundary_val(const int dom,
                               const Scalar& field,
                               const int bound = INNER_BC) {
  auto this_domain = field(dom).get_domain();
  auto npts = this_domain->get_nbr_points();
  Index pos(npts);
  auto dim = pos.get_ndim();

  switch (bound) {
    case INNER_BC:
      break;
    case OUTER_BC:
      pos.set(0) = npts(0) - 1;
      pos.set(1) = npts(1) - 1;
      if (dim == 3)
        pos.set(2) = npts(2) - 1;
      break;
    case EQUI:
      pos.set(0) = npts(0) - 1;
      pos.set(1) = npts(1) - 1;
      break;
    case INNER_EQUI:
      pos.set(1) = npts(1) - 1;
      break;
    default:
      std::cout << "Unknown bound sent to get_boundary_val: " << bound
                << std::endl;
      break;
  }
  return field(dom)(pos);
}

/**
 * bco_utils::set_radius
 * Quickly get the basic outer radius of a given domain to update the corresponding
 * radius in the config file based on the provides indexes
 *
 * @tparam space_t type of numerical space
 * @param[input]  dom: index of domain of interest
 * @param[input]  space: numerical space
 * @param[input]  bconfig: configuration file
 * @param[input]  idxs: 1 or more indexes related to the config file entry to update
 */
template <typename config_t, typename space_t, typename... idx_t>
void set_radius(const int& dom,
                const space_t& space,
                config_t& bconfig,
                const idx_t... idxs) {
  bconfig.set(idxs...) = get_radius(space.get_domain(dom), OUTER_BC);
}

/**
 * bco_utils::set_NS_bounds
 *
 * Set boundaries for a NS companion based on the config file parameters
 *
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds to update
 * @param[input]  bconfig: configuration file
 * @param[input]  bco index: of the NS in the configuration file
 */
template <typename ary_t, typename config_t, typename... idx_t>
void set_NS_bounds(ary_t& bounds, config_t& bconfig, idx_t... bco) {
  //sizes
  const int size = bounds.size();
  const int ninshells = (!std::isnan(bconfig.set(NINSHELLS, bco...)))
                            ? bconfig(NINSHELLS, bco...)
                            : 0;
  const int nshells = (!std::isnan(bconfig.set(NSHELLS, bco...)))
                          ? bconfig(NSHELLS, bco...)
                          : 0;

  //indexes
  const int rin = 0;
  const int rout = size - 1;
  const int r = rin + ninshells + 1;

  bounds[rin] = bconfig(RIN, bco...);
  bounds[r] = bconfig(RMID, bco...);
  bounds[rout] = bconfig(ROUT, bco...);

  double lower = bounds[rin] + (bounds[r] - bounds[rin]) * 0.8;
  double delta_r = (bounds[r] - lower) / (ninshells + 1.);

  for (int shell = rin + 1; shell <= ninshells; ++shell) {
    bounds[shell] = delta_r * shell + lower;
  }

  delta_r = (bounds[rout] - bounds[r]) / (nshells + 1.);
  for (int i = 1, shell = r + 1; i <= nshells; ++shell, ++i) {
    bounds[shell] = bounds[r] + delta_r * i;
  }
}

/**
 * bco_utils::set_isolated_BH_bounds
 *
 * Set boundaries for an isolated  BH based on the config file parameters
 *
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds to update
 * @param[input]  bconfig: configuration file
 */
template <typename ary_t, typename config_t>
void set_isolated_BH_bounds(ary_t& bounds, config_t& bconfig) {

  // Configurator index of companion compact object
  const int size = bounds.size();

  //indexes
  const int rin = 0;
  const int r = rin + 1;
  int rout = size - 1;

  bounds[rin] = bconfig(RIN);
  bounds[r] = bconfig(RMID);

  double delrBH = bconfig(ROUT) - bconfig(RMID);

  const int shells = bconfig(NSHELLS);
  double delr = delrBH / (bconfig(NSHELLS) + 1);

  for (int i = 0, b = r + 1; i < shells; ++b, ++i)
    bounds[b] = bounds[r] + delr * (i + 1);
  //    bounds[b] = bounds[r] + bconfig(ROUT) / M_PI * atan((i + 1)/bconfig(NSHELLS));

  bounds[rout] = bconfig(ROUT);
}

/**
 * bco_utils::set_BH_bounds -- deprecated
 *
 * Set boundaries for a BH companion based on the config file parameters
 * and it's companion
 *
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds to update
 * @param[input]  bconfig: configuration file
 * @param[input]  bco: index of the BH in the configuration file
 */
template <typename ary_t, typename config_t>
void set_BH_bounds(ary_t& bounds,
                   config_t& bconfig,
                   const int bco,
                   const bool adapt_shells = false) {

  // Configurator index of companion compact object
  const int bco2 = (bco == BCO1) ? BCO2 : BCO1;
  const int size = bounds.size();

  //indexes
  const int rin = 0;
  const int r = rin + 1;
  int rout = size - 1;

  bounds[rin] = bconfig(RIN, bco);
  bounds[r] = bconfig(RMID, bco);
  bounds[rout] = bconfig(ROUT, bco);

  const int shells = bconfig(NSHELLS, bco);
  double sum = 0.;
  ary_t new_ary;
  new_ary.push_back(bounds[rin]);
  new_ary.push_back(bounds[r]);

  // fixme - need to recall how this limit was obtained
  const double limit = bconfig(ROUT, bco) * invpsisq;

  // needed to determine proper shall spacing to not have problems in kadath import?
  const double scale_fact = bco_utils::gold_ratio;

  // y_intercept of new_bound == the location of the first shell
  const double y_intercept = 2. * scale_fact * bounds[r];

  // the intercept_fac ensures the y_intercept is just that
  const double intercept_fac = -scale_fact * std::log(y_intercept);

  // relation that provides shell locations
  auto new_bound = [&](int n) {
    return std::exp((n - intercept_fac) / scale_fact);
  };

  // define new shells bounds
  int N = 0;
  double next_bound = new_bound(N);
  while (next_bound < limit && next_bound < bconfig(ROUT, bco)) {
    new_ary.push_back(next_bound);
    N++;
    next_bound = new_bound(N);
  }
  // in case NSHELLS are put in by hand - we add additional ones based on a
  // naive even distribution
  if (new_ary.size() - 2 < bconfig(NSHELLS, bco)) {
    // determine how many shells have already been defined
    auto current_shells = new_ary.size() - 2;
    // how many shells still need to defined
    auto remaining_shells = shells - current_shells;

    // radius of the last defined shell
    auto last_shell_r = new_ary.back();
    // radial distance between last shell and ROUT
    double delrBH = bconfig(ROUT, bco) - last_shell_r;
    // determine the equal spacing between reminaing shells
    double delr = delrBH / (remaining_shells + 1);

    for (int i = 0; i < remaining_shells; ++i)
      new_ary.push_back(last_shell_r + delr * (i + 1));
  }
  // define the bound corresponding to ROUT
  new_ary.push_back(bounds[rout]);
  // update the number of shells in the event shells have been automatically
  // added
  bconfig(NSHELLS, bco) = new_ary.size() - 3;
  bounds = std::move(new_ary);
}

/**
 * bco_utils::gen_shell_bound_radius
 *
 * Recursive routine to Generate the next shell boundary radius
 * using the input field.
 * By default, field is intended to be \partial_r^2(Psi).
 * This has not been tested with other input fields.
 *
 * @tparam T return type (double most likely)
 * @param[input]  field: e.g. \partial_r^2(Psi)
 * @param[input]  r0: fixed inner radius
 * @param[input]  r1: outer radius guess
 * @return r1: new shell radius
 */
template <class T>
T gen_shell_bound_radius(Scalar& field,
                         T r0,
                         T r1,
                         T xc = 0.,
                         T threshold = 0.9,
                         T fac = 1.) {
  auto reldiff = [](auto ref, auto cmp) {
    return 1. - cmp / ref;
  };

  // Setup points relative to coordinate center
  Kadath::Point pt_r0(3);
  pt_r0.set(1) = xc + r0;

  Kadath::Point pt_r1(3);
  pt_r1.set(1) = xc + r1;
  // end point setup

  // relative difference of ddrP at both radii
  auto ddrP_rel_diff =
      std::abs(reldiff(field.val_point(pt_r0), field.val_point(pt_r1)));

  // relative difference to the desired threshold
  auto relth = std::abs(reldiff(threshold, ddrP_rel_diff));

#ifdef DEBUG
  std::cout << "r0: " << r0 << "\t ddrP0: " << field.val_point(pt_r0)
            << "\t P: " << pt_r0 << '\n'
            << "r1: " << r1 << "\t ddrP1: " << field.val_point(pt_r1)
            << "\t ddrP_rel_diff: " << ddrP_rel_diff << "\t P: " << pt_r1
            << '\n'
            << "\t relth: " << relth << "\t r1 / r0: " << r1 / r0 << '\n';
  std::cout << pt_r0 << "," << pt_r1 << '\n';
#endif

  auto dr = r1 - r0;

  if (ddrP_rel_diff > threshold) {
    auto new_fac = 0.5 * fac;
    r1 = r0 + dr * new_fac;
    r1 = gen_shell_bound_radius(field, r0, r1, xc, threshold, new_fac);
  } else {
    // If the radii are far enough apart we found a good radius
    if (r1 / r0 > bco_utils::gold_ratio)
      return r1;
    else {
      return r0 * bco_utils::gold_ratio;
    }
  }
  return r1;
}

/**
 * bco_utils::set_arb_bounds
 *
 * Generate bounds for the local shells around compact objects using
 * \partial_r^2(Psi) for binary initial data
 *
 * @tparam config_t Configurator type
 * @param[input]  bconfig: Configurator
 * @param[input]  bco: index of compact object data (BCO1/BCO2)
 * @param[input]  ddrconf_sol: Estimate or true solution of \partial_r^2(Psi)
 * @param[input]  adapted_dom_sol: domain index of the inner adapted boundary
 * @param[input]  threshold: the threshold to meet when finding a shell radius
 * @return bounds: vector of bounds
 */
template <typename config_t>
std::vector<double> set_arb_bounds(config_t& bconfig,
                                   const int bco,
                                   Scalar& ddrconf_sol,
                                   const int adapted_dom_sol,
                                   const double threshold = 0.95) {

  // Generate a new vector for bounds
  // FIXME: this ignores shells inside of NS'
  std::vector<double> bounds;
  bounds.push_back(bconfig(RIN, bco));
  bounds.push_back(bconfig(RMID, bco));

  auto adapted_dom(ddrconf_sol.get_space().get_domain(adapted_dom_sol));

  // Coordinate centered point
  Kadath::Point pt(adapted_dom->get_center());
  double xc{pt(1)};

  // Initial radii starting at the inner adapted radius on the
  // equitorial plane since this should be the largest radius
  auto r0{get_radius(adapted_dom, INNER_EQUI)};
  // We only add shells out to ROUT
  auto r1{bconfig(ROUT, bco)};
  auto Rout{r1};

  // some upper bound that should never be hit!
  auto max_shells = std::ceil((Rout - r0) / (3. * r0));
  for (auto i = 0; i < max_shells; ++i) {

    pt.set(1) = xc + r0;

    // stop looking for shells once \partial_r^2 \Psi
    // is roughly flat
    // FIXME:? this will fail at local extrema!
    if (ddrconf_sol.val_point(pt) < 0.01)
      break;

    // get new shell boundary outer radius
    r1 = gen_shell_bound_radius(ddrconf_sol, r0, r1, xc, threshold);

    // check for minimum spacing between last shell boundary
    // and Rout
    if (Rout / r1 > bco_utils::gold_ratio) {
#ifdef DEBUG
      cout << "Rout/r1: " << Rout / r1 << ", Success.\n";
#endif
      bounds.push_back(r1);
    } else {
#ifdef DEBUG
      cout << "Rout/r1: " << Rout / r1 << ", failed.\n";
#endif
      break;
    }

    // update variables for finding next shell
    r0 = r1;
    r1 = Rout;
    if (i == max_shells - 1) {
      std::cerr << "Max iterations hit in setting CO bounds. Something went "
                   "very wrong!\n";
      std::cerr << "r0: " << r0 << '\n'
                << "r1: " << r1 << '\n'
                << "max_shells: " << (Rout - r0) / (3. * r0) << '\n';
      std::_Exit(EXIT_FAILURE);
    }
  }
  // Add ROUT as the last bound
  bounds.push_back(Rout);
  bconfig(NSHELLS, bco) = bounds.size() - 3;
  return bounds;
}

/**
 * bco_utils::set_arb_boundsv3
 *
 * Generate bounds for the local shells around compact objects using
 * \partial_r^2(Psi) for binary initial data
 *
 * @tparam config_t Configurator type
 * @param[input]  bconfig: Configurator
 * @param[input]  ddrconf_sol: Estimate or true solution of \partial_r^2(Psi)
 * @param[input]  adapted_dom_sol: domain index of the inner adapted boundary
 * @param[input]  threshold: the threshold to meet when finding a shell radius
 * @return bounds: vector of bounds
 */
template <typename config_t, class... idx_t>
std::vector<double> set_arb_boundsv3(config_t& bconfig,
                                     Scalar& field,
                                     const int adapted_dom_sol,
                                     idx_t... BCOidx) {
  auto reldiff = [](auto ref, auto cmp) {
    return 1. - cmp / ref;
  };

  // Generate a new vector for bounds
  // FIXME: this ignores shells inside of NS'
  std::vector<double> bounds;
  bounds.push_back(bconfig(RIN, BCOidx...));
  bounds.push_back(bconfig(RMID, BCOidx...));

  auto adapted_dom(field.get_space().get_domain(adapted_dom_sol));

  // Coordinate centered point
  Kadath::Point center_pt(adapted_dom->get_center());
  double xc{center_pt(1)};

  // Initial radii starting at the inner adapted radius on the
  // equitorial plane since this should be the largest radius
  double r_init{get_radius(adapted_dom, INNER_EQUI)};
  double field_init{get_boundary_val(adapted_dom_sol, field)};

  // FIXME need a threashold here to ensure exp doesn't get unbounded
  double scale_fac = 2.0 - std::exp(field_init);

  double min_spacing_fac =
      (std::isnan(bconfig.set(BCO_PARAMS::MIN_SHELL_DR, BCOidx...)))
          ? 1.855
          : bconfig(BCO_PARAMS::MIN_SHELL_DR, BCOidx...);
  double r0 = r_init * scale_fac * min_spacing_fac;

  // We only add shells out to ROUT
  auto Rout{bconfig(ROUT, BCOidx...)};
  if (Rout / r0 > min_spacing_fac) {
    bounds.push_back(r0);

    // Setup points relative to coordinate center
    double r1{r0};
    Kadath::Point pt_r1(center_pt);
    pt_r1.set(1) = xc + r1;
    double field_r1 = field.val_point(pt_r1);
    // end point setup

    // some upper bound that should never be hit!
    auto max_shells =
        std::ceil((Rout - r_init) / (bco_utils::gold_ratio * r_init));

    for (auto i = 1; i < max_shells; ++i) {
      r0 = r1;
      double field_r0 = field_r1;

      //r1  = r0 * pow(min_spacing_fac, i);
      scale_fac = 2.0 - std::exp(field_r0);
      r1 = r0 * scale_fac * min_spacing_fac;
      pt_r1.set(1) = xc + r1;

      field_r1 = field.val_point(pt_r1);

      // stop looking for shells once \partial_r^2 \Psi
      // is roughly flat
      // Note: this will fail at local extrema!
      if (std::fabs(field_r1) <= 1e-5 || (Rout / r1 < bco_utils::gold_ratio)) {
        // cout << field_r1 << ", " << r1 << ", " << Rout / r1 << endl;
        break;
      }

      bounds.push_back(r1);

      if (i == max_shells - 1) {
        std::cerr << "Max iterations hit in setting CO bounds. Something went "
                     "very wrong!\n";
        std::cerr << "r0: " << r0 << '\n'
                  << "r1: " << r1 << '\n'
                  << "max_shells: " << max_shells << '\n';
      }
    }
  }
  // Add ROUT as the last bound
  bounds.push_back(Rout);
  bconfig(NSHELLS, BCOidx...) = bounds.size() - 3;
  return bounds;
}

/**
 * bco_utils::print_bounds
 *
 * helper function for printing domain boundaries in the readers - requires foreach compatible
 * container
 *
 * @tparam space_t type of numerical space
 * @param[input]  bounds: array of bounds
 * @param[input]  name: string with the name you want associated with the bounds (e.g. NS1)
 */
template <typename ary_t>
void print_bounds(std::string name, const ary_t& bary) {
  std::cout << name << ": ";
  for (auto& e : bary) {
    std::cout << e << " ";
  }
  std::cout << std::endl;
}

/**
 * bco_utils::mirr_from_mch
 *
 * helper function to compute Mirr from MCH
 *
 * @param[input]  chi: dimensionless spin parameter
 * @param[input]  mch: Christodoulou mass of the BH
 */
inline double mirr_from_mch(const double chi, const double mch) {
  return std::sqrt((1 + sqrt(1 - chi * chi)) / 2.) * mch;
}

/**
 * bco_utils::syst_mch
 *
 * helper function to compute MCH from a system of equations
 *
 * @tparam space_t: space type
 * @param[input]  syst: system of equations
 * @param[input]  space: numerical space
 * @param[input]  eq: equation string to lookup for spin equation (e.g. Sint1, Sint2)
 * @param[input]  dom: domain index to evalute the integral (i.e. inner homothetic)
 */
template <typename space_t>
double syst_mch(System_of_eqs& syst,
                const space_t& space,
                const std::string eq,
                const int dom) {
  Val_domain integS(syst.give_val_def(eq.c_str())()(dom));
  double S = space.get_domain(dom)->integ(integS, INNER_BC);

  Val_domain integMsq(syst.give_val_def("intMsq")()(dom));
  double Mirrsq = space.get_domain(dom)->integ(integMsq, INNER_BC);
  double Mirr = std::sqrt(Mirrsq);
  return std::sqrt(Mirrsq + S * S / 4. / Mirrsq);
}

/**
 * bco_utils::com_estimate
 *
 * helper function to compute the Newtonian COM as an initial estimate
 *
 * @param[input]  distance: separation distance
 * @param[input]  M1: mass of bco1
 * @param[input]  M2: mass of bco2
 */
inline double com_estimate(const double distance,
                           const double M1,
                           const double M2) {
  double half_dist = distance / 2.;
  double mass_sum = M1 + M2;
  return half_dist * (M1 - M2) / mass_sum;
}

/**
 * bco_utils::update_config_NS_radii
 *
 * helper function to update NS radii in the config file
 *
 * @tparam space_t: space type
 * @tparam config_t: Configurator type
 * @tparam Idx: parameter pack for Configurator secondary index.
 * @param[input]  space: numerical space
 * @param[input]  bconfig: Configurator file reference
 * @param[input]  dom: domain to pull min/max radii from
 * @param[input]  idx: Either empty or BCO1/BCO2
 */
template <typename space_t, typename config_t, typename... Idx>
void update_config_NS_radii(space_t& space,
                            config_t& bconfig,
                            const size_t dom,
                            Idx... idx) {

  auto [r_min, r_max] = bco_utils::get_rmin_rmax(space, dom);
  bconfig.set(RIN, idx...) = 0.5 * r_min;
  bconfig.set(RMID, idx...) = r_max;
  bconfig.set(ROUT, idx...) = gold_ratio * r_max;

  if (std::isnan(bconfig.set(BCO_PARAMS::MIN_SHELL_DR, idx...))) {
    bconfig.set(BCO_PARAMS::MIN_SHELL_DR, idx...) = 1.855;
  }
}

/**
 * bco_utils::update_config_BH_radii
 *
 * helper function to update BH radii in the config file
 *
 * @tparam space_t: space type
 * @tparam config_t: Configurator type
 * @tparam Idx: parameter pack for Configurator secondary index.
 * @param[input]  space: numerical space
 * @param[input]  bconfig: Configurator file reference
 * @param[input]  dom: domain to pull min/max radii from
 * @param[input]  conf: conformal factor scalar field
 * @param[input]  idx: Either empty or BCO1/BCO2
 */
template <typename space_t, typename config_t, typename... Idx>
void update_config_BH_radii(space_t& space,
                            config_t& bconfig,
                            const size_t dom,
                            const Scalar& conf,
                            Idx... idx) {

  auto [rmin, trmax] = bco_utils::get_rmin_rmax(space, dom);

  // estimate how small the inner radius should be based on relation
  // between conformal factor and numerical radius.
  // see https://arxiv.org/pdf/0805.4192, eq(64)
  double conf_inner = bco_utils::get_boundary_val(dom + 1, conf, INNER_BC);
  double conf_i_sq = conf_inner * conf_inner;
  double est_r_div2 = bconfig(MCH, idx...) / conf_i_sq;
  bconfig.set(RIN, idx...) = est_r_div2;

  // update config RMID based on AH Surface radius
  bco_utils::set_radius(dom, space, bconfig, RMID, idx...);
}

/**
 * bco_utils::set_decay
 *
 * helper function to initialize decay weight when superimposing two
 * BCOs.  Default is set to their separation distance.
 *
 * @tparam config_t: configurator file type
 * @param[input] bconfig: binary configurator file
 * @param[input] bco: BCO configurator index
 * @return: returns 1/weight^4
 */
template <typename config_t>
double set_decay(config_t& bconfig, const size_t bco) {
  if (std::isnan(bconfig.set(DECAY, bco)))
    bconfig.set(DECAY, bco) = bconfig(DIST) / 2.;
  const double weight4 = std::pow(bconfig(DECAY, bco), 4.);
  return 1. / weight4;
}

/**
 * bco_utils::compute_kerr_mirr
 *
 * Compute analytical Kerr irreducible mass
 *
 * @tparam config_t: configurator file type
 * @tparam coIdx: optional index of the BH in the binary config file
 * @param[input] bconfig: configurator file
 * @param[input] bco: BCO configurator index
 * @return: returns Mirr
 */
template <typename config_t, typename... coIdx>
double compute_kerr_mirr(config_t& bconfig, coIdx... bco) {
  // S := CHI * MCH^2
  const double MCHsq = bconfig(MCH, bco...) * bconfig(MCH, bco...);
  const double S = bconfig(CHI, bco...) * MCHsq;
  const double Mirrsq = (MCHsq + std::sqrt(MCHsq * MCHsq - S * S)) / 2;
  return std::sqrt(Mirrsq);
}

/**
 * bco_utils::get_bound_filled_field
 *
 * Fills each domain of a scalar field with the value at the respective
 * boundary
 *
 * @param[input] field: input scalar field
 * @param[input] bound: boundary to evaluate (e.g. INNER_BC, OUTER_BC)
 * @return: returns a Scalar field
 */
inline Scalar get_bound_filled_field(Scalar const& field, int const bound) {
  int const ndom = field.get_nbr_domains();
  Scalar out(field, false);
  for (auto d = 0; d < ndom; ++d) {
    auto npts = field(d).get_conf().get_dimensions();
    Index pos(npts);
    Index pos_b(pos);
    switch (bound) {
      case INNER_BC:
        break;
      case OUTER_BC:
        pos_b.set(1) = npts(1) - 1;
        break;
    }
    do {
      pos_b.set(2) = pos(2);
      pos_b.set(3) = pos(3);
      out.set_domain(d).set(pos) = field(d)(pos_b);
    } while (pos.inc());
    out.set_domain(d).set_base() = field(d).get_base();
  }
  return out;
}

/**
 * bco_utils::get_bound_filled_field_from_one_dom
 *
 * Fills all domains of a scalar field with the value at the respective
 * boundary from the specified domain only
 *
 * @param[input] field: input scalar field
 * @param[input] bound: boundary to evaluate (e.g. INNER_BC, OUTER_BC)
 * @param[input] dom: domain to copy to all other domains
 * @return: returns a Scalar field
 */
inline Scalar get_bound_filled_field_from_one_dom(Scalar const& field,
                                                  int const bound,
                                                  int const dom) {
  auto bound_field(get_bound_filled_field(field, bound));
  int const ndom = field.get_nbr_domains();
  Scalar out(bound_field);
  for (auto d = 0; d < ndom; ++d) {
    out.set_domain(d) = bound_field(dom);
    out.set_domain(d).set_base() = bound_field(dom).get_base();
  }
  return out;
}

/**
 * bco_utils::print_bounds_from_space
 *
 * Utility to print boundaries to cout from space
 *
 * @tparam space_t: numerical space type
 * @param[input] space: numerical space
 * @param[input] bound: boundary to evaluate (e.g. INNER_BC, OUTER_BC)
 */
template <class space_t>
void print_bounds_from_space(space_t const& space, int bound = OUTER_BC) {
  for (int i = 0; i < space.get_nbr_domains(); ++i)
    std::cout << bco_utils::get_radius(space.get_domain(i), bound) << " ";
  std::cout << std::endl;
};

/**
 * bco_utils::print_constant_space_resolution
 *
 * Utility to print resolution of the input space assuming
 * constant resolution amongst all domains
 *
 * @tparam space_t: numerical space type
 * @param[input] space: numerical space
 */
template <class space_t>
void print_constant_space_resolution(space_t const& space) {
  auto dom = space.get_domain(0);
  auto ndim = dom->get_ndim();
  std::array<std::string, 3> directions{"r", "theta", "phi"};
  for (auto i = 0; i < ndim; ++i)
    std::cout << dom->get_nbr_points()(i) << " (" << directions[i] << ")     ";
  std::cout << "\n";
};

/**
 * @brief Determine the next resolution that is reasonable to use with FFTW3
 *
 * @param res Current resolution
 * @return int Next resolution
 */
inline int next_resolution(int const res) {
  std::vector<int> ress{9, 11, 13, 17, 21, 25, 33};
  auto res_it = std::find_if(ress.begin(), ress.end(),
                             [&res](auto& n) { return n == res; });
  if (res_it == ress.end()) {
    std::stringstream msg;
    msg << "Resolution " << res << " not found in approved list.\n";
    std::__throw_runtime_error(msg.str().c_str());
  } else if (*res_it == ress.back()) {
    std::stringstream msg;
    msg << "Resolution cannot be increased beyond "<< ress.back() << " .\n";
    std::__throw_runtime_error(msg.str().c_str());
  }
  res_it++;
  return *res_it;
}

/** @}*/
}  // namespace bco_utils
}  // namespace Kadath
