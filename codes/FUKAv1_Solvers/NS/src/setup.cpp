/*
 * Copyright 2021
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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
 *
 */
#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <type_traits>
#include "Configurator/config_bco.hpp"
#include "EOS/EOS.hh"
#include "EOS/standalone/tov.hh"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_polar.hpp"
using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

template <typename eos_t, typename config_t>
int norot_3dsetup(config_t bconfig);

int main(int argc, char** argv) {
  // class to load the configuration
  kadath_config_boost<BCO_NS_INFO> bconfig;

  // check if configuration is given
  // or if a standard setup should be generated
  if (argc < 2) {
    std::cout << "Using default Togashi setup.\n";
    bconfig.set_defaults();
  } else {
    std::string input_filename = std::string{argv[1]};
    std::cout << "Creating setup from: " << input_filename << "\n";

    bconfig.set_filename(input_filename);
    bconfig.open_config();
  }

  // set new base filename for new setup
  std::string base_filename = "./initns.info";
  bconfig.set_filename(base_filename);

  /* Loading of the configuration parameters.
   *
   * See $HOME_KADATH/include/Configurator/config_enums.hpp for
   * possible BCO_PARAMS and EOS_PARAMS index names 
   */
  const double h_cut = bconfig.eos<double>(HCUT);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

  if (eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut);
    bconfig.set(HC) = EOS<eos_t, PRESSURE>::h_cold__rho(bconfig(NC));

    // call setup routine, setting up the actual numerical domains
    int res = norot_3dsetup<eos_t>(bconfig);
  } else if (eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0)
                               ? 2000
                               : bconfig.eos<int>(INTERP_PTS);

    EOS<eos_t, PRESSURE>::init(eos_file, h_cut, interp_pts);
    bconfig.set(HC) = EOS<eos_t, PRESSURE>::h_cold__rho(bconfig(NC));

    // call setup routine, setting up the actual numerical domains
    int res = norot_3dsetup<eos_t>(bconfig);
  } else {
    std::cerr << eos_type << " is not recognized.\n";
    std::_Exit(EXIT_FAILURE);
  }

  std::cout << bconfig << std::endl;

  return 0;
}

// Here we take a TOV solution and create a 1D interpolator for the relevant quantities
// as a function of the isotropic radius
template <typename tov_t>
auto setup_interpolation(tov_t& tov) {
  const size_t max_iter = tov->state.size();
  const size_t start = 0;  //number of iterations to skip
  const size_t num_pts = max_iter - start;

  std::unique_ptr<double[]> radius_lin_ptr{new double[num_pts]};
  std::unique_ptr<double[]> conf_lin_ptr{new double[num_pts]};
  std::unique_ptr<double[]> lapse_lin_ptr{new double[num_pts]};
  std::unique_ptr<double[]> rho_lin_ptr{new double[num_pts]};

  for (auto j = 0; j < max_iter; ++j) {
    //lapse est.
    auto phi = tov->state[j][tov->PHI];

    lapse_lin_ptr[j] = std::exp(phi);

    conf_lin_ptr[j] = tov->state[j][tov->CONF];

    radius_lin_ptr[j] = tov->state[j][tov->RISO];

    rho_lin_ptr[j] = tov->state[j][tov->RHOB];
  }

  // linear interpolation of conf(r_isotropic), lapse(r_isotropic), and rho(r_isotropic)
  linear_interp_t<double, 3> ltp(num_pts, std::move(radius_lin_ptr),
                                 std::move(lapse_lin_ptr),
                                 std::move(rho_lin_ptr),
                                 std::move(conf_lin_ptr));
  return ltp;
}

// this function takes a configuration and handles creating the corresponding numerical setup
template <typename eos_t, typename config_t>
int norot_3dsetup(config_t bconfig) {
  enum ltpQ { LAPSE = 0, RHO, CONF };

  // Find 1D TOV solution for a given ADM Mass
  auto tov = std::make_unique<MargheritaTOV<eos_t>>();
  tov->solve_for_MADM(bconfig(MADM));

  // update surface radius estimate
  bconfig.set(RMID) = tov->state.back()[tov->RISO];

  // interpolate TOV solution
  auto lintp = setup_interpolation(tov);

  // number of dimensions
  const int dim = bconfig(DIM);

  // collocation point type and number of collocation points per domain
  int type_coloc = CHEB_TYPE;
  Dim_array res(dim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = bconfig(BCO_RES);
  res.set(2) = bconfig(BCO_RES) - 1;

  // physical center of the system, arbitrary
  Point center(dim);
  for (int i = 1; i <= dim; i++)
    center.set(i) = 0;

  if (!bconfig.control(USE_CONFIG_VARS)) {
    bconfig.set(RIN) = 0.5 * bconfig(RMID);
    bconfig.set(ROUT) = 1.5 * bconfig(RMID);
  }

  // domain boundaries, i.e. radii for each domain of the single star space
  int ndom = 4;
  Array<double> bounds(ndom - 1);
  bounds.set(0) = bconfig(RIN);
  bounds.set(1) = bconfig(RMID);
  bounds.set(2) = bconfig(ROUT);

  // generate a full single star space including compactification to infinity
  Space_spheric_adapted space(type_coloc, center, res, bounds);

  // use a basis of cartesian type
  Base_tensor basis(space, CARTESIAN_BASIS);

  // generate a radius field for use with the 1D interpolators
  CoordFields<Space_spheric_adapted> cfg(space);
  Scalar r_field(cfg.radius());

  // start to create the final fields

  // H = log(h), the logarithm of the enthalpy
  Scalar logh(space);
  logh.annule_hard();

  // conformal factor is initialized to one, representing flat space as an initial guess
  Scalar conf(space);
  conf = 1.;

  // same for the lapse
  Scalar lapse(conf);

  // update the fields based on TOV solution for a given domain
  auto update_fields = [&](const size_t dom) {
    Index pos(space.get_domain(dom)->get_nbr_points());
    do {
      double rval = r_field(dom)(pos);
      auto all_ltp = lintp.interpolate_all(rval);
      auto h = EOS<eos_t, DENSITY>::h_cold__rho(all_ltp[ltpQ::RHO]);
      if (pos(0) == 0 && pos(1) == 0 && pos(2) == 0)
        bconfig.set(HC) = h;
      logh.set_domain(dom).set(pos) = (std::log(h) < 0) ? 0. : std::log(h);
      lapse.set_domain(dom).set(pos) = all_ltp[ltpQ::LAPSE];
      conf.set_domain(dom).set(pos) = all_ltp[ltpQ::CONF];
    } while (pos.inc());
  };
  for (int i = 0; i < 3; ++i)
    update_fields(i);

  for (int d = 2; d < ndom; ++d)
    logh.set_domain(d).annule_hard();

  // shift is zero in the beginning and in general for TOV solutions
  Vector shift(space, CON, basis);
  shift.annule_hard();
  shift.std_base();
  logh.std_base();
  conf.std_base();
  lapse.std_base();
  // finished create the fields

  // write space and config to files on one processor
  bconfig.set_stage(PRE) = false;
  bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh);
  return EXIT_SUCCESS;
}
