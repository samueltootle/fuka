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

// FUKA includes
#include "Configurator/config_bco.hpp"
#include "Solvers/solver_startup.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "EOS/standalone/tov.hh"

// Kadath includes
#include "kadath.hpp"
#include "kadath_adapted_polar.hpp"
#include "kadath_adapted.hpp"

// C++ includes
#include <fstream>
#include <string>
#include <type_traits>
#include <memory>
#include <algorithm>

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Solvers;

template<class eos_t, class config_t>
int norot_2dsetup (config_t& bconfig);

template<class config_t>
void set_2dfields(config_t& bconfig);

int main(int argc, char** argv) {
  int rank = 0;
  using config_t = kadath_config_boost<BCO_NS_INFO>;
  using InitSolver = Initialize_Solver<config_t>;
  
  // Initialize static member variables
  InitSolver::input_configname = "initial_ns2d.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;
  
  if(InitSolver::example_setup) {
    if(rank == 0) {
      // Generate <example name>.info and terminate      
      bconfig.set_defaults();
      bconfig.control(CONTROLS::SEQUENCES) = InitSolver::setup_first;
      bconfig.set(BCO_PARAMS::DIM) = 2;
      set_2dfields(bconfig);
      bconfig.set(BCO_PARAMS::DIFF_LAWQ) = 1;
      bconfig.set(BCO_PARAMS::DIFF_ARATIO) = 1.;
      bconfig.set(BCO_PARAMS::DIFF_RRATIO) = 0.875;
      bconfig.set_stage(STAGES::DIFF_ROT) = true;
      bconfig.write_config();
    }
  } else {
    bconfig.open_config();

    /* Loading of the configuration parameters.
    *
    * See $HOME_KADATH/include/Configurator/config_enums.hpp for
    * possible BCO_PARAMS and EOS_PARAMS index names 
    */
    const double h_cut = bconfig.eos<double>(HCUT);
    const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      bconfig.set(HC) = EOS<eos_t,PRESSURE>::h_cold__rho(bconfig(NC));
      
      // call setup routine, setting up the actual numerical domains
      int res = norot_2dsetup<eos_t>(bconfig);
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      bconfig.set(HC) = EOS<eos_t,PRESSURE>::h_cold__rho(bconfig(NC));
      
      // call setup routine, setting up the actual numerical domains
      int res = norot_2dsetup<eos_t>(bconfig);
    }
    else { 
      std::cerr << eos_type << " is not recognized.\n";
      std::_Exit(EXIT_FAILURE);
    }
  }

  return 0;
}


// Here we take a TOV solution and create a 1D interpolator for the relevant quantities
// as a function of the isotropic radius
template<typename tov_t>
auto setup_interpolation(tov_t& tov ) {
  const size_t max_iter = tov->state.size();
  const size_t start = 0; //number of iterations to skip
  const size_t num_pts = max_iter - start;
  
  std::unique_ptr<double[]> radius_lin_ptr{new double[num_pts]};
  std::unique_ptr<double[]> conf_lin_ptr{new double[num_pts]};
  std::unique_ptr<double[]> lapse_lin_ptr{new double[num_pts]};
  std::unique_ptr<double[]> rho_lin_ptr{new double[num_pts]};
  
  for(auto j=0; j < max_iter; ++j) {
    //lapse est.
    auto phi = tov->state[j][tov->PHI];
    
    lapse_lin_ptr[j] = std::exp(phi);
    
    conf_lin_ptr[j] = tov->state[j][tov->CONF];

    radius_lin_ptr[j] = tov->state[j][tov->RISO];

    rho_lin_ptr[j] = tov->state[j][tov->RHOB];
  }
  
  // linear interpolation of conf(r_isotropic), lapse(r_isotropic), and rho(r_isotropic)
  linear_interp_t<double,3> ltp(num_pts, std::move(radius_lin_ptr), std::move(lapse_lin_ptr), std::move(rho_lin_ptr), std::move(conf_lin_ptr)); 
  return ltp; 
}

template<class config_t>
void set_2dfields(config_t& bconfig) {
  for(auto i = 0; i < BCO_FIELDS::NUM_BCO_FIELDS; i++)
    bconfig.set_field(i) = false;

  bconfig.set_field(BCO_FIELDS::LAPSE) = true;
  bconfig.set_field(BCO_FIELDS::LAP_ATERM) = true;
  bconfig.set_field(BCO_FIELDS::LOGH)  = true;
}

template<class eos_t, class config_t>
int norot_2dsetup(config_t& bconfig) {

  enum ltpQ { LAPSE=0, RHO, CONF }; 
  // Find 1D TOV solution for a given ADM Mass
  auto tov = std::make_unique<MargheritaTOV<eos_t>>();
  tov->solve_for_MADM(bconfig(MADM));

  // update surface radius estimate
  bconfig.set(RMID) = tov->state.back()[tov->RISO];
  bconfig.set(RIN) = 0.5 * bconfig(RMID);
  bconfig.set(ROUT) = bco_utils::gold_ratio * bconfig(RMID);

  // interpolate TOV solution
  auto lintp = setup_interpolation(tov); 
  
  const int dim = 2;
  bconfig.set(DIM) = dim;

  int type_coloc = CHEB_TYPE;
  Dim_array res(dim);
  res.set(0) = bconfig(BCO_RES);
  res.set(1) = bconfig(BCO_RES);

  Point center(dim);
  for (int i = 1; i <= dim; i++)
    center.set(i) = 0;

  int ndom      = 4 + bconfig(NSHELLS);
  Array<double> bounds(ndom - 1);
  bounds.set(0) = bconfig(RIN);
  bounds.set(1) = bconfig(RMID);
  bounds.set(2) = bconfig(ROUT);

  for(int shell = 1, b = 3; shell <= bconfig(NSHELLS); ++shell, ++b) {
    bounds.set(b) = bounds(b-1) * 2;
  }

  Space_polar_adapted space(type_coloc, center, res, bounds);

  Scalar logh(space);
  logh.annule_hard();

  Scalar conf(space);
  conf = 1.;

  // same for the lapse
  Scalar lapse(conf);

  Scalar r_field(space);
  r_field.annule_hard();
  for (auto i = 0; i < ndom; ++i)
    r_field.set_domain(i) = space.get_domain(i)->get_radius();
  auto npts = space.get_domain(ndom-1)->get_nbr_points();
  Index pos_c(npts);
  pos_c.set(0) = npts(0) - 1; // outer_bc

  for(auto i = 0; i < npts(1); ++i) {
    pos_c.set(1) = i;
    r_field.set_domain(ndom-1).set(pos_c) = 1e10;
  }
  pos_c.set_start();

  // update the fields based on TOV solution for a given domain
  auto update_fields= [&](const size_t dom) {
    Index pos(space.get_domain(dom)->get_nbr_points());
    do {
      double rval = r_field(dom)(pos);
      auto all_ltp = lintp.interpolate_all(rval);
      auto h = EOS<eos_t,DENSITY>::h_cold__rho(all_ltp[ltpQ::RHO]);
      if(dom == 0 && pos(0) == 0 && pos(1) == 0)
        bconfig.set(HC) = h;
      logh.set_domain(dom).set(pos) = (std::log(h) <= 0) ? 0. : std::log(h); 
      lapse.set_domain(dom).set(pos) = all_ltp[ltpQ::LAPSE];
      conf.set_domain(dom).set(pos) = all_ltp[ltpQ::CONF];
    }while(pos.inc());
  };
  for(int i = 0; i < ndom; ++i)
    update_fields(i);

  for (int d = 2; d < ndom; d++)
    logh.set_domain(d).annule_hard();

  // Fix compactified domain metric variables
  conf.set_domain(ndom-1) = 1 + (conf(ndom-1)(pos_c) - 1) * r_field(ndom-1)(pos_c) / r_field(ndom-1);
  lapse.set_domain(ndom-1) = 1 + (lapse(ndom-1)(pos_c) - 1) * r_field(ndom-1)(pos_c) / r_field(ndom-1);

  Scalar A(conf * conf);
  Scalar nu(log(lapse));
  Scalar nulogA(nu + log(A));
  Scalar nulogB(nulogA); // non-rotating limit B = A
  
  nu.std_base();
  nulogA.std_base();
  nulogB.std_base();
  logh.std_base();

  set_2dfields(bconfig);
  bconfig.set_filename("initns");
  bco_utils::save_to_file(space, bconfig, nulogA, nu, logh, nulogB);
  return EXIT_SUCCESS;
}
