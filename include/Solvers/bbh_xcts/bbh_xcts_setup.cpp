#pragma once
#include "mpi.h"
#include "bbh_xcts_solver.hpp"
#include<array>
#include<string>

/**
 * \addtogroup BBH_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
namespace bco_u = ::Kadath::bco_utils;
/**
 * @brief 
 * 
 * @tparam config_t 
 * @param bconfig 
 */

template<class config_t>
void bbh_xcts_setup_bin_config(config_t& bconfig){
  check_dist(bconfig(BIN_PARAMS::DIST), 
    bconfig(BCO_PARAMS::MCH, NODES::BCO1), bconfig(BCO_PARAMS::MCH, NODES::BCO2));

  // Binary Parameters
  bconfig.set(BIN_PARAMS::REXT) = 2 * bconfig(BIN_PARAMS::DIST);
  bconfig.set(BIN_PARAMS::Q) = bconfig(BCO_PARAMS::MCH, NODES::BCO2) 
                             / bconfig(BCO_PARAMS::MCH, NODES::BCO1);
  
  // classical Newtonian estimate
  bconfig.set(BIN_PARAMS::COM) = bco_u::com_estimate(bconfig(BIN_PARAMS::DIST), 
    bconfig(BCO_PARAMS::MCH, NODES::BCO1), bconfig(BCO_PARAMS::MCH, NODES::BCO2));
  
  // obtain 3PN estimate for the global, orbital omega
  bco_u::KadathPNOrbitalParams(bconfig, \
        bconfig(BCO_PARAMS::MCH, NODES::BCO1), bconfig(BCO_PARAMS::MCH, NODES::BCO2));

  // delete ADOT, this can always be recalculated during
  // the eccentricity reduction stage
  bconfig.reset(BIN_PARAMS::ADOT);
}

template<class config_t>
void bbh_xcts_setup_space (config_t& bconfig) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<int,2> bcos{NODES::BCO1, NODES::BCO2};
  std::array<std::string, 2> filenames;

  for(int i = 0; i < 2; ++i)
    filenames[i] = solve_BH_from_binary(bconfig, bcos[i]);  

  if(rank == 0)
    bbh_xcts_superimposed_import(bconfig, filenames);
  MPI_Barrier(MPI_COMM_WORLD);

  bconfig.open_config();
  bconfig.control(CONTROLS::SEQUENCES) = false;
}

template<class config_t>
void bbh_xcts_superimposed_import(config_t& bconfig, std::array<std::string, 2> BHfilenames) {
  // load single BH configuration
  std::string bh1filename{BHfilenames[0]};
  kadath_config_boost<BCO_BH_INFO> BH1config(bh1filename);
  
  std::string bh2filename{BHfilenames[1]};
  kadath_config_boost<BCO_BH_INFO> BH2config(bh2filename);

  bbh_xcts_setup_boosted_3d(BH1config, BH2config, bconfig);
}

inline void bbh_xcts_setup_boosted_3d(
  kadath_config_boost<BCO_BH_INFO>& BH1config, 
  kadath_config_boost<BCO_BH_INFO>& BH2config,
  kadath_config_boost<BIN_INFO>& bconfig) {
  
  std::string in_spacefile = BH1config.space_filename();
  FILE *ff1 = fopen(in_spacefile.c_str(), "r");
  Space_adapted_bh spacein1(ff1) ;
  Scalar       confin1 (spacein1, ff1) ;
	Scalar       lapsein1(spacein1, ff1) ;
	Vector       shiftin1(spacein1, ff1) ;
  fclose(ff1) ;
  bco_u::update_config_BH_radii(spacein1, BH1config, 1, confin1);
  
  in_spacefile = BH2config.space_filename();
  FILE *ff2 = fopen(in_spacefile.c_str(), "r");
  Space_adapted_bh spacein2(ff2) ;
  Scalar       confin2 (spacein2, ff2) ;
	Scalar       lapsein2(spacein2, ff2) ;
	Vector       shiftin2(spacein2, ff2) ;
  fclose(ff2) ;
  bco_u::update_config_BH_radii(spacein2, BH2config, 1, confin2);
  // end opening old solutions
  
  // update binary parameters - but save shell input
  const int nshells1 = bconfig(BCO_PARAMS::NSHELLS, NODES::BCO1);
  const int nshells2 = bconfig(BCO_PARAMS::NSHELLS, NODES::BCO2);
  for(int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i) {
    bconfig.set(i, NODES::BCO1) = BH1config.set(i) ;
    bconfig.set(i, NODES::BCO2) = BH2config.set(i) ;
  }
  bconfig.set(BCO_PARAMS::NSHELLS, NODES::BCO1) = nshells1;
  bconfig.set(BCO_PARAMS::NSHELLS, NODES::BCO2) = nshells2;
  
  const double r_max_tot = (bconfig(BCO_PARAMS::RMID, NODES::BCO1) > bconfig(RMID, NODES::BCO2)) ? \
    bconfig(BCO_PARAMS::RMID, NODES::BCO1) : bconfig(BCO_PARAMS::RMID, NODES::BCO2);
  
  const double rout_sep_est = (bconfig(BIN_PARAMS::DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  const double rout_max_est = 2. * bco_u::gold_ratio * r_max_tot;
  
  bconfig.set(BCO_PARAMS::ROUT, NODES::BCO1) = (rout_sep_est > rout_max_est) ? rout_max_est : rout_sep_est;
  bconfig.set(BCO_PARAMS::ROUT, NODES::BCO2) = bconfig(BCO_PARAMS::ROUT, NODES::BCO1);

  std::vector<double> out_bounds(1+bconfig(BIN_PARAMS::OUTER_SHELLS));
std::vector<double> BH1_bounds;
  {
    auto ddrPsi(compute_ddrPsi(
      spacein1, 
      confin1, 
      Metric_flat(spacein1, shiftin1.get_basis()), 
      {0,1}
    ));
    BH1_bounds = bco_u::set_arb_bounds(bconfig, BCO1, ddrPsi, 2, 0.9);
  }
  std::cout << "Bound1 done\n";
  bco_u::print_bounds("bh1",BH1_bounds);
  std::vector<double> BH2_bounds;
  {
    auto ddrPsi(compute_ddrPsi(
      spacein2, 
      confin2, 
      Metric_flat(spacein2, shiftin2.get_basis()), 
      {0,1}
    ));
    BH2_bounds = bco_u::set_arb_bounds(bconfig, BCO2, ddrPsi, 2, 0.9);
  }
  std::cout << "Bound2 done\n";
  bco_u::print_bounds("bh2",BH2_bounds);

  // for out_bounds.size > 1 - add equi-distant shells
  for(int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(BIN_PARAMS::REXT) + e * 0.25 * bconfig(BIN_PARAMS::REXT);

  bco_u::set_BH_bounds(BH1_bounds, bconfig, NODES::BCO1, true);
  bco_u::set_BH_bounds(BH2_bounds, bconfig, NODES::BCO2, false);

  // create space containing the domain decomposition
  int type_coloc = CHEB_TYPE;
  Space_bin_bh space(type_coloc, bconfig(BIN_PARAMS::DIST), BH1_bounds, BH2_bounds, out_bounds, bconfig(BIN_PARAMS::BIN_RES));
  Base_tensor basis(space, CARTESIAN_BASIS);

  auto interp_BH_fields = [&](auto& bhconf, auto& bhlapse, auto& bhshift, auto& bhspacein) {
    const Domain_shell_outer_homothetic* old_bh_outer = 
      dynamic_cast<const Domain_shell_outer_homothetic*>(bhspacein.get_domain(1));
    //Update BH fields based to help with interpolation later
    bco_u::update_adapted_field(bhconf, 2, 1, old_bh_outer, OUTER_BC);
    bco_u::update_adapted_field(bhlapse, 2, 1, old_bh_outer, OUTER_BC);
    for(int i = 1; i <=3; ++i)
      bco_u::update_adapted_field(bhshift.set(i), 2, 1, old_bh_outer, OUTER_BC);

  };
  
  interp_BH_fields(confin1, lapsein1, shiftin1, spacein1);
  interp_BH_fields(confin2, lapsein2, shiftin2, spacein2);
  
  double xc1 = bco_u::get_center(space, space.BH1);
  double xc2 = bco_u::get_center(space, space.BH2);
  
  Scalar conf(space);
  conf.annule_hard();

  Scalar lapse(space);
  lapse.annule_hard();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();
  
  const double bh1_invw4 = bco_u::set_decay(bconfig, NODES::BCO1);
  const double bh2_invw4 = bco_u::set_decay(bconfig, NODES::BCO2);
  
  // start importing the fields from the single star
  int ndom = space.get_nbr_domains();
  for(int dom = 0; dom < ndom; dom++)
  {
    // get an index in each domain to iterate over all colocation points
		Index new_pos(space.get_domain(dom)->get_nbr_points());

		do {
      if(dom <= space.BH1+1 || dom == space.BH2 || dom == space.BH2+1)
        continue;
      // get cartesian coordinates of the current colocation point
			double x = space.get_domain(dom)->get_cart(1)(new_pos);
			double y = space.get_domain(dom)->get_cart(2)(new_pos);
			double z = space.get_domain(dom)->get_cart(3)(new_pos);

      // define a point shifted suitably to the stellar centers in the binary
			Point absol1(3);
			absol1.set(1) = (x - xc1);
			absol1.set(2) = y;
			absol1.set(3) = z;
      double r2 = y * y + z * z; 
      double r2_1 = (x - xc1) * (x - xc1) + r2;
      double r4_1  = r2_1 * r2_1;
      double r4_invw4_1 = r4_1 * bh1_invw4;
      double decay_1 = std::exp(-r4_invw4_1);

      Point absol2(3);
			absol2.set(1) = (x - xc2);
			absol2.set(2) = y;
			absol2.set(3) = z;
      double r2_2 = (x - xc2) * (x - xc2) + r2;
      double r4_2  = r2_2 * r2_2;
      double r4_invw4_2 = r4_2 * bh2_invw4;
      double decay_2 = std::exp(-r4_invw4_2);

      if (dom < ndom - 1) {
        conf .set_domain(dom).set(new_pos) = 1. + decay_1 * (confin1.val_point(absol1) - 1.) \
                                           + decay_2 * (confin2.val_point(absol2) - 1.);
        lapse.set_domain(dom).set(new_pos) = 1. + decay_1 * (lapsein1.val_point(absol1) - 1.) \
                                           + decay_2 * (lapsein2.val_point(absol2) - 1.);
        for (int i = 1; i <= 3; i++)
          shift.set(i).set_domain(dom).set(new_pos) = decay_1 * shiftin1(i).val_point(absol1) \
                                                    + decay_2 * shiftin2(i).val_point(absol2);
   
      } else {
        // We have to set the compactified domain manually since the outer collocation point is always
        // at inf which is undefined numerically
				conf .set_domain(dom).set(new_pos) = 1.;
				lapse.set_domain(dom).set(new_pos) = 1.;
      }
      // loop over all colocation points
		} while(new_pos.inc());
	} //end importing fields
  
  // safety to ensure excision region is empty
  auto clear_excision_region = [&] (int nuc) {
    for(int d = nuc; d < nuc+2; ++d) {
      conf.set_domain(d).annule_hard();
      lapse.set_domain(d).annule_hard();
      for (int i = 1; i <= 3; i++)
        shift.set(i).set_domain(d).annule_hard();
    }
  };
  clear_excision_region(space.BH1);
  clear_excision_region(space.BH2);
    
	conf.std_base();
	lapse.std_base();
  shift.std_base();

  bco_u::save_to_file(space, bconfig, conf, lapse, shift);
}
/** @}*/
}}