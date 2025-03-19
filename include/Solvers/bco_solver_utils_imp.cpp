/*
 * Copyright 2022
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author:
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
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
#include <iostream>
#include "Solvers/bh_3d_xcts/bh_3d_xcts_driver.hpp"
#include "bco_solver_utils.hpp"
#include "bco_utilities.hpp"
#include "bh_3d_xcts/bh_3d_xcts_solver.hpp"
#include "mpi.h"
#include "ns_3d_xcts/ns_3d_xcts_driver.hpp"
#include "ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "ns_3d_xcts/NS_ISO_to_XCTS_convert.hpp"
#include "ns_isotropic/ns_isotropic_driver.hpp"
#include "sequences/parameter_sequence.hpp"
#include "solvers.hpp"

/**
 * \addtogroup Solver_utils
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {

template <typename config_t>
std::string solve_NS_from_binary(config_t& bconfig, const size_t bco) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string output_path = (bconfig.control(CONTROLS::SAVE_COS))
                                ? get_cos_path()
                                : bconfig.config_outputdir() + "/COs";
  fs::create_directory(output_path);

  kadath_config_boost<BCO_NS_INFO> nsconfig;
  nsconfig.set_defaults();

  // Tells the NS driver to initialize the numerical space and fields
  nsconfig.control(CONTROLS::SEQUENCES) = true;

  // copy parameters from binary configuration
  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    nsconfig.set(i) = bconfig.set(i, bco);

  for (int i = 0; i < EOS_PARAMS::NUM_EOS_PARAMS; ++i)
    nsconfig.set_eos(i) = bconfig.set_eos(i, bco);

  for (int i = 0; i < CONTROLS::NUM_CONTROLS; ++i)
    nsconfig.control(i) = bconfig.control(i);

  for (int i = 0; i < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++i)
    nsconfig.seq_setting(i) = bconfig.seq_setting(i);

  Parameter_sequence resolution("res", BCO_PARAMS::BCO_RES);
  resolution.set(9, 9, nsconfig(BCO_PARAMS::BCO_RES));

  nsconfig.set_filename("initns");
  nsconfig.set_outputdir(output_path);

  if (bconfig.control(CONTROLS::USE_BOOSTED_CO))
    nsconfig.set_stage(STAGES::BIN_BOOST) = true;

  const int ninshells = bconfig(BCO_PARAMS::NINSHELLS, bco);
  const int nshells = bconfig(BCO_PARAMS::NSHELLS, bco);

  nsconfig(BCO_PARAMS::NINSHELLS) = 0;
  nsconfig(BCO_PARAMS::NSHELLS) =
      (bconfig.control(CONTROLS::CO_USE_SHELLS)) ? nshells : 0;

  int err = ns_3d_xcts_binary_boost_driver(nsconfig, resolution, output_path,
                                           bconfig, bco);
  MPI_Barrier(MPI_COMM_WORLD);

  // update binary parameters
  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    bconfig.set(i, bco) = nsconfig.set(i);
  // put shells back for the binary
  bconfig.set(BCO_PARAMS::NINSHELLS, bco) = ninshells;
  bconfig.set(BCO_PARAMS::NSHELLS, bco) = nshells;

  return nsconfig.config_filename_abs();
}

template <typename config_t>
std::string solve_BH_from_binary(config_t& bconfig, const size_t bco) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string output_path = (bconfig.control(CONTROLS::SAVE_COS))
                                ? get_cos_path()
                                : bconfig.config_outputdir() + "/COs";
  fs::create_directory(output_path);

  kadath_config_boost<BCO_BH_INFO> bhconfig;
  bhconfig.set_defaults();

  // Tells the NS driver to initialize the numerical space and fields
  bhconfig.control(CONTROLS::SEQUENCES) = true;

  // copy parameters from binary configuration
  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    bhconfig.set(i) = bconfig.set(i, bco);

  for (int i = 0; i < CONTROLS::NUM_CONTROLS; ++i)
    bhconfig.control(i) = bconfig.control(i);

  for (int i = 0; i < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++i)
    bhconfig.seq_setting(i) = bconfig.seq_setting(i);

  Parameter_sequence resolution("res", BCO_PARAMS::BCO_RES);
  resolution.set(9, 9, bhconfig(BCO_PARAMS::BCO_RES));

  bhconfig.set_filename("initbh");
  bhconfig.set_outputdir(output_path);

  if (bconfig.control(CONTROLS::USE_BOOSTED_CO))
    bhconfig.set_stage(STAGES::BIN_BOOST) = true;

  const int nshells = bconfig(NSHELLS, bco);
  bhconfig(BCO_PARAMS::NSHELLS) =
      (bconfig.control(CONTROLS::CO_USE_SHELLS)) ? nshells : 0;

  int err = bh_3d_xcts_binary_boost_driver(bhconfig, resolution, output_path,
                                           bconfig, bco);
  MPI_Barrier(MPI_COMM_WORLD);

  // update binary parameters
  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    bconfig.set(i, bco) = bhconfig.set(i);
  bconfig.set(BCO_PARAMS::NSHELLS, bco) = nshells;

  MPI_Barrier(MPI_COMM_WORLD);
  return bhconfig.config_filename_abs();
}

inline void check_dist(double dist,
                       double M1,
                       double M2,
                       double garbage_factor) {
  auto M = M1 + M2;
  auto const garbage_dist = garbage_factor * M;
  auto const recommended_dist = 8. * M;
  if (dist <= garbage_dist) {
    std::cerr << "Distance is set to (" << dist
              << ") which may not give results. \nSet to (" << recommended_dist
              << ") for something reasonable.\n";
    std::_Exit(EXIT_FAILURE);
  }
}

template <class space_t, class metric_t>
Scalar compute_ddrPsi(space_t& space,
                      Scalar& conf,
                      metric_t metric,
                      std::vector<int> excluded_doms,
                      int bound_dom = -1) {
  std::string const def_drP{"drP = dr(P)"};
  std::string const def_drdrP{"ddrP = dr(drP)"};
  auto ndom{space.get_nbr_domains()};
  System_of_eqs syst(space);
  metric.set_system(syst, "f");
  syst.add_cst("P", conf);

  if (excluded_doms.size() == 0) {
    syst.add_def(def_drP.c_str());
    syst.add_def(def_drdrP.c_str());
  } else {
    auto last_dom = (bound_dom == -1) ? ndom : bound_dom;
    for (auto dom = 0; dom < last_dom; ++dom) {
      auto res = std::find(excluded_doms.begin(), excluded_doms.end(), dom);
      if (res == std::end(excluded_doms)) {
        syst.add_def(dom, def_drP.c_str());
        syst.add_def(dom, def_drdrP.c_str());
      }
    }
  }
  Scalar field(syst.give_val_def("ddrP"));
  field.std_base();
  return field;
}

template <class space_t, class metric_t>
Scalar compute_drPsi(space_t& space,
                     Scalar& conf,
                     metric_t metric,
                     std::vector<int> excluded_doms,
                     int bound_dom = -1) {
  std::string const def_drP{"drP = dr(P)"};
  auto ndom{space.get_nbr_domains()};
  System_of_eqs syst(space);
  metric.set_system(syst, "f");
  syst.add_cst("P", conf);

  if (excluded_doms.size() == 0) {
    syst.add_def(def_drP.c_str());
  } else {
    auto last_dom = (bound_dom == -1) ? ndom : bound_dom;
    for (auto dom = 0; dom < last_dom; ++dom) {
      auto res = std::find(excluded_doms.begin(), excluded_doms.end(), dom);
      if (res == std::end(excluded_doms)) {
        syst.add_def(dom, def_drP.c_str());
      }
    }
  }
  Scalar field(syst.give_val_def("drP"));
  field.std_base();
  return field;
}

template <typename config_t>
std::string solve_NS_ISO_from_XCTS_config(config_t& bconfig, ns_sequence const & seq) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string output_path = (bconfig.control(CONTROLS::SAVE_COS))
                                ? get_cos_path()
                                : bconfig.config_outputdir() + "/COs";
  fs::create_directory(output_path);

  kadath_config_boost<BCO_ISO_NS_INFO> nsconfig;
  nsconfig.set_defaults();

  // Tells the NS driver to initialize the numerical space and fields
  nsconfig.control(CONTROLS::SEQUENCES) = true;

  // Activate all controls that are active in the seq config
  for (int i = 0; i < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++i)
    nsconfig.seq_setting(i) = bconfig.seq_setting(i);

  // copy parameters from binary configuration
  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    nsconfig.set(i) = bconfig.set(i);
  nsconfig.set(BCO_PARAMS::DIM) = 2;

  for (int i = 0; i < CONTROLS::NUM_CONTROLS; ++i)
    nsconfig.control(i) = bconfig.control(i);

  update_eos_parameters(bconfig, nsconfig);
  update_diffrot_parameters(bconfig, nsconfig);

  Parameter_sequence resolution("res", BCO_PARAMS::BCO_RES);
  resolution.set(9, 9, nsconfig(BCO_PARAMS::BCO_RES));

  nsconfig.set_filename("initns");
  nsconfig.set_outputdir(output_path);

  if (bconfig.control(CONTROLS::USE_BOOSTED_CO))
    nsconfig.set_stage(STAGES::BIN_BOOST) = true;

  auto ns_iso_sol_config = ns_isotropic_sequence(nsconfig, seq, resolution, output_path);

  const std::string eos_type = bconfig.template eos<std::string>(EOSTYPE);
  EOS_Function_Dispatcher::dispatch<NS_ISO_to_XCTS_convert>(ns_iso_sol_config, eos_type,
                                                            ns_iso_sol_config, output_path);
  bconfig.set_filename("initns_xcts.info");
  bconfig.set_outputdir(output_path);
  MPI_Barrier(MPI_COMM_WORLD);
  bconfig.open_config();

  return bconfig.config_filename_abs();
}
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath