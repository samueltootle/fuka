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
#pragma once
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "coord_fields.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "name_tools.hpp"
#include <cstdlib>
#include <string>
#include <filesystem>
namespace fs = std::filesystem;

/** \addtogroup Solver_base
  * \ingroup FUKA
  * Various FUKA initial data solvers
  * @{
  */

using namespace ::Kadath::FUKA_Config;
using namespace ::Kadath::FUKA_Config_Utils;
namespace Kadath {
namespace FUKA_Solvers {

/**
 * @brief Get the global path to saved compact object solutions (COs)
 * 
 * @return std::string 
 */
inline std::string get_cos_path() {
  const std::string home_kadath{std::getenv("HOME_KADATH")};
  std::string central_abs{home_kadath+"/COs/"};
  return central_abs;
}

#define FORMAT std::setw(13) << std::left << std::showpos 
const int RELOAD_FILE = 2;
const int RUN_BOOST = 3;

/**
 * @brief Solver base for FUKA solvers
 * 
 * @tparam config_t Configurator type
 * @tparam space_t Numerical space type
 */
template<typename config_t, typename space_t>
class Solver {
  public:
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;
 
  protected:
  space_t& space;
	config_t& bconfig;
  const int ndom;
  STAGES solver_stage;

  public:
  Solver(config_t& config_in, space_t& space_in) 
    : space(space_in), bconfig(config_in), ndom(space_in.get_nbr_domains()) {}
  virtual ~Solver() = default;
  
  protected:
  virtual void print_diagnostics(const System_of_eqs& syst, const int ite, 
    const double conv) const = 0;
  virtual std::string converged_filename(const std::string stage) const = 0;
  virtual void save_to_file() const = 0;
  
  /**
   * @brief Consistent interface for writing a checkpoint
   * 
   * @param termination_chkpt Toggle writing to stdout for termination
   */
  void checkpoint(bool termination_chkpt = false) const  {
    // Backup activated stages
    auto const final_stages{bconfig.return_stages()};

    // Clear stages and only set the current stage as active
    auto & stages = bconfig.return_stages();
    stages.fill(false);
    stages[this->solver_stage] = true;

    // Save to file
    save_to_file();

    // Reset to original stages
    stages = final_stages;
    
    if(termination_chkpt) {
      std::cerr << "***Writing early termination chkpt " 
                << bconfig.config_filename() << "***\n";
      std::_Exit(EXIT_FAILURE);
    }
  }

  /**
   * @brief Branch when maximum iterations are exceeded
   * 
   * @param rank MPI rank
   * @param ite Iteration
   * @param conv Current convergence
   */
  void check_max_iter_exceeded(const int& rank, const int& ite, const double& conv) const {
    bool exceeded = (ite > bconfig.seq_setting(MAX_ITER)) && conv >= bconfig.seq_setting(PREC);
    if(exceeded && \
       conv < 10. * bconfig.seq_setting(PREC)) {
      if(rank == 0)
        std::cout << "Max iterations exceeded at precison, " << conv
                  << "\nFinishing since precision < 10. * PREC....\n"
                  << "Running at higher resolution may help.\n";
    }
    else if(exceeded) {
      if(rank == 0)
        std::cout << "Max iterations exceeded at precison, " << conv;
    }
    else {
      return;
    }
    auto s = converged_filename("termination_chkpt");
    bconfig.set_filename(s);
    if(rank == 0)
      checkpoint(true);
    
    // FIXME this should be handled better than hard termination
    std::_Exit(EXIT_FAILURE);
  }

  /**
   * @brief Branch for handling finding a previous solution
   * 
   * @param last_stage Stage name used when writing to file
   * @return bool
   */
  bool solution_exists(std::string last_stage="") {
    bool exists = false;
    std::string prev_name{converged_filename(last_stage)};
    std::string prev_abs{bconfig.config_outputdir()+prev_name};
    
    const std::string home_kadath{std::getenv("HOME_KADATH")};
    std::string central_abs{home_kadath+"/COs/"+prev_name};
    
    auto check_and_update_config = [&](auto p) {
      exists = true;
      base_config_t old_solution(p+".info");
      if(bconfig.set(BCO_PARAMS::NSHELLS) != old_solution.set(BCO_PARAMS::NSHELLS))
        return false;
      
      // make sure we copy stages, controls, and settings over
      for(auto idx = 0; idx < STAGES::NUM_STAGES; ++idx)
        old_solution.set_stage(idx) = bconfig.set_stage(idx);
      for(auto idx = 0; idx < CONTROLS::NUM_CONTROLS; ++idx)
        old_solution.control(idx) = bconfig.control(idx);
      for(auto idx = 0; idx < SEQ_SETTINGS::NUM_SEQ_SETTINGS; ++idx)
        old_solution.seq_setting(idx) = bconfig.seq_setting(idx);
      
      // Deactivate current stage since we found solution
      old_solution.set_stage(solver_stage) = false;
      bconfig = old_solution;
      return exists;
    };

    /// Check based on current output directory
    if(fs::exists(prev_abs+".info") && fs::exists(prev_abs+".dat")){
      exists = check_and_update_config(prev_abs);
    }
    /// Check based on global output directory
    else if(bconfig.control(SAVE_COS) && fs::exists(central_abs+".info") && fs::exists(central_abs+".dat")) {
      exists = check_and_update_config(central_abs);
    }
    
    return exists;
  }

  /**
   * @brief Extract EOS name from the Config object
   * 
   * @tparam bco_idx optional index type for binary Configs
   * @param bco optional index for binary Configs
   * @return std::string 
   */
  template<typename... bco_idx>
  std::string extract_eos_name(bco_idx... bco) const {
    const std::string eos_file_abs = bconfig.template eos<std::string>(EOSFILE, bco...);
    const std::string eos_file = extract_filename(eos_file_abs);
    return eos_file.substr(0, eos_file.find("."));
  }

  public:
  void set_solver_stage(STAGES const _stage) { this->solver_stage = _stage; }
};

/**
 * @brief Solver base for FUKA solvers
 * 
 * @tparam config_t Configurator type
 * @tparam space_t Numerical space type
 */
template<typename config_t, typename space_t>
class XCTS_Solver : public Solver<config_t, space_t> {
  public:
  using cfgen_t  = CoordFields<space_t>;
  using cfary_t  = std::array<std::optional<Vector>, NUM_VECTORS>;
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;
 
  protected:
  using Solver<config_t, space_t>::space;
  using Solver<config_t, space_t>::bconfig;
  using Solver<config_t, space_t>::solver_stage;

  Base_tensor& basis;
  cfgen_t cfields;
  cfary_t coord_vectors;

  public:
  XCTS_Solver(config_t& config_in, space_t& space_in, Base_tensor& base_in) 
    : Solver<config_t, space_t>(config_in, space_in), cfields(space), basis(base_in) {}
};
/** @}*/
}}