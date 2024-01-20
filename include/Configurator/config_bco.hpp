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
#include "configurator_boost.hpp"
#include <cmath>
#include <variant>
namespace Kadath {
namespace FUKA_Config {
/**
 * \addtogroup Containers
 * @ingroup Configurator
 * @{*/

/**
 * BCO_INFO
 * Each BCO retains it's own configuration data.
 * This includes the tree that it can update which is later written to file.
 *
 */
class BCO_INFO {

private:
  std::string node_t{"bco"}; ///< string containing node type

protected:
  using Array = std::array<double, NUM_BCO_PARAMS>;
  using Map   = std::map<std::string, BCO_PARAMS>;

  /**
   * std::map for mapping configuration strings to enum indexes
   */
  const Map bco_map{MBCO_PARAMS};
  std::map<std::string, STAGES> bco_stages{MSTAGE};
  
  /** 
   *Array to store parameters
   */
  Array bco_params{};

public:
  /**
   * Constructor.  Everything is set to nan.
   * this is a feature - not a bug
   */
  BCO_INFO() : bco_stages(MSTAGE) {
    bco_params.fill(std::nan("1"));
  }

  /**
   * BCO_INFO copy constructor
   */
  BCO_INFO(const BCO_INFO& b) : bco_params(b.bco_params), bco_stages(b.bco_stages) {}

  /**
   * BCO_INFO move constructor
   */
  BCO_INFO(BCO_INFO&& b) noexcept 
    : bco_params(std::move(b.bco_params)), bco_stages(std::move(b.bco_stages)) {}

  /**
   * BCO_INFO assignment operator
   */
  BCO_INFO& operator=(const BCO_INFO& b) {
    this->bco_params = b.bco_params;
    this->bco_stages = b.bco_stages;
    return *this;
  }

   /**
   * BCO_INFO move assignment operator
   */
  BCO_INFO& operator=(BCO_INFO&& b) noexcept {
    this->bco_params = std::move(b.bco_params);
    this->bco_stages = std::move(b.bco_stages);
    return *this;
  }

  /**
   * BCO_INFO::return_branch
   * Build a branch based on non-nan parameters
   *
   * param[output] branch branch containing BCO paramters
   */
  virtual Tree return_branch() {
    return build_branch<Tree>(MBCO_PARAMS, bco_params);
  }
  
  /**
   * BCO_INFO::get_name_string
   * Returns BCO type with an added int, 1 or 2.  Needed
   * for config_binary.hpp 
   *
   * @param[input]  N integer to concatenate with type
   * @param[output] string concatenated string
   */
  virtual std::string get_name_string(const int N = 0) {
    return (N == 0) ? node_t.data() : node_t.data() + std::to_string(N);
  }

  /**
   * BCO_INFO::print_me
   * Debug function to print stored node_type
   */
  virtual void print_me() {
    std::cout << node_t << std::endl;
    std::cout << "-----------------------" << std::endl;
  }

  /**
   * BCO_INFO::read_params
   * Store parameters from an input branch
   *
   * param[input] branch branch containing BCO paramters
   */
  virtual void read_params(Tree &branch) {
    read_keys(bco_map, bco_params, branch);
  }

  /**
   * BCO_INFO::get_type
   * Returns BCO type
   *
   * @param[output] node_t node type
   */
  virtual std::string get_type() const { return node_t; }

  /**
   * BCO_INFO::get_map
   * Returns BCO map which contains the mapping from configuration
   * strings to enumerator indexes.  See config_enum.hpp
   *
   * @param[output] bco_map Parameter map
   */
  Map const &get_map() const { return bco_map; }
  
  /* BIN_INFO::get_stage_map
   * Returns binary parameter stages map
   *
   * @param[output] bin_map binary parameter map
   */
  virtual std::map<std::string, STAGES> const & get_stage_map() const { 
    return bco_stages; 
  }

  /**
   * BCO_INFO::operator()
   * Returns a given BCO parameter based on a given parameter
   * index contained in bco_map
   *
   * @param[output] bco_params[idx] requested BCO parameter
   */
  auto &operator()(const int idx) { return bco_params[idx]; }

  friend std::ostream &operator<<(std::ostream &, const BCO_INFO &);
};

/**
 * BCO_BH_INFO
 * Child class containing parameters for a basic BH.
 * We also store default parameters which is used for an initial setup.
 */

class BCO_BH_INFO : public BCO_INFO {
private:
  std::string node_t = "bh"; ///< string containing node type


public:
  BCO_BH_INFO() : BCO_INFO() {
    bco_stages = MBHSTAGE;
  }
  /**
   * BCO_INFO copy constructor
   */
  BCO_BH_INFO(const BCO_BH_INFO& b) : BCO_INFO(b) {}

  /**
   * BCO_INFO move constructor
   */
  BCO_BH_INFO(BCO_BH_INFO&& b) noexcept : BCO_INFO(std::move(b)) {}

  /**
   * BCO_INFO assignment operator
   */
  BCO_BH_INFO& operator=(BCO_BH_INFO& b) {
    this->bco_params = b.bco_params;
    this->bco_stages = b.bco_stages;
    return *this;
  }

   /**
   * BCO_INFO move assignment operator
   */
  BCO_BH_INFO& operator=(const BCO_BH_INFO&& b) noexcept {
    this->bco_params = std::move(b.bco_params);
    this->bco_stages = std::move(b.bco_stages);
    return *this;
  }

  /**
   * BCO_BH_INFO::get_name_string
   * Returns BCO type with an added int, 1 or 2.  Needed
   * for config_binary.hpp 
   *
   * @param[input]  N integer to concatenate with type
   * @param[output] string concatenated string
   */
  virtual std::string get_name_string(const int N = 0) override {
    return (N == 0) ? node_t.data() : node_t.data() + std::to_string(N);
  }

  /**
   * BCO_BH_INFO::print_me
   * Debug function to print stored node_type
   */
  virtual void print_me() override {
    std::cout << node_t << std::endl;
    std::cout << "-----------------------" << std::endl;
  }

  /**
   * BCO_BH_INFO::get_type
   * Returns BCO type
   *
   * @param[output] node_t node type
   */
  virtual std::string get_type() const override { return node_t; }

  /** 
   * BCO_BH_INFO::set_defaults
   * Allow the setting of default configurator values for a base BH setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_defaults(config_t& bconfig) {
    // start - set BH properties in config file
    // Resolution of the initial setup
    bconfig.set(BCO_PARAMS::BCO_RES) = 9;
    bconfig.set(BCO_PARAMS::DIM)     = 3;
    
    // Units of the system - 4 * PI * G
    bconfig.set(BCO_PARAMS::BCO_QPIG) = 4 * M_PI;
    
    // Nucleus Radius 
    bconfig.set(BCO_PARAMS::RIN)  = 0.1;
    
    // Initial guess for BH - radius
    bconfig.set(BCO_PARAMS::RMID) = 0.3;
    
    // Radius of outer boundary - matched with compactified domain
    bconfig.set(BCO_PARAMS::ROUT) = 5 * bconfig(RMID) ;
    
    // Additional shells can be added between RMID and ROUT for refined resolution
    bconfig.set(BCO_PARAMS::NSHELLS) = 0;
    
    // Initial dimensionless spin
    bconfig.set(BCO_PARAMS::CHI) = 0;
    
    // Initial guess of omega
    bconfig.set(BCO_PARAMS::OMEGA) = 0;

    // MIRR and MCH are fixing parameters...this is what the resulting BH will be
    bconfig.set(BCO_PARAMS::MIRR)  = 0.5 ;
    bconfig.set(BCO_PARAMS::MCH)   = 0.5 ;

    bconfig.set(BCO_PARAMS::BVELX) = 0.;
    bconfig.set(BCO_PARAMS::BVELY) = 0.;
    
    /** 
     * fixed lapse is used for the PRE stage only.  
     * Once the system is solved using the Neumann lapse condition
     * this value is updated in the config file in the standard solver */
    bconfig.set(BCO_PARAMS::FIXED_LAPSE)  = .3;
    //end   - set BH parameters

    //start - set BH stages in config file
    bconfig.set_stage(STAGES::TOTAL_BC) = true;
    //end   - set BH stages

    //start - set BH fields in config file
    bconfig.set_field(BCO_FIELDS::CONF)  = true;
    bconfig.set_field(BCO_FIELDS::LAPSE) = true;
    bconfig.set_field(BCO_FIELDS::SHIFT) = true;
    //end   - set BH fields

    bconfig.seq_setting(SEQ_SETTINGS::INIT_RES) = 9;
    bconfig.control(CONTROLS::SAVE_COS) = false;
  }

    /** 
   * BCO_BH_INFO::set_minimal_defaults
   * Allow the setting of default configurator values for a base BH setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_minimal_defaults(config_t& bconfig) {
    //start - set BH properties in config file
    //Resolution of the initial setup
    bconfig.set(BCO_PARAMS::BCO_RES) = 9;
    
    //Initial dimensionless spin
    bconfig.set(BCO_PARAMS::CHI) = 0;
    bconfig.set(BCO_PARAMS::MCH) = 0.5 ;
    //end   - set BH parameters
    bconfig.set_stage(STAGES::TOTAL_BC) = true;
  }
};

/**
 * BCO_NS_INFO
 * Child class containing parameters for a NS
 * We also store default parameters which is used for an initial setup.
 *
 * vars_t is a std::variant container that allows the storage of various
 * types while being type safe.  This is extremely important when
 * considering the EOS could be simply a file name for a table, or this
 * can be used to store information relation to a piecewise polytrope EOS.
 * This leaves a lot of room for modification later while only needing to
 * modify config_enums for the added parameters.
 *
 * EOSMap is an additional map that only handles EOS parameters.
 */

class BCO_NS_INFO : public BCO_INFO {
public:  
  using vars_t  = std::variant<double, int, std::string>;
  using EOSArray = std::array<vars_t, NUM_BCO_PARAMS>;
  using EOSMap   = std::map<std::string, EOS_PARAMS>;
  using DIFFROT_ary = std::array<vars_t, DIFFROT_PARAMS::NUM_DIFFROT_PARAMS>;
  using DIFFROT_map = std::map<std::string, DIFFROT_PARAMS>;

protected:
  std::string node_t{"ns"}; ///< node type

  EOSArray eos_params{}; ///< Array storing EOS parameters

  //std::map for mapping configuration strings to EOS enum indexes
  const EOSMap eos_map{MEOS_PARAMS};

  DIFFROT_ary diffrot_params{};
  DIFFROT_map const diffrot_map{MDIFFROT_PARAMS};

public:
  /** 
   * BCO_NS_INFO::BCO_NS_INFO
   * In addition to the parent constructor, we initialize the
   * the EOS array to NaN as well. 
   */
  BCO_NS_INFO() : BCO_INFO() {
    bco_stages = MNSSTAGE;
    eos_params.fill(std::nan("1"));
    diffrot_params.fill(std::nan("1"));
  }
  
  /**
   * BCO_NS_INFO copy constructor
   */
  BCO_NS_INFO(const BCO_NS_INFO& b) 
    : BCO_INFO(b), 
      eos_params{b.eos_params}, diffrot_params(b.diffrot_params) 
      { }

  /**
   * BCO_NS_INFO move constructor
   */
  BCO_NS_INFO(BCO_NS_INFO&& b) noexcept 
    : BCO_INFO(std::move(b)), 
      eos_params(std::move(b.eos_params)),
      diffrot_params(std::move(b.diffrot_params))
      { }
  
  /**
   * BCO_NS_INFO assignment operator
   */
  BCO_NS_INFO& operator=(const BCO_NS_INFO& b) {
    this->bco_params = b.bco_params;
    this->bco_stages = b.bco_stages;
    this->eos_params = b.eos_params;
    this->diffrot_params = b.diffrot_params;
    return *this;
  }
  
  /**
   * BCO_NS_INFO move assignment operator
   */
  BCO_NS_INFO& operator=(BCO_NS_INFO&& b) noexcept {
    this->bco_params = std::move(b.bco_params);
    this->bco_stages = std::move(b.bco_stages);
    this->eos_params = std::move(b.eos_params);
    this->diffrot_params = std::move(b.diffrot_params);
    return *this;
  }

  /**
   * BCO_NS_INFO::return_eos_params
   * Return a reference to the EOS array.  Mainly for testing.
   * Recommended to use get/set_eos_param for safety
   * param[out] eos_params
   */
  const EOSArray& return_eos_params() const { return eos_params; }
  /**
   * BCO_NS_INFO::return_eos_map
   * Return a reference to the EOS parameter map
   * param[out] eos_params
   */
  EOSMap const& get_eos_map() const { return eos_map; }

  /**
   * BCO_NS_INFO::return_diffrot_params
   * Return a reference to the differential rotation array.  Mainly for testing.
   * Recommended to use get_diffort_param for safety
   * param[out] diffrot_params
   */
  const DIFFROT_ary& return_diffrot_params() const { return diffrot_params; }
  /**
   * BCO_NS_INFO::return_diffrot_map
   * Return a reference to the differential rotation parameter map
   * param[out] map_diffrot_params
   */
  DIFFROT_map const& get_diffrot_map() const { return diffrot_map; }

  /**
   * BCO_NS_INFO::get_name_string
   * Returns BCO type with an added int, 1 or 2.  Needed
   * for config_binary.hpp 
   *
   * @param[input]  N integer to concatenate with type
   * @param[output] string concatenated string
   */
  virtual std::string get_name_string(const int N = 0) override {
    return (N == 0) ? node_t.data() : node_t.data() + std::to_string(N);
  }

  /**
   * BCO_NS_INFO::print_me
   * Debug function to print stored node_type
   */
  virtual void print_me() override {
    std::cout << node_t << std::endl;
    std::cout << "-----------------------" << std::endl;
  }

  /**
   * BCO_NS_INFO::read_params
   * Reads the mapped parameters from a given branch
   * 
   * @param[input] branch input branch to read bco and eos params from
   */
  virtual void read_params(Tree &branch) override {
    read_keys(bco_map, bco_params, branch);
    read_keys(eos_map, eos_params, branch);
    
    if(branch.find("differential_rotation") != branch.not_found())
      read_keys(diffrot_map, diffrot_params, read_branch(branch, "differential_rotation"));
    // if(!is_storage_all_nan(diffrot_params)) {
    //   print_params(diffrot_map, diffrot_params);
    // }
  }

  /**
   * BCO_NS_INFO::return_branch
   * Builds a single branch from the stored bco and eos parameters
   * 
   * @param[output] branch branch containing bco and eos params
   */
  //
  virtual Tree return_branch() override {
    Tree branch(build_branch<Tree>(bco_map, bco_params));
    Tree eos_childs = build_branch<Tree>(eos_map, eos_params);
    for (auto child : eos_childs)
      branch.push_back(child);
    
    if(!is_storage_all_nan(diffrot_params)) {
      Tree diffrot_branch(build_branch<Tree>(diffrot_map, diffrot_params));
      branch.add_child("differential_rotation", diffrot_branch);
    }    
    return branch;
  }

  /**
   * BCO_NS_INFO::get_type
   * Returns BCO type
   *
   * @param[output] node_t node type
   */
  virtual std::string get_type() const override { return node_t; }
  
  /**
   * BCO_NS_INFO::set_eos_param
   * Returns reference to EOS param allowing value assignment
   * Note type is vars_t which is a std::variant
   *
   * @param[output] eos_param[idx]
   */
  auto& set_eos_param(const int idx) { return eos_params[idx]; }

  /**
   * BCO_NS_INFO::get_eos_param
   * For accessing values from a std::variant, we need to specify
   * the expected type to access it. Returns are value only.
   *
   * @tparam T type of EOS parameter - See vars_t for options
   * @param[output] eos_param[idx]
   */
  template<typename T>
  const T get_eos_param(const int idx) const { return std::get<T>(eos_params[idx]); }

/**
   * BCO_NS_INFO::set_diffrot_param
   * Returns reference to EOS param allowing value assignment
   * Note type is vars_t which is a std::variant
   *
   * @param[output] eos_param[idx]
   */
  auto& set_diffrot_param(const int idx) { return diffrot_params[idx]; }

  /**
   * BCO_NS_INFO::get_diffrot_param
   * For accessing values from a std::variant, we need to specify
   * the expected type to access it. Returns are value only.
   *
   * @tparam T type of EOS parameter - See vars_t for options
   * @param[output] eos_param[idx]
   */
  template<typename T>
  const T get_diffrot_param(const int idx) const { return std::get<T>(diffrot_params[idx]); }

  /** 
   * BCO_NS_INFO::set_defaults
   * Allow the setting of default configurator values for a base NS setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_defaults(config_t& bconfig) {
    // start - set NS properties in config file
    bconfig.set_eos(EOS_PARAMS::EOSFILE)    = "togashi.lorene";
    bconfig.set_eos(EOS_PARAMS::EOSTYPE)    = "Cold_Table";
    bconfig.set_eos(EOS_PARAMS::HCUT)       = 0.0;
    bconfig.set_eos(EOS_PARAMS::INTERP_PTS) = 2000;
    
    bconfig.set(BCO_PARAMS::HC) = 1.26;
    bconfig.set(BCO_PARAMS::NC) = 1.37e-3;

    // Resolution of the initial setup
    bconfig.set(BCO_PARAMS::BCO_RES) = 9;
    bconfig.set(BCO_PARAMS::DIM)     = 3;
    
    // Units of the system - 4 * PI * G
    bconfig.set(BCO_PARAMS::BCO_QPIG) = 4 * M_PI;
    
    // Initial guess for NS - radius
    bconfig.set(BCO_PARAMS::RMID)     = 6.2;
    
    // Nucleus Radius 
    bconfig.set(BCO_PARAMS::RIN)      = 0.5 * bconfig(RMID);

    // Radius of outer boundary - matched with compactified domain
    bconfig.set(BCO_PARAMS::ROUT)     = 1.5 * bconfig(RMID) ;
    
    // Additional shells between RIN and RMID
    bconfig.set(BCO_PARAMS::NINSHELLS) = 0;
    // Additional shells between RMID and ROUT
    bconfig.set(BCO_PARAMS::NSHELLS)   = 0;
    
    // Initial dimensionless spin
    bconfig.set(BCO_PARAMS::CHI)   = 0;
    
    // Initial guess of omega
    bconfig.set(BCO_PARAMS::OMEGA) = 0;

    // We initialize based on fixed MADM
    bconfig.set(BCO_PARAMS::MADM)    = 1.4 ;
    bconfig.set(BCO_PARAMS::QLMADM)  = 1.4 ;
    bconfig.set(BCO_PARAMS::MB)      = 1.55 ;
    // end   - set NS parameters

    // start - set NS stages in config file
    bconfig.set_stage(STAGES::PRE)      = false;
    bconfig.set_stage(STAGES::NOROT_BC) = true;
    bconfig.set_stage(STAGES::TOTAL_BC) = true;
    // end   - set NS stages

    // start - set NS fields in config file
    bconfig.set_field(BCO_FIELDS::SHIFT)    = true;
    bconfig.set_field(BCO_FIELDS::LAPSE)    = true;
    bconfig.set_field(BCO_FIELDS::CONF)     = true;
    bconfig.set_field(BCO_FIELDS::LOGH)     = true;
    // end   - set NS fields

    bconfig.seq_setting(SEQ_SETTINGS::INIT_RES) = 9;
    bconfig.control(CONTROLS::SAVE_COS) = false;
  }

  /** 
   * BCO_NS_INFO::set_minimal_defaults
   * Allow the setting of default configurator values for a base NS setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_minimal_defaults(config_t& bconfig) {
    // start - set NS properties in config file
    bconfig.set_eos(EOS_PARAMS::EOSFILE)    = "togashi.lorene";
    bconfig.set_eos(EOS_PARAMS::EOSTYPE)    = "Cold_Table";

    // Resolution of the initial setup
    bconfig.set(BCO_PARAMS::BCO_RES) = 9;

    // Initial dimensionless spin
    bconfig.set(BCO_PARAMS::CHI)     = 0;

    // We initialize based on fixed MADM
    bconfig.set(BCO_PARAMS::MADM)    = 1.4 ;
    // end   - set NS parameters

    // start - set NS stages in config file
    bconfig.set_stage(STAGES::NOROT_BC) = true;
    bconfig.set_stage(STAGES::TOTAL_BC) = true;
    // end   - set NS stages
  }
};

class BCO_ISO_NS_INFO : public BCO_NS_INFO {
  public:
  /** 
   * BCO_NS_INFO::BCO_NS_INFO
   * In addition to the parent constructor, we initialize the
   * the EOS array to NaN as well. 
   */
  BCO_ISO_NS_INFO() : BCO_NS_INFO() {
    bco_stages = M2DNSSTAGE;
  }

  /** 
   * BCO_NS_INFO::set_defaults
   * Allow the setting of default configurator values for a base NS setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_defaults(config_t& bconfig) {
    // start - set NS properties in config file
    bconfig.set_eos(EOS_PARAMS::EOSFILE)    = "togashi.lorene";
    bconfig.set_eos(EOS_PARAMS::EOSTYPE)    = "Cold_Table";
    bconfig.set_eos(EOS_PARAMS::HCUT)       = 0.0;
    bconfig.set_eos(EOS_PARAMS::INTERP_PTS) = 2000;
    
    bconfig.set(BCO_PARAMS::HC) = 1.26;
    bconfig.set(BCO_PARAMS::NC) = 1.37e-3;

    // Resolution of the initial setup
    bconfig.set(BCO_PARAMS::BCO_RES) = 9;
    bconfig.set(BCO_PARAMS::DIM)     = 2;
    
    // Units of the system - 4 * PI * G
    bconfig.set(BCO_PARAMS::BCO_QPIG) = 4 * M_PI;

    // Additional shells between RMID and ROUT
    bconfig.set(BCO_PARAMS::NSHELLS)   = 0;
    
    // Initial dimensionless spin
    bconfig.set(BCO_PARAMS::CHI)   = 0;
    
    // Initial guess of omega
    bconfig.set(BCO_PARAMS::OMEGA) = 0;

    // We initialize based on fixed MADM
    bconfig.set(BCO_PARAMS::MADM)    = 1.4 ;
    bconfig.set(BCO_PARAMS::MB)      = 1.5 ;
    // end   - set NS parameters

    // start - set NS stages in config file
    bconfig.set_stage(STAGES::NOROT_BC) = true;
    bconfig.set_stage(STAGES::UNIFORM_ROT) = true;
    // end   - set NS stages

    // start - set NS fields in config file
    
    // Make sure all other fields are disabled
    // Starting from scratch, we only need three fields for the
    // non-rotating case.
    for(auto i = 0; i < BCO_FIELDS::NUM_BCO_FIELDS; i++)
      bconfig.set_field(i) = false;
    
    // Document fields that will be stored
    bconfig.set_field(BCO_FIELDS::LOGH)   = true;
    bconfig.set_field(BCO_FIELDS::NU)     = true;
    bconfig.set_field(BCO_FIELDS::LAP_ATERM) = true;
    bconfig.set_field(BCO_FIELDS::LAP_BTERM) = true;
    // end   - set NS fields

    bconfig.seq_setting(SEQ_SETTINGS::INIT_RES) = 9;
    bconfig.seq_setting(SEQ_SETTINGS::MAX_ITER) = 20;
    bconfig.control(CONTROLS::SAVE_COS) = false;
  }

  /** 
   * BCO_NS_INFO::set_minimal_defaults
   * Allow the setting of default configurator values for a base NS setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_minimal_defaults(config_t& bconfig) {
    // start - set NS properties in config file
    bconfig.set_eos(EOS_PARAMS::EOSFILE)    = "gam2.polytrope";
    bconfig.set_eos(EOS_PARAMS::EOSTYPE)    = "Cold_PWPoly";

    // Resolution of the initial setup
    bconfig.set(BCO_PARAMS::BCO_RES) = 9;

    // Initial omega
    bconfig.set(BCO_PARAMS::OMEGA)   = 0;

    // Currently only HC is reliable...
    bconfig.set(BCO_PARAMS::MADM)    = 1.4 ;
    bconfig.set(BCO_PARAMS::HC) = 1.26;
    // end   - set NS parameters

    // start - set NS stages in config file
    bconfig.set_stage(STAGES::NOROT_BC) = true;
    bconfig.set_stage(STAGES::UNIFORM_ROT) = true;
    // end   - set NS stages

    bconfig.seq_setting(SEQ_SETTINGS::MAX_ITER) = 20;
  }
};

/**
 * BCO_KSBH_INFO
 * Child class containing parameters for a KerrSchild BH.
 * We also store default parameters which is used for an initial setup.
 */

class BCO_KSBH_INFO : public BCO_INFO {
private:
  std::string node_t = "bh"; ///< string containing node type

public:
  BCO_KSBH_INFO() : BCO_INFO() {
    bco_stages = MKSBHSTAGE;
  }
  BCO_KSBH_INFO(const BCO_KSBH_INFO& b) : BCO_INFO(b) {}
  BCO_KSBH_INFO(BCO_KSBH_INFO&& b) noexcept : BCO_INFO(std::move(b)) {}
  BCO_KSBH_INFO& operator=(BCO_KSBH_INFO& b) {
    this->bco_params = b.bco_params;
    this->bco_stages = b.bco_stages;
    return *this;
  }
  BCO_KSBH_INFO& operator=(const BCO_KSBH_INFO&& b) noexcept {
    this->bco_params = std::move(b.bco_params);
    this->bco_stages = std::move(b.bco_stages);
    return *this;
  }

  /**
   * BCO_KSBH_INFO::get_name_string
   * Returns BCO type with an added int, 1 or 2.  Needed
   * for config_binary.hpp 
   *
   * @param[input]  N integer to concatenate with type
   * @param[output] string concatenated string
   */
  virtual std::string get_name_string(const int N = 0) override {
    return (N == 0) ? node_t.data() : node_t.data() + std::to_string(N);
  }

  /**
   * BCO_KSBH_INFO::print_me
   * Debug function to print stored node_type
   */
  virtual void print_me() override {
    std::cout << node_t << std::endl;
    std::cout << "-----------------------" << std::endl;
  }

  /**
   * BCO_KSBH_INFO::get_type
   * Returns BCO type
   *
   * @param[output] node_t node type
   */
  virtual std::string get_type() const override { return node_t; }

  /** 
   * BCO_KSBH_INFO::set_defaults
   * Allow the setting of default configurator values for a base BH setup - do not modify 
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_defaults(config_t& bconfig) {
    //start - set BH properties in config file
    //Resolution of the initial setup
    bconfig.set(BCO_RES)  = 9;
    bconfig.set(DIM)      = 3;
    
    //Units of the system - 4 * PI * G
    bconfig.set(BCO_QPIG) = 4 * M_PI;
    
    //Nucleus Radius 
    bconfig.set(RIN)      = 1.;
    
    //Initial guess for BH - radius
    bconfig.set(RMID)     = 2.;
    
    //Radius of outer boundary - matched with compactified domain
    bconfig.set(ROUT)     = 2. * bconfig(RMID) ;
    
    //Additional shells can be added between RMID and ROUT for refined resolution
    bconfig.set(NSHELLS)  = 1;
    
    //Initial dimensionless spin
    bconfig.set(CHI)      = 0;
    
    //Initial guess of omega
    bconfig.set(OMEGA)    = 0;

    //MIRR and MCH are fixing parameters...this is what the resulting BH will be
    bconfig.set(MIRR)     = 1. ;
    bconfig.set(MCH)      = 1. ;
    bconfig.set(KERR_MCH) = 1. ;
    bconfig.set(KERR_CHI) = 0. ;

    bconfig.set(BVELX)    = 0.;
    bconfig.set(BVELY)    = 0.;
    //end   - set BH parameters

    //start - set BH stages in config file
    bconfig.set_stage(TOTAL_BC) = true;
    //end   - set BH stages

    //start - set BH fields in config file
    bconfig.set_field(CONF)  = true;
    bconfig.set_field(LAPSE) = true;
    bconfig.set_field(SHIFT) = true;
    bconfig.set_field(KS_K)  = true;
    bconfig.set_field(KS_LAPSE)  = true;
    bconfig.set_field(KS_METRIC)  = true;
    //end   - set BH fields

    bconfig.seq_setting(INIT_RES) = 9;
    bconfig.control(SAVE_COS) = true;
  }
};
/**
 * @}*/

std::ostream &operator<<(std::ostream &out, const BCO_INFO &BCO);
}}