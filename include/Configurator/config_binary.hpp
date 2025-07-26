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
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include "config_bco.hpp"
#include "configurator_boost.hpp"

namespace Kadath {
namespace FUKA_Config {
/**
 * \addtogroup Containers
 * @ingroup Configurator
 * @{*/

/**
 * BIN_INFO
 * Locally store the parameters related to the binary
 * Store pointers to related BCOs which store their own parameters.
 *
 */
class BIN_INFO {

 private:
  const std::string node_t{"binary"};  ///< node type

  /**
   * BIN_INFO::init_bco
   * Static switch to determine the type of BCO to point to.  If one intends
   * to implement new BCOs, this needs to be updated along with MBCO map.
   *
   * @param[input]  idx Index of the bco to be stored
   * @param[input]  bco_type String with the bco type (e.g NS, BH)
   */

  void init_bco(const int idx, const std::string bco_type) {
    switch (MBCO.at(bco_type)) {
      case NODES::BH:
        BCOS[idx] = std::make_unique<BCO_BH_INFO>();
        break;
      case NODES::NS:
        BCOS[idx] = std::make_unique<BCO_NS_INFO>();
        break;
      default:
        std::cerr << "Static switch not defined for " << bco_type << std::endl;
        std::_Exit(EXIT_FAILURE);
        break;
    }
  }

 protected:
  using Array = std::array<double, NUM_BPARAMS>;
  using BArray = std::array<std::unique_ptr<BCO_INFO>, 2>;
  using Map = std::map<std::string, BIN_PARAMS>;

 public:
  BArray BCOS{};             ///< Array of pointers to BCOs
  Array bin_params{};        ///<Array storing Binary parameters
  Map bin_map{MBIN_PARAMS};  ///<Map of parameter strings to enum indexes
  std::map<std::string, STAGES> bin_stages{
      MSTAGE};  ///<Map of only relevant stages

  /**
   * Default constructor initializing pointers to 0
   * and parameters to NaN.  This is a feature, not a
   * bug.
   */
  BIN_INFO() { bin_params.fill(std::nan("1")); }

  /**
   * Copy constructor. Slightly messy due to handling of pointers to
   * BCO parameter containers.
   */
  BIN_INFO(const BIN_INFO& b)
      : bin_params{b.bin_params}, bin_stages{b.bin_stages} {
    for (int i = 0; i < 2; ++i) {
      if (b.BCOS[i]) {
        auto b_type = b.BCOS[i]->get_type();
        if (b_type == "bh") {
          auto child_ptr = dynamic_cast<BCO_BH_INFO*>(b.BCOS[i].get());
          BCOS[i] = std::make_unique<BCO_BH_INFO>(*child_ptr);
        } else if (b_type == "ns") {
          auto child_ptr = dynamic_cast<BCO_NS_INFO*>(b.BCOS[i].get());
          BCOS[i] = std::make_unique<BCO_NS_INFO>(*child_ptr);
        }
      }
    }
  }

  /**
   * Move contructor
   */
  BIN_INFO(BIN_INFO&& b) noexcept
      : bin_params(std::move(b.bin_params)),
        bin_stages(std::move(b.bin_stages)) {
    for (int i = 0; i < 2; ++i)
      BCOS[i] = std::move(b.BCOS[i]);
  }

  /**
   * Assignment operator
   */
  BIN_INFO& operator=(const BIN_INFO& b) {
    if (this == &b)
      return *this;

    BIN_INFO tmp(b);
    *this = std::move(tmp);

    return *this;
  }

  /**
   * Move assignment operator
   */
  BIN_INFO& operator=(BIN_INFO&& b) noexcept {
    this->bin_params = std::move(b.bin_params);
    for (int i = 0; i < 2; ++i)
      this->BCOS[i] = std::move(b.BCOS[i]);
    this->bin_stages = std::move(b.bin_stages);
    return *this;
  }

  /* BIN_INFO::init_binary
   * initialize the pointers of the BCOs of the binary
   *
   * @param[input] bco_types vector of strings containing bco types
   */
  void init_binary(std::array<std::string, 2> bco_types) {
    int idx = 0;
    for (auto& b : bco_types) {
      init_bco(idx, b);
      ++idx;
    }
    initialize_bin_stage_map(bco_types);
  }

  void initialize_bin_stage_map(std::array<std::string, 2> bco_types) {
    if (bco_types[0] == bco_types[1]) {  //check for BBH or BNS
      if (bco_types[0] == "ns")
        bin_stages = MBNSSTAGE;
      if (bco_types[0] == "bh")
        bin_stages = MBBHSTAGE;
    } else if (bco_types[0] == "ns" && bco_types[1] == "bh")
      bin_stages = MBHNSSTAGE;
  }

  /* BIN_INFO::read_params
   * A tree containing all the binary and bco parameters is received.
   * The parameters related to the binary are stored locally.
   * The BCO types are determined based on the branch nodes
   *
   * @param[input] bin_tree Tree containing all parameters of the binary and BCOs
   */
  void read_params(Tree& bin_tree) {
    // Read-in strictly binary related parameters
    read_keys(bin_map, bin_params, bin_tree);

    // Determine BCO types based on branch node names
    std::array<std::string, 2> tmp_ary;
    get_branch_nodes(MBCO, tmp_ary, bin_tree, true);

    // Temporary index used to iterate over BCO pointers
    int tidx = 0;

    // Initialize BCO pointer array based on the BCO types from the Tree
    // and send the branch to BCO_INFO child for reading in the parameters
    for (auto& bco : tmp_ary) {
      init_bco(tidx, bco.substr(0, 2));
      BCOS[tidx]->read_params(bin_tree.get_child(bco.data()));
      bco = bco.substr(0, 2);
      ++tidx;
    }

    initialize_bin_stage_map(tmp_ary);
  }

  /**
   * BIN_INFO::get_type
   * Returns binary type
   *
   * @param[output] node_t node type
   */
  std::string get_type() const { return node_t; }

  /**
   * BIN_INFO::print_me
   * Debug function to print stored node_type
   */
  void print_me() const {
    std::cout << node_t << std::endl;
    std::cout << "-----------------------" << std::endl;
  }

  /* BIN_INFO::return_branch
   * Build a tree containing all the binary and bco parameters.
   * The parameters related to the binary are built locally.
   * The BCO branches are built by BCO_INFO and added to the tree
   *
   * @param[output] branch Tree containing all parameters of the binary and BCOs
   */
  Tree return_branch() {
    Tree branch;
    branch.put_child(node_t, build_branch<Tree>(MBIN_PARAMS, bin_params));
    branch.add_child(node_t + "." + get_name_string(BCO1),
                     BCOS[BCO1]->return_branch());
    branch.add_child(node_t + "." + get_name_string(BCO2),
                     BCOS[BCO2]->return_branch());
    return branch;
  }

  /* BIN_INFO::return_params
   * This function will return the raw arrays, but this is discouraged.
   * It is best to use the overloaded operators to access information
   *
   * @param[output] bin_params Array of binary parameters
   */
  Array& return_params() { return bin_params; }

  /* BIN_INFO::return_bcos
   * This function will return the raw arrays, but this is discouraged.
   * It is best to use the overloaded operators to access information
   *
   * @param[output] BCOs array of BCO pointers
   */
  BArray& return_bcos() { return BCOS; }

  /* BIN_INFO::get_map
   * Returns binary parameter map
   *
   * @param[output] bin_map binary parameter map
   */
  Map const& get_map() const { return bin_map; }

  /* BIN_INFO::get_stage_map
   * Returns binary parameter stages map
   *
   * @param[output] bin_map binary stages map
   */
  const auto& get_stage_map() const { return bin_stages; }

  void set_stage_map(const std::map<std::string, STAGES>& _stages) {
    bin_stages = _stages;
  }

  /* BIN_INFO::get_map
   * Returns BCO parameter map
   *
   * @param[input]  BOCidx BCO pointer index
   * @param[output] bco_map binary parameter map
   */
  auto& get_map(const int BCOidx) const { return BCOS[BCOidx]->get_map(); };

  auto& get_eos_map(const int BCOidx = BCO1) const {
    if (auto child_ptr = dynamic_cast<BCO_NS_INFO*>(BCOS[BCOidx].get())) {
      return child_ptr->get_eos_map();
    }
    throw std::invalid_argument(
        "\nInvalid EOS Parameter indices for assignment\n");
  };

  std::string get_name_string(const int BCOidx) const {
    return BCOS[BCOidx]->get_name_string(BCOidx + 1);
  }

  std::string get_name_string() const { return get_type(); }

  /* BIN_INFO::operator()
   * This operator returns the requested parameter for a given BCO
   *
   * @param[input]  idx Index of the Paramter of interest - see config_enum
   * @param[input]  BCOidx Index of BCO of interest
   * @param[output] *BCOS[BCOidx])(idx) referene to BCO parameter requested
   */
  auto& operator()(const int idx, const int BCOidx) {
    return (*BCOS[BCOidx])(idx);
  }

  /* BIN_INFO::operator()
   * This operator returns the requested parameter for the binary
   *
   * @param[input]  idx Index of the Paramter of interest - see config_enum
   * @param[output] bin_params[idx] referene to binary parameter requested
   */
  auto& operator()(const int idx) { return bin_params[idx]; }

  /* BIN_INFO::set_eos_param
   * Returns a reference to a BCO's EOS parameter to set.  Cannot be used
   * to read the value directly since parameters are stored as std::variant.
   *
   * @param[input]  idx Index of the Paramter of interest - see config_enum
   * @param[input]  BCOidx Index of BCO of interest
   * @param[output] eos_param reference to eos parameter to assign
   * @throws std::invalid_argument Throws when dynamic_cast fails
   */
  auto& set_eos_param(const int idx, const int BCOidx) const {
    if (auto child_ptr = dynamic_cast<BCO_NS_INFO*>(BCOS[BCOidx].get())) {
      return child_ptr->set_eos_param(idx);
    }
    throw std::invalid_argument(
        "\nInvalid EOS Parameter indices for assignment\n");
  }

  /* BIN_INFO::get_eos_param
   * Returns a given EOS parameter
   *
   * @tparam T Parameter indicating the type of the EOS parameter
   * @param[input]  idx Index of the Paramter of interest - see config_enum
   * @param[input]  BCOidx Index of BCO of interest
   * @param[output] eos_param eos parameter value - not assignable.
   * @throws std::invalid_argument Throws when dynamic_cast fails
   */
  template <typename T>
  constexpr T get_eos_param(const int idx, const int BCOidx) const {
    if (auto child_ptr = dynamic_cast<BCO_NS_INFO*>(BCOS[BCOidx].get())) {
      return child_ptr->template get_eos_param<T>(idx);
    }
    throw std::invalid_argument(
        "\nInvalid EOS Parameter indices for reading\n");
  }

  friend std::ostream& operator<<(std::ostream&, const BIN_INFO&);

  /**
   * BIN_INFO::set_defaults
   * Allow the setting of default configurator values for base binary setup
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_defaults(config_t& bconfig) {
    bool includes_matter = false;
    std::array<double, 2> Ms{};
    // copy Compact Object defaults
    for (auto& bco : {NODES::BCO1, NODES::BCO2}) {
      const auto bcotype = BCOS[bco]->get_type();
      if (bcotype == "ns") {
        kadath_config_boost<BCO_NS_INFO> nsconfig;
        nsconfig.set_defaults();
        for (auto i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
          bconfig.set(i, bco) = nsconfig.set(i);
        for (auto i = 0; i < EOS_PARAMS::NUM_EOS_PARAMS; ++i)
          bconfig.set_eos(i, bco) = nsconfig.set_eos(i);
        includes_matter = true;
        Ms[bco] = bconfig(BCO_PARAMS::MADM, bco);
      } else if (bcotype == "bh") {
        kadath_config_boost<BCO_BH_INFO> bhconfig;
        bhconfig.set_defaults();
        for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
          bconfig.set(i, bco) = bhconfig.set(i);
        Ms[bco] = bconfig(BCO_PARAMS::MCH, bco);
      }
    }
    // Set binary defaults
    bconfig.set(BIN_PARAMS::BIN_RES) = 9.;
    bconfig.set(BIN_PARAMS::DIST) = 10. * Ms[0] + Ms[1];
    bconfig.set(BIN_PARAMS::REXT) = 2. * bconfig(BIN_PARAMS::DIST);
    bconfig.set(BIN_PARAMS::QPIG) = bconfig(BCO_PARAMS::BCO_QPIG, NODES::BCO1);
    bconfig.set(BIN_PARAMS::Q) = Ms[1] / Ms[0];
    bconfig.set(BIN_PARAMS::COM) = 0.;
    bconfig.set(BIN_PARAMS::COMY) = 0.;
    bconfig.set(BIN_PARAMS::OUTER_SHELLS) = 0;
    bconfig.set(BIN_PARAMS::GOMEGA) = 0.;

    // Set default fields
    bconfig.set_field(BCO_FIELDS::SHIFT) = true;
    bconfig.set_field(BCO_FIELDS::LAPSE) = true;
    bconfig.set_field(BCO_FIELDS::CONF) = true;
    if (includes_matter) {
      bconfig.set_field(BCO_FIELDS::LOGH) = true;
      bconfig.set_field(BCO_FIELDS::PHI) = true;
    }

    // The following controls are required for an
    // initial solution using the v2 Solvers
    bconfig.control(CONTROLS::SEQUENCES) = true;
    bconfig.control(CONTROLS::USE_BOOSTED_CO) = true;
    bconfig.control(CONTROLS::FIXED_GOMEGA) = true;
    bconfig.control(CONTROLS::NEW_ID) = true;

    bconfig.set_stage(STAGES::TOTAL_BC) = true;
    bconfig.set_stage(STAGES::ECC_RED) = true;

    bconfig.seq_setting(SEQ_SETTINGS::INIT_RES) = 9;
  }

  /**
   * BIN_INFO::set_minimal_defaults
   * Allow the setting of default configurator values for base binary setup
   *
   * @tparam config_t configuration file type
   * @param bconfig reference to configuration file to be modified
   */
  template <typename config_t>
  void set_minimal_defaults(config_t& bconfig) {
    bool includes_matter = false;
    std::array<double, 2> Ms{};
    // copy Compact Object defaults
    for (auto& bco : {NODES::BCO1, NODES::BCO2}) {
      const auto bcotype = BCOS[bco]->get_type();
      if (bcotype == "ns") {
        kadath_config_boost<BCO_NS_INFO> nsconfig;
        nsconfig.set_minimal_defaults();
        for (auto i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
          bconfig.set(i, bco) = nsconfig.set(i);
        for (auto i = 0; i < EOS_PARAMS::NUM_EOS_PARAMS; ++i)
          bconfig.set_eos(i, bco) = nsconfig.set_eos(i);
        includes_matter = true;
        Ms[bco] = bconfig(BCO_PARAMS::MADM, bco);
      } else if (bcotype == "bh") {
        kadath_config_boost<BCO_BH_INFO> bhconfig;
        bhconfig.set_minimal_defaults();
        for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
          bconfig.set(i, bco) = bhconfig.set(i);
        Ms[bco] = bconfig(BCO_PARAMS::MCH, bco);
      }
    }
    // Set binary defaults
    bconfig.set(BIN_PARAMS::BIN_RES) = 9.;
    bconfig.set(BIN_PARAMS::DIST) = 10. * Ms[0] + Ms[1];

    bconfig.set(BIN_PARAMS::OUTER_SHELLS) = 0;

    bconfig.set_stage(STAGES::TOTAL_BC) = true;
    bconfig.set_stage(STAGES::ECC_RED) = true;

    bconfig.control(CONTROLS::SAVE_COS) = false;
    bconfig.control(CONTROLS::NEW_ID) = true;
  }
};

/**
 * @}*/

std::ostream& operator<<(std::ostream& out, const BIN_INFO& BIN);
}  // namespace FUKA_Config
}  // namespace Kadath
