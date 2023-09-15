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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/optional.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include "config_enums.hpp"
#include <stdexcept>
#include <variant>
#include <iostream>
#include <iomanip>
#include <string>
#include <typeinfo>
#include <type_traits>
#include <sstream>
#include <cmath>
namespace pt = boost::property_tree;

/** @defgroup Configurator_Utils
  * @ingroup Configurator
	* Utility functions to enable parsing and saving Configurator files currently
	* based on the BOOST proprietary INFO format.  In principle this can be changed easily
	* by changing the include to another BOOST parser, however, this is untested.
  * @{
  */
namespace Kadath {
namespace FUKA_Config_Utils {
/**
  * check_for_nan
  *
  * Checks to determine if a variable is NaN.  Execution fails if variable is NaN
	* or not found in the map at all.
  * 
  * @tparam map_t map type
	* @tparam T var type
	* @param[input] map: std::map containing enum/string maps - used only for error deduction
  * @param[input] var: variable to determine if it's NaN or not found
  * @param[input] idx: index of parameter - used only for error deduction
  * @return bool
  */
template<typename map_t, typename T>
bool check_for_nan(const map_t& map, const T& var, const int idx);

/**
  * read_branch
  *
  * Given an input tree, return a branch based on a give node name
  * 
  * @tparam tree_t tree type
  * @param[input] tree: tree to read branch from
  * @param[input] node: node string branch should come from
  * @return branch
  */
template<typename tree_t>
tree_t read_branch(const tree_t& tree, std::string node);

/**
  * get_branch_nodes
  *
	* store node names into storage array based on expected node names map.
	* rsuffix = remove suffix - necessary to remove 1 or 2 from bh1, etc since we map only the pure names (bh, ns, etc)
  * 
  * @tparam ary_t array type
  * @tparam map_t map type
  * @tparam tree_t tree type
  * @param[input]  storage_map: map of keys that we know how to store
  * @param[output] storage array: to store read-in keys into
  * @param[input]  tree: tree to read branch from
  * @param[input]  rsuffix: whether we need to remove suffix (e.g. NS1 -> NS) 
  */
template<typename ary_t, typename map_t, typename tree_t>
void get_branch_nodes(const map_t& storage_map, ary_t& storage, const tree_t& branch, bool rsuffix = false);

/**
  * read_keys
  *
  * from a tree/branch, read all keys and store those found in storage_map into storage
	* -Unknown parameters are ignored.
  * -Boolean parameters not found are set to false
  * 
  * @tparam ary_t array type
  * @tparam map_t map type
  * @tparam tree_t tree type
  * @param[input]  storage_map: map of keys that we know how to store
  * @param[output] storage array: to store read-in keys into
  * @param[input]  tree: tree to read branch from
  */
template<typename ary_t, typename map_t, typename tree_t>
void read_keys(const map_t& storage_map, ary_t& storage, const tree_t& tree) noexcept;

/**
  * print_params
  *
	* useful print tool that shows only the saved parameters in a stored Array given the proper Map
	* i.e. parameter != FALSE && != NaN
	*
	* @tparam ary_t type of storage array
	* @tparam map_t type of string->index map
	* @param[input] storage_map: std::map containing enum/string of known parameters
  * @param[input] storage: parameter storage array
  * @param[output] storage array: to store read-in keys into
  */
template <typename map_t, typename ary_t>
void print_params(const map_t& storage_map, const ary_t& storage, std::ostream& out=std::cout);

/**
  * build_branch
  *
  * from a storage array, prepare a branch for being written to file
  * 
  * @tparam ary_t array type
  * @tparam map_t map type
  * @tparam tree_t branch type
  * @param[input] storage_map: map of keys that we know how to store
  * @param[input] storage array: to store read-in keys into
  * @param[input] inc_off: control whether to print disable parameters
  * @return       branch branch built from non-nan/false storage values
  */
template <typename tree_t, typename map_t, typename ary_t>
tree_t build_branch(const map_t& storage_map, const ary_t& storage, const bool inc_off=false);

/**
  * get_last_enabled_no_throw
  *
  * from an array of bools, determine the last enabled, true, in the array
  * and determine the associated mapping.  Returns a tuple with the name and index
  * 
  * @tparam ary_t array type
  * @tparam map_t map type
  * @param[input] storage_map: map of boolean parameters
  * @param[input] storage array: of boolean parameters
  * @return       std::tuple tuple of index and associated mapped string
  */
template <typename ary_t, typename map_t>
std::tuple<std::string, int> get_last_enabled(map_t enum_map, ary_t toggle); 

/**
  * get_last_enabled
  *
  * Same as get_last_enabled_no_throw only an exception will throw if an
  * active stage is not found
  * 
  * @tparam ary_t array type
  * @tparam map_t map type
  * @param[input] storage_map: map of boolean parameters
  * @param[input] storage array: of boolean parameters
  * @return       std::tuple tuple of index and associated mapped string
  */
template <typename ary_t, typename map_t>
std::tuple<std::string, int> get_last_enabled(map_t enum_map, ary_t toggle); 

/**
 * append_map
 *
 * Here we check to see if there is an enabled value in the boolean storage array that corresponds\n 
 * to a key/value pair that doesn't already exist in the parial_map.  IFF one is found, we add it\n 
 * to a new map and return it.\n 
 * This is primarily useful for developers so they can have non-standard stages enabled in the\n 
 * config file which aren't deleted after a run.
 *
 * @tparam map_t: type of the map <string, enum>
 * @tparam boolary_t: array of bools <bool, enumsize>
 * @param[input] full_map: Full map of possible elements
 * @param[input] partial_map: subset of elements from full_map
 * @param[input] storage: array of toggles based on full_map
 * @return map: possibly an appended version of partial_map.
 */
template <typename map_t, typename boolary_t>
map_t append_map (const map_t& full_map, const map_t& partial_map, const boolary_t& storage);


template <class tree_t>
auto find_leaf(tree_t const & tree, std::string key) {
  
  std::string branch_name{};

  auto search_branch = [&](auto& branch) -> std::tuple<std::string, std::string> {
    for(auto const & node : branch) {
      if(node.second.empty())
        if(node.first == key) {
          return std::make_tuple(node.first, node.second.data());
        }
    }
    return {};
  };
  std::function<std::tuple<std::string,std::string,std::string>(tree_t const &)> recursive_search;
  recursive_search = [&](auto& branch) -> std::tuple<std::string,std::string,std::string> {
    auto tmp = search_branch(branch);

    if(!std::get<0>(tmp).empty()) {  
      #ifdef DEBUG
      std::cout << branch_name << ", " << std::get<0>(tmp) << "," << std::get<1>(tmp) << std::endl;
      #endif
      return std::make_tuple(branch_name, std::get<0>(tmp), std::get<1>(tmp));
    }
    for(auto const & node : branch) {
      if(!node.second.empty()) {
        branch_name = node.first;
        tree_t t_branch = read_branch(branch, branch_name);
        auto tmp_tuple = recursive_search(t_branch);
        if(!std::get<0>(tmp_tuple).empty()) {
          return tmp_tuple;
        }
      }
    }
    return {};
  };
  auto res = recursive_search(tree);
  
  return res;
}

/**
 * @brief Checks if all storage elements are std::nan.  Works
 * for fundamental and variant types
 * 
 * @tparam ary_t Template argument for storage array
 * @param storage storage array
 * @return true if all elements are std::nan
 * @return false 
 */
template<class ary_t>
constexpr inline bool is_storage_all_nan(ary_t& storage);

/**
 * @brief Recursively add branch data and update the data
 * path
 * 
 * @tparam tree_t Boost tree type by default
 * @param tree Tree to add data to
 * @param branch Branch to read data from
 * @param path base path to add data to
 */
template <typename tree_t>
void add_branch_data(tree_t& tree, tree_t& branch, std::string path = "");

/** @} end config_utils group */

#include "config_utils_boost_imp.cpp"
}}
