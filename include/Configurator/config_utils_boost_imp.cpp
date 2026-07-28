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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <variant>
#include "config_enums.hpp"
#include "config_utils_boost.hpp"

/** \addtogroup Configurator_Utils
  * @{ */

template <typename map_t, typename T>
bool check_for_nan(const map_t& map, const T& var, const int idx) {
    if (!std::isnan(var)) {
        return false;
    } else {
        auto var_name =
            std::find_if(map.begin(), map.end(), [idx](const auto& pair) {
                return pair.second == idx;
            });

        std::string msg;
        if (var_name == map.end())
            msg = "No var found matching index: " + std::to_string(idx);
        else
            msg = "Var \"" + var_name->first + "\" is undefined.";

        throw std::runtime_error(msg.c_str());
    }
    return true;
}

template <typename tree_t>
tree_t read_branch(const tree_t& tree, std::string node) {
    if (tree.find(node.data()) != tree.not_found()) {
        tree_t branch = tree.get_child(node.data());
        return branch;
    } else {
        std::cerr << node << " not found " << std::endl;
        return tree_t{};
    }
}

template <typename ary_t, typename map_t, typename tree_t>
void get_branch_nodes(const map_t& storage_map,
                      ary_t& storage,
                      const tree_t& branch,
                      bool rsuffix) {
    size_t idx = 0;
    for (const auto& node : branch) {
        if (!node.second.empty()) {
            std::string tstr = node.first;
            tstr = (rsuffix) ? tstr.substr(0, 2) : tstr;
            const auto& it = storage_map.find(tstr);
            if (it != storage_map.end()) {
                if (idx < storage.size()) {
                    storage[idx] = node.first;
                    idx++;
                } else {
                    std::cerr
                        << "Err: index out of range while gathering nodes."
                        << std::endl;
                    std::_Exit(EXIT_FAILURE);
                }
            } else {
                std::cerr << "Err: " << tstr << " not a recognized node"
                          << std::endl;
                std::_Exit(EXIT_FAILURE);
            }
        }
    }
}

template <typename ary_t, typename map_t, typename tree_t>
void read_keys(const map_t& storage_map,
               ary_t& storage,
               const tree_t& tree) noexcept {
    for (const auto& ele : storage_map) {
        if (tree.find(ele.first.data()) != tree.not_found()) {
            auto key = tree.get_child(ele.first.data()).data();
            int idx = ele.second;
            using var_t = decltype(storage[idx]);
            if (key == "on" || key == "off") {
                if (std::is_assignable<var_t, bool>::value)
                    storage[idx] = (key == "on") ? true : false;
                else
                    std::cerr
                        << "Boolean parameter " << ele.first
                        << " encountered, but cannot assign to storage array\n";
                continue;
            }
            try {
                double td = std::stod(key);
                int ti = std::stoi(key);

                //prefer to store int over double if assignable to array
                if ((td - ti) == 0 && std::is_assignable<var_t, int>::value)
                    storage[idx] = ti;
                else if (std::is_assignable<var_t, double>::value)
                    storage[idx] = td;
                else
                    std::cerr << "int/double parameter " << ele.first
                              << " encountered, but cannot assign to storage "
                                 "array \n";
                continue;
            } catch (...) {}  //catch typecasting error from stod/stoi
            if constexpr (std::is_assignable<var_t, std::string>::value) {
                std::string val = key;
                storage[idx] = key;
                continue;
            } else
                std::cerr
                    << "String parameter " << ele.first
                    << " encountered, but cannot assign to storage array \n";
        }
    }
    return;
}

template <typename ary_t, typename map_t, typename tree_t>
void read_keys_as_string(const map_t& storage_map,
                         ary_t& storage,
                         const tree_t& tree) noexcept {
    for (const auto& ele : storage_map) {
        if (tree.find(ele.first.data()) != tree.not_found()) {
            auto key = tree.get_child(ele.first.data()).data();
            int idx = ele.second;

            std::string val = key;
            storage[idx] = key;
        }
    }
    return;
}

template <typename map_t, typename ary_t>
void print_params(const map_t& storage_map,
                  const ary_t& storage,
                  std::ostream& out) {
    using constvar_t =
        typename std::remove_pointer<decltype(storage.begin())>::type;
    using var_t = typename std::remove_cv<constvar_t>::type;
    for (const auto& a : storage_map) {
        auto print = [&](auto&& arg) {
            using CONSTTYPE =
                typename std::remove_reference<decltype(arg)>::type;
            using ARGTYPE = typename std::remove_cv<CONSTTYPE>::type;
            if constexpr (std::is_same<ARGTYPE, double>::value) {
                if (std::isnan(arg)) {
                    return;
                }
            } else if constexpr (std::is_same_v<var_t, bool>) {
                if (arg) {
                    out << std::setw(20) << a.first << ": on" << std::endl;
                }
                return;
            }
            out << std::setw(20) << a.first << ": " << arg << std::endl;
        };
        const auto& ele = storage.at(a.second);
        //Checks if the type stored in Array is fundamental - i.e. not a std::variant
        if constexpr (!std::is_fundamental<var_t>::value &&
                      !std::is_same_v<var_t, std::string>)
            std::visit(print, ele);
        else
            print(ele);
    }
    return;
}

template <typename tree_t, typename map_t, typename ary_t>
tree_t build_branch(const map_t& storage_map,
                    const ary_t& storage,
                    const bool inc_off) {
    using constvar_t =
        typename std::remove_pointer<decltype(storage.begin())>::type;
    using var_t = typename std::remove_cv<constvar_t>::type;
    tree_t branch;
    for (const auto& a : storage_map) {
        // Lambda expression is necessary for use with std::variant since std::visit does
        // not allow passing parameters other than variants.  Therefore, a Lambda is used
        // to be able to capture 'a' by reference without having to pass it explicitly
        auto add_key = [&](auto&& arg) -> void {
            using CONSTTYPE =
                typename std::remove_reference<decltype(arg)>::type;
            using ARGTYPE = typename std::remove_cv<CONSTTYPE>::type;
            if constexpr (std::is_same<ARGTYPE, double>::value) {
                if (!std::isnan(arg))
                    branch.put(a.first, arg);
                return;
            } else if constexpr (std::is_same<ARGTYPE, bool>::value) {
                if (arg)
                    branch.put(a.first, "on");
                else if (!arg && inc_off)
                    branch.put(a.first, "off");
                return;
            } else
                branch.put(a.first, arg);
        };
        //Checks if the type stored in Array is fundamental - i.e. not a std::variant
        if constexpr (!std::is_fundamental<var_t>::value &&
                      !std::is_same_v<var_t, std::string>)
            std::visit(add_key, storage[a.second]);
        else
            add_key(storage[a.second]);
    }
    return branch;
}

template <typename map_t, typename val_t>
std::tuple<std::string, int> get_key_val_pair_from_val(map_t enum_map,
                                                       val_t val) {

    auto name = std::find_if(enum_map.begin(),
                             enum_map.end(),
                             [val](const auto& kvpair) {
                                 return kvpair.second == val;
                             });

    if (name != enum_map.end()) {
        return std::make_tuple(name->first, name->second);
    }

    throw std::invalid_argument("\nNo valid key-value pair found\n)");
}

template <typename ary_t, typename map_t>
std::tuple<std::string, int> get_last_enabled_no_throw(map_t enum_map,
                                                       ary_t toggle) {
    auto last_enabled =
        std::distance(std::find(toggle.rbegin(), toggle.rend(), 1),
                      toggle.rend() - 1);
    auto name = std::find_if(enum_map.begin(),
                             enum_map.end(),
                             [last_enabled](const auto& map) {
                                 return map.second == last_enabled;
                             });
    if (name != enum_map.end()) {
        return std::make_tuple(name->first, last_enabled);
    }
    return {};
}

template <typename ary_t, typename map_t>
std::tuple<std::string, int> get_last_enabled(map_t enum_map, ary_t toggle) {
    auto res_tuple = get_last_enabled_no_throw(enum_map, toggle);
    if (!std::get<0>(res_tuple).empty()) {
        return res_tuple;
    }
    throw std::invalid_argument("\nNo valid stages enabled\n)");
}

template <typename map_t, typename boolary_t>
map_t append_map(const map_t& full_map,
                 const map_t& partial_map,
                 const boolary_t& storage) {
    map_t map{partial_map};

    for (const auto& [str, i] : full_map) {
        if (map.find(str) == map.end() && storage[i])
            map.insert({str, i});
    }
    return map;
}

template <class ary_t>
constexpr inline bool is_storage_all_nan(ary_t& storage) {
    using constvar_t =
        typename std::remove_pointer<decltype(storage.begin())>::type;
    using var_t = typename std::remove_cv<constvar_t>::type;
    auto check = [&](auto&& arg) {
        using CONSTTYPE = typename std::remove_reference<decltype(arg)>::type;
        using ARGTYPE = typename std::remove_cv<CONSTTYPE>::type;
        if constexpr (std::is_same<ARGTYPE, double>::value) {
            if (std::isnan(arg)) {
                return true;
            }
        }
        return false;
    };
    for (auto& ele : storage) {
        bool isnan = true;
        if constexpr (!std::is_fundamental<var_t>::value) {
            isnan = std::visit(check, ele);
        } else if constexpr (std::is_fundamental<var_t>::value) {
            isnan = check(ele);
        }
        if (!isnan)
            return false;
    }
    return true;
}

template <typename tree_t>
void add_branch_data(tree_t& tree, tree_t& branch, std::string path) {
    for (auto& key : branch) {
        // Update path
        std::string news = path + "." + key.first;
        if (key.second.empty())
            tree.put(news.c_str(), key.second.data());
        else
            // new_tree.put_child(news.c_str(), key.second);
            add_branch_data(tree, key.second, news);
    }
}

/**
  * @} */
