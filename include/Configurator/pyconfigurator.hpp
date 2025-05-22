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
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <string>
#include "config_bco.hpp"
#include "config_binary.hpp"

namespace Kadath {
namespace FUKA_pyTools {
using namespace Kadath::FUKA_Config;

/** @defgroup pyConfigurator
 * @ingroup Configurator
 * Utility functions to enable parsing and saving Configurator files currently
 * based on the BOOST proprietary INFO format.  In principle this can be changed easily
 * by changing the include to another BOOST parser, however, this is untested.
 * @{
 */

/**
 * @brief Build configuration dictionary for python
 *
 * @tparam bconfig_t Configurator type
 * @tparam idx_t Optional pack
 * @param bconfig the Configurator
 * @param bco Optional index for a subobject
 * @return boost::python::dict
 */
template <class bconfig_t, typename... idx_t>
boost::python::dict build_config_dict(bconfig_t& bconfig, idx_t... bco) {
  auto map = bconfig.get_map(bco...);
  boost::python::dict config;
  for (auto& tup : map) {
    auto name{tup.first};
    auto index{tup.second};
    if (!std::isnan(bconfig.set(index, bco...)))
      config[name] = bconfig(index, bco...);
  }
  return config;
}

/**
 * @brief Build EOS configuration dictionary for python
 *
 * @tparam bconfig_t Configurator type
 * @tparam idx_t Optional pack
 * @param bconfig the Configurator
 * @param bco Optional index for a subobject
 * @return boost::python::dict
 */
template <class bconfig_t, typename... idx_t>
void add_eos_todict(bconfig_t& bconfig,
                    boost::python::dict& config,
                    idx_t... bco) {
  // map of <std::string, enum>
  auto map = bconfig.get_eos_map(bco...);
  for (auto& tup : map) {
    auto name{tup.first};
    auto index{tup.second};

    auto assign_val = [&](auto&& v) {
      using v_t = std::decay_t<decltype(v)>;
      if constexpr (std::is_same_v<v_t, double>) {
        if (!std::isnan(v))
          config[name] = v;
      } else {
        config[name] = v;
      }
    };

    // EOS Params are stored in a std::variant
    std::visit(assign_val, bconfig.set_eos(index, bco...));
  }
  return;
}

/**
 * @brief Build Differential rotation configuration dictionary for python
 *
 * @tparam bconfig_t Configurator type
 * @tparam idx_t Optional pack
 * @param bconfig the Configurator
 * @param bco Optional index for a subobject
 * @return boost::python::dict
 */
template <class bconfig_t, typename... idx_t>
void add_diffrot_todict(bconfig_t& bconfig,
                        boost::python::dict& config,
                        idx_t... bco) {
  // map of <std::string, enum>
  boost::python::dict diffrot;
  auto map = bconfig.get_diffrot_map(bco...);
  for (auto& tup : map) {
    auto name{tup.first};
    auto index{tup.second};

    auto assign_val = [&](auto&& v) {
      using v_t = std::decay_t<decltype(v)>;
      if constexpr (std::is_same_v<v_t, double>) {
        if (!std::isnan(v))
          diffrot[name] = v;
      } else {
        diffrot[name] = v;
      }
    };

    // Differential Rotation Params are stored in a std::variant
    std::visit(assign_val, bconfig.set_diffrot(index, bco...));
  }
  if (boost::python::len(diffrot) > 0)
    config["differential_rotation"] = diffrot;
  return;
}

/**
 * @brief Base class for the python Configurator reader
 *
 * @tparam config_t
 */
template <class config_t>
struct Configurator_reader_t {
  Configurator_reader_t(std::string const filename_)
      : filename(filename_), bconfig(config_t(filename)) {
    add_dict();
  }

  // Keep this public to make it easy to access
  boost::python::dict config;

 protected:
  /**
   * @brief Add a dictionary to this->config
   *
   * @tparam NS If a NS is present we extract EOS info
   * @tparam idx_t index pack
   * @param bco optional index (for binary configs)
   */
  template <bool NS = false, typename... idx_t>
  void add_dict(idx_t... bco) {
    std::string dict_name{bconfig.get_name_string(bco...)};
    boost::python::dict tmp = build_config_dict(bconfig, bco...);
    if constexpr (NS)
      add_eos_todict(bconfig, tmp, bco...);
    config[dict_name] = tmp;
  };

  std::string filename{};
  config_t bconfig;
};

using bin_config_t = kadath_config_boost<BIN_INFO>;

/**
 * @brief Binary specialization Configurator reader
 *
 */
struct bin_Configurator_reader_t : Configurator_reader_t<bin_config_t> {
  bin_Configurator_reader_t(std::string const filename_)
      : Configurator_reader_t(filename_) {
    // Base constructor only extracts binary parameters
    // Here we add the seperate configurations for each compact object
    if (bconfig.get_name_string(BCO1)[0] == 'n')
      add_dict<true>(BCO1);
    else
      add_dict(BCO1);

    if (bconfig.get_name_string(BCO2)[0] == 'n')
      add_dict<true>(BCO2);
    else
      add_dict(BCO2);
  }
};

using ns_config_t = kadath_config_boost<BCO_NS_INFO>;

/**
 * @brief Neutron star specialization Configurator reader
 *
 */
struct ns_Configurator_reader_t : Configurator_reader_t<ns_config_t> {
  ns_Configurator_reader_t(std::string const filename_)
      : Configurator_reader_t(filename_) {
    // The base class only extracts the object parameters
    // Here we extract the dictionary and add the EOS information
    auto dict_name = bconfig.get_name_string();
    boost::python::dict dict =
        boost::python::extract<boost::python::dict>(config[dict_name]);

    add_eos_todict(bconfig, dict);
    add_diffrot_todict(bconfig, dict);
    config[dict_name] = dict;
  }
};

// dummy constructor function, defining readers through boost python
template <typename reader_t>
void constructPythonConfigurator(std::string reader_name) {
  using namespace boost::python;

  auto reader = class_<reader_t>(reader_name.c_str(), init<std::string>());
  reader.def_readonly("config", &reader_t::config);
}
}  // namespace FUKA_pyTools
}  // namespace Kadath

/** @} end pyConfigurator group */
