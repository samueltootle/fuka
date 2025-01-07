/*
 * This file is part of the KADATH library.
 * Copyright (C) 2024, Samuel Tootle
 *                     <tootle@th.physik.uni-frankfurt.de>
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
#include <memory>
#include "Configurator/config_enums.hpp"
#include "EOS/EOS.hh"
#include "standalone/cold_pwpoly.hh"
#include "standalone/cold_pwpoly_implementation.hh"
#include "standalone/cold_table.hh"
#include "standalone/cold_table_implementation.hh"
#include "standalone/setup_cold_table.cc"
#include "standalone/setup_polytrope.cc"
#include "system_of_eqs.hpp"

#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#include "ghl_eos_helpers/ghl_hybrid_helpers.hpp"
#endif

namespace Kadath {
namespace FUKA_EOS {

/**
 * @brief Setup the EOS operators in System_of_eqs
 *
 * @tparam eos_t FUKA_EOS_Wrapper type
 * @param syst System of equations
 * @param p Parameter container (unused, but required by add_ope)
 */
template <class eos_t>
struct set_eos_ope_struct {
  void operator()(System_of_eqs& syst, Param& p) {
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  }
};

/**
 * @brief Wrapper struct for calling a templated functor
 *
 * @tparam T class type for the Functor
 */
template <class T>
struct functor_wrapper {
  /**
   * @brief Functor that calls a templated Functor F with optional arguments
   * Args...
   *
   * @tparam F templated Functor
   * @tparam Args... type(s) of optional parameter pack of Functor arguments
   * @param args optional parameter pack of Functor arguments
   */
  template <template <typename> class F, typename... Args>
  auto operator()(Args&&... args) const {
    return F<T>()(std::forward<Args>(args)...);
  }
};

using namespace Kadath::FUKA_Config;

/**
 * @brief Dispatcher class: Encapsulates all decision-making logic for various
 * EOS'
 */
struct EOS_Function_Dispatcher {
  /**
   * @brief Launches a templated Functor that requires knowledge of the
   * EOS type (module) being used to reduce code duplication.
   *
   * @tparam F templated Functor
   * @tparam Args... type(s) of optional parameter pack of Functor
   * arguments
   * @param args optional parameter pack of Functor arguments
   */
  template <template <typename> class F, class config_t, typename... Args>
  static inline auto dispatch(config_t& bconfig,
                              const std::string eos_type,
                              Args&&... args) {
    if (eos_type == "Cold_PWPoly") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, margherita_pwp>;
      functor_wrapper<eos_t> wrapper;
      return wrapper.template operator()<F>(std::forward<Args>(args)...);
    } else if (eos_type == "Cold_Table") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, margherita_1d>;
      functor_wrapper<eos_t> wrapper;
      return wrapper.template operator()<F>(std::forward<Args>(args)...);
    }
#ifdef WITH_GRHAYL_EOS
    else if (eos_type == "grhayl_eos_tabulated") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, ghl_eos_tabulated>;
      functor_wrapper<eos_t> wrapper;
      return wrapper.template operator()<F>(std::forward<Args>(args)...);
    } else if (eos_type == "grhayl_eos_hybrid") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, ghl_eos_hybrid>;
      functor_wrapper<eos_t> wrapper;
      return wrapper.template operator()<F>(std::forward<Args>(args)...);
    }
#endif
    throw std::invalid_argument("\nInvalid EOS type\n");
  }
};

// Simple struct to simplify EOS initialization and remove code duplication.
struct EOS_initialize {
  template <class config_t, class... BCO_t>
  static inline auto init(config_t& bconfig, BCO_t... bco) {
    using namespace ::Kadath::FUKA_Config;
    using namespace ::Kadath::FUKA_EOS;
    using namespace ::Kadath::Margherita;

    const double h_cut = bconfig.template eos<double>(EOS_PARAMS::HCUT, bco...);
    std::string filename =
        bconfig.template eos<std::string>(EOS_PARAMS::EOSFILE, bco...);
    const std::string eos_type =
        bconfig.template eos<std::string>(EOS_PARAMS::EOSTYPE, bco...);

    auto get_default_path = [&]() {
      std::string default_path{"./"};
      const std::string kadath_environment_var{"HOME_KADATH"};
      if (std::getenv(kadath_environment_var.c_str())) {
        std::string const home_kadath{
            std::getenv(kadath_environment_var.c_str())};
        default_path = home_kadath + "/eos/";
      }
      return default_path;
    };
    std::string const default_path{get_default_path()};

    // if no path is given, we set the default EOS diretory to look for the
    // relevant table/polytrope
    if (filename.rfind("/") == std::string::npos) {
      filename = default_path + filename;
    }

    if (eos_type == "Cold_PWPoly") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, margherita_pwp>;
      return Margherita_setup_polytrope(filename);
    } else if (eos_type == "Cold_Table") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, margherita_1d>;
      const int interp_pts =
          (bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS, bco...) == 0)
              ? 2000
              : bconfig.template eos<int>(EOS_PARAMS::INTERP_PTS, bco...);

      setup_Cold_Table(filename, interp_pts, h_cut);
      // std::cout << Kadath::Margherita::Cold_Table::rhomin << std::endl;
      // std::cout << Kadath::Margherita::Cold_Table::rhomax << std::endl;
      // std::cout << Kadath::Margherita::Cold_Table::press_min << std::endl;
      // std::cout << Kadath::Margherita::Cold_Table::press_max << std::endl;
      // std::cout << Kadath::Margherita::Cold_Table::hmin << std::endl;
      // std::cout << Kadath::Margherita::Cold_Table::hmax << std::endl;

      return;
    }
#ifdef WITH_GRHAYL_EOS
    else if (eos_type == "grhayl_eos_tabulated") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, ghl_eos_tabulated>;

      using namespace ::Kadath::GHL_EOS;
      auto eos_params = ghl_setup_table(filename);
      eos_t::ghl_eos_params = std::move(eos_params);

      return;
    } else if (eos_type == "grhayl_eos_hybrid") {
      using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, ghl_eos_hybrid>;
      using namespace ::Kadath::GHL_EOS;
      auto eos_params = populate_ghl_polytrope(filename);
      eos_t::ghl_eos_params = std::move(eos_params);
      return;
    }
#endif
    throw std::invalid_argument("\nInvalid EOS type\n");
  }
};
}  // namespace FUKA_EOS
}  // namespace Kadath