/*
 * Copyright 2023
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
#include "Configurator/configurator_base.hpp"
#include "Configurator/configurator_boost.hpp"
#include <vector>
#include <functional>
/**
 * \addtogroup Sequences
 * \ingroup FUKA
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

// Forward declarations for ostream operators
template<class... Ts>
struct Parameter_sequence ;

template<class... Ts>
std::ostream& operator<<(std::ostream& out, const Parameter_sequence<Ts...>& Seq);

/**
 * @brief Base class for storing sequences for different parameters
 * 
 */
struct Parameter_sequence_base {
  protected:
  /// Config parameter string (P)
  std::string parameter_str{};
  /// Config parameter value (V)
  double parameter_val{std::nan("1")};
  /// Sequence initial value (P_init)
  double seqinit{std::nan("1")};
  /// Sequence final value (P_final)
  double seqfinal{std::nan("1")};

  public:
  Parameter_sequence_base() = default;
  Parameter_sequence_base(Parameter_sequence_base const &) = default;
  Parameter_sequence_base(std::string _str) : parameter_str(_str) {}
  /// Determine if a sequence has been initialized
  bool is_set() const {
      return !std::isnan(seqinit) && !std::isnan(seqfinal);
  }
  /// Determine if a default value is initialized
  bool is_default_set() const {
      return !std::isnan(parameter_val);
  }

  /// Getters
  double const & init() const { return seqinit;}
  double const & final() const { return seqfinal;}
  double const & default_val() const { return parameter_val;}
  std::string const & str() const { return parameter_str; }
};

/**
 * @brief Parameter sequence for a given set of Config indices
 * 
 * @tparam Ts Parameter pack storing Config index types
 */
template<class... Ts>
struct Parameter_sequence : public Parameter_sequence_base {
  protected:
  /// Config indices for the intended parameter to sequence
  std::tuple<Ts...> parameter_indices;
  /// Conditional function to determine loop criteria
  std::function<bool(double const &, double const &)> conditional;
  /// Number of solutions in sequence
  int N{1};
  
  public:
  
  Parameter_sequence() = default;
  Parameter_sequence(Parameter_sequence const &) = default;
  /**
   * @brief Construct a new Parameter_sequence object from a string and tuple objects
   * 
   * @param _str Input parameter string
   * @param _t Tuple containing Config indicies for the desired parameter
   */
  Parameter_sequence(std::string _str, std::tuple<Ts...> _t) : 
    Parameter_sequence_base(_str), parameter_indices(_t) {}
  
  /**
   * @brief Construct a new Parameter_sequence object from a string and a pack of indices
   * 
   * @param _str Input parameter string
   * @param _ts Parameter pack consisting of the desired Config indices
   */
  Parameter_sequence(std::string _str, Ts... _ts) : 
    Parameter_sequence_base(_str), parameter_indices(std::make_tuple(_ts...)) {}
  
  /**
   * @brief Set default, initial, and final values
   * 
   * @param _val Default value - ignored if is_set()
   * @param _init Initial value
   * @param _final Final value
   */
  void set(double _val, double _init, double _final) {
    parameter_val = _val;
    seqinit = _init;
    seqfinal = _final;
    
    /// Determine the conditional based on whether the
    /// sequence needs to increase or decrease
    if(is_set()) {
      if (seqinit < seqfinal) 
        conditional = std::less_equal<double>{};
      else
        conditional = std::greater_equal<double>{};
    }
  }

  /// Set number of sequences
  void set_N(int _N) { 
    assert(N > 0);
    N = _N;
  }

  /// Compute sequence step_size
  double step_size() const {
    if(is_set())
      return (seqfinal-seqinit) / N;
    return std::nan("1");
  }

  /// Getters
  int const & iterations() const { return N; }
  std::tuple<Ts...> const & get_indices() const { return parameter_indices; }
  /// Evaluate conditional to determine, e.g. if a loop should end
  bool loop_condition(double const & val) const { return conditional(val, seqfinal); }  
  /// Formatted output
  friend std::ostream &operator<< <Ts...>(std::ostream &, const Parameter_sequence<Ts...> &);
};

/**
 * @brief Resolution sequence
 * 
 * @tparam ndom Number of domains
 * @tparam Ts Parameter pack storing Config index types
 */
template<int ndom, class... Ts>
struct Resolution_sequence : public Parameter_sequence_base {
  using Pary_t = std::array<int, 3>;
  using Cary_t = std::array<Pary_t, ndom>;
  protected:
  std::tuple<Ts...> parameter_indices;
  Cary_t resolution;

  public:
  Resolution_sequence() = default;

  /**
   * @brief Construct a new Parameter_sequence object from a string and tuple objects
   * 
   * @param _str Input parameter string
   * @param _t Tuple containing Config indicies for the desired parameter
   */
  Resolution_sequence(std::string _str, std::tuple<Ts...> _t) : 
    Parameter_sequence_base(_str), parameter_indices(_t) {}

  /**
   * @brief Construct a new Parameter_sequence object from a string and a pack of indices
   * 
   * @param _str Input parameter string
   * @param _ts Parameter pack consisting of the desired Config indices
   */
  Resolution_sequence(std::string _str, Ts... _ts) : 
    Parameter_sequence_base(_str), parameter_indices(std::make_tuple(_ts...)) {}

  /// Getter
  std::tuple<Ts...> const & get_indices() const { return parameter_indices; }
  /// Formatted output
  friend std::ostream &operator<< <ndom, Ts...>(std::ostream &, const Resolution_sequence<ndom, Ts...> &);
};
/** @}*/
}}