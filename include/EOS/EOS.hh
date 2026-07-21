/*
 * This file is part of the KADATH library.
 * Copyright (C) 2019, Samuel Tootle
 *                     <tootle@th.physik.uni-frankfurt.de>
 * Copyright (C) 2019, Ludwig Jens Papenfort
 *                      <papenfort@th.physik.uni-frankfurt.de>
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
#include <cstdlib>
#include <string>
#include "FUKA_EOS_Wrapper.hh"
#include "name_tools.hpp"
#include "standalone/cold_pwpoly.hh"
#include "standalone/cold_pwpoly_implementation.hh"
#include "standalone/cold_table.hh"
#include "standalone/cold_table_implementation.hh"
#include "standalone/setup_cold_table.cc"
#include "standalone/setup_polytrope.cc"

#include <scalar.hpp>
#include <term_eq.hpp>
#include <val_domain.hpp>

using namespace Kadath;

/**
 * The various hydrodynamic quantities that can be obtained from this
 * interface for a given EOS
 */
enum eos_var_t { PRESSURE, EPSILON, DENSITY, DHDRHO };

/**
 * This interface provides the user defined OPEs to provide hydrodynamic
 * quantities as a function of the specific enthalpy (\b h) to the Kadath
 * System_of_eqs framework.\n The equation of state is managed by a modified,
 * standalone version of Margherita: https://github.com/fil-grmhd/Margherita-EOS
 *
 * @tparam eos Margherita EOS type (e.g. Cold_PWPoly)
 * @tparam var EOS variable to update when action() is executed. see eos_var_t
 */
template <typename eos, eos_var_t var>
class EOS {
   private:
    /**
   * EOS::term_by_term_variation
   *
   * calculate the numerical variation of the dependent quantity (var), term by
   * term, as a function of a scalar quantity (currently h).
   *
   * @param [input] dom: domain to update
   * @param [input] so: previous numerical variatioon of the variable field (h)
   * @param [input] scalar: current value of the variable field (h)
   */
    static inline Val_domain term_by_term_variation(int dom,
                                                    const Val_domain so,
                                                    const Val_domain scalar) {
        if (so.check_if_zero()) {
            return so;
        }

        // need to work in configuration space
        so.coef_i();
        Val_domain res(so.get_domain());
        res.allocate_conf();
        Index pos(so.get_conf().get_dimensions());

        do {
            double eps_cold = 0.0;
            double h = scalar(pos);
            double dh = so(pos);

            double rho = eos::rho__h_cold(h);
            double dpdrho = eos::dpress_cold_drho__rho(rho);
            double drho = rho * 1. / dpdrho * dh;

            if constexpr (var == EPSILON) {
                double pressure = eos::press_cold_eps_cold__rho(eps_cold, rho);
                res.set(pos) = pressure / pow(rho, 2) * drho;

            } else if constexpr (var == PRESSURE) {
                res.set(pos) = eos::dpress_cold_drho__rho(rho) * drho;

            } else if constexpr (var == DENSITY) {
                res.set(pos) = drho;

            } else if constexpr (var == DHDRHO) {
                res.set(pos) = 0.0;  //-1.0 / rho / rho * dpdrho * drho;
            } else
                std::cerr << "Ill-defined variable in EOS class, please check."
                          << std::endl;
        } while (pos.inc());

        res.set_base() = so.get_base();
        return res;
    }

    /**
   * EOS::term_by_term
   *
   * calculate the value the dependent quantity (var), term by term,
   * as a function of a scalar quantity (currently h).
   *
   * @param [input] dom: domain to update
   * @param [input] so: current value of the variable field (h)
   */
    static inline Val_domain term_by_term(int dom, const Val_domain so) {
        if (so.check_if_zero()) {
            return so;
        }

        // need to work in configuration space
        so.coef_i();
        Val_domain res(so.get_domain());
        res.allocate_conf();
        Index pos(so.get_conf().get_dimensions());

        do {
            double eps_cold = 0.0;
            double h = so(pos);
            double rho = eos::rho__h_cold(h);
            double pressure = eos::press_cold_eps_cold__rho(eps_cold, rho);

            if constexpr (var == EPSILON)
                res.set(pos) = eps_cold;
            else if constexpr (var == DENSITY)
                res.set(pos) = rho;
            else if constexpr (var == PRESSURE)
                res.set(pos) = pressure;
            else if constexpr (var == DHDRHO) {
                res.set(pos) = 1. / h * eos::dpress_cold_drho__rho(rho);
            } else
                std::cerr << "Ill-defined variable in EOS class, please check."
                          << std::endl;
        } while (pos.inc());

        res.set_base() = so.get_base();
        return res;
    }

   public:
    /**
   * EOS::h_cold__rho
   *
   * compute specific enthalpy from a given density.  Used primarily in analysis
   * codes.
   *
   * @param [input] rho: density
   */
    static double h_cold__rho(double rho) {
        double eps_cold = 0.;
        double pressure = eos::press_cold_eps_cold__rho(eps_cold, rho);
        double h = 1. + eps_cold + pressure / rho;

        return h;
    }

    /**
   * EOS::get
   *
   * for a given specific enthalpy, get the dependent quantity(var)
   *
   * @param [input] h: specific enthalpy
   */
    static double get(double h) {
        double eps_cold = 0.0;
        double rho = eos::rho__h_cold(h);
        double pressure = eos::press_cold_eps_cold__rho(eps_cold, rho);

        if constexpr (var == EPSILON)
            return eps_cold;
        else if constexpr (var == DENSITY)
            return rho;
        else if constexpr (var == PRESSURE)
            return pressure;
    }

    /**
   * EOS::action
   *
   * This is called by the System of equations in order to generate the
   * corresponding Scalar field and its variation based on a definition using a
   * user defined OPE.
   *
   * @param [input] term: term to get information from (i.e. specific enthalpy).
   * @param [input] p: Kadath parameter.  Not used, but required for user
   * defined OPEs
   */
    static Term_eq action(const Term_eq& term, Param* p) {
        Term_eq target(term);

        int dom = term.get_dom();
        if (target.get_type_data() != TERM_T) {
            std::cerr << "EOS only defined with respect for a tensor"
                      << std::endl;
            std::_Exit(EXIT_FAILURE);
        }

        Scalar scalar(target.get_val_t());
        Scalar result(scalar, false);

        Val_domain value(scalar(dom));
        if (value.check_if_zero())
            result.set_domain(dom).set_zero();
        else
            result.set_domain(dom) = EOS::term_by_term(dom, value);

        // Check if a derivative has been defined (has the pointer been initialized)
        // if so, calculate the variance for new value_domain
        if (target.get_p_der_t() != 0x0) {
            Scalar scalar_var(target.get_der_t(), true);
            Scalar result_var(scalar_var, false);

            Val_domain value_var(scalar_var(dom));
            if (value.check_if_zero())
                result_var.set_domain(dom).set_zero();
            else
                result_var.set_domain(dom) =
                    EOS::term_by_term_variation(dom, value_var, value);

            return Term_eq(dom, result, result_var);
        } else {
            return Term_eq(dom, result);
        }
    }
};
