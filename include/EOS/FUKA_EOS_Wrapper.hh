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
#include "standalone/Margherita_EOS.h"
#include "standalone/cold_pwpoly.hh"
#include "standalone/cold_pwpoly_implementation.hh"
#include "standalone/cold_table.hh"
#include "standalone/cold_table_implementation.hh"
#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#endif

namespace Kadath {
namespace FUKA_EOS {

typedef enum {
  margherita_pwp,
  margherita_1d,
  ghl_eos_simple,
  ghl_eos_hybrid,
  ghl_eos_tabulated
} fuka_eos_t;

template <class FUKA_eos_t, FUKA_eos_t eos>
struct FUKA_EOS_Wrapper {
#ifdef WITH_GRHAYL_EOS
  static std::unique_ptr<ghl_eos_parameters> ghl_eos_params;
#endif

  static inline double rho__press_cold(const double P_in) {
    double rho;
    if constexpr (eos == margherita_pwp) {
      double press = P_in;
      Kadath::Margherita::Cold_PWPoly::error_t error;
      rho = Kadath::Margherita::Cold_PWPoly::rho__press_cold(press, error);
    } else if constexpr (eos == margherita_1d) {
      double press = P_in;
      Kadath::Margherita::Cold_Table::error_t error;
      rho = Kadath::Margherita::Cold_Table::rho__press_cold(press, error);
    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      rho = ghl_hybrid_compute_rho_cold_from_P_cold(ghl_eos_params.get(), P_in);
    } else if constexpr (eos == ghl_eos_tabulated) {
      rho = ghl_tabulated_compute_rho_from_P(ghl_eos_params.get(), P_in);
    }
#endif
    return rho;
  }

  static inline double rho_energy__rho_cold(const double rhoB_in) {
    double eps;
    if constexpr (eos == margherita_pwp) {
      double press;
      double rho = rhoB_in;
      Kadath::Margherita::Cold_PWPoly::error_t error;
      press =
          Kadath::Margherita::Cold_PWPoly::press_cold_eps_cold__rho(eps, rho,
                                                                    error);
    } else if constexpr (eos == margherita_1d) {
      double press;
      double rho = rhoB_in;
      Kadath::Margherita::Cold_Table::error_t error;
      press = Kadath::Margherita::Cold_Table::press_cold_eps_cold__rho(eps, rho,
                                                                       error);
    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      double P;
      ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos_params.get(), rhoB_in, &P,
                                             &eps);
    } else if constexpr (eos == ghl_eos_tabulated) {
      eps = ghl_tabulated_compute_eps_from_rho(ghl_eos_params.get(), rhoB_in);
    }
#endif
    return rhoB_in * (1.0 + eps);
  }

  static inline double dpress_cold_drho__rho(const double rhoB_in) {
    double dpress_cold_drho;
    if constexpr (eos == margherita_pwp) {
      double rho = rhoB_in;
      Kadath::Margherita::Cold_PWPoly::error_t error;
      dpress_cold_drho =
          Kadath::Margherita::Cold_PWPoly::dpress_cold_drho__rho(rho, error);
    } else if constexpr (eos == margherita_1d) {
      double rho = rhoB_in;
      Kadath::Margherita::Cold_Table::error_t error;
      dpress_cold_drho =
          Kadath::Margherita::Cold_Table::dpress_cold_drho__rho(rho, error);
    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      double K, Gamma, eps, press_cold;
      ghl_hybrid_get_K_and_Gamma(ghl_eos_params.get(), rhoB_in, &K, &Gamma);

      ghl_hybrid_compute_P_cold(ghl_eos_params.get(), rhoB_in, &press_cold);
      dpress_cold_drho = Gamma * press_cold / rhoB_in;
    } else if constexpr (eos == ghl_eos_tabulated) {
      dpress_cold_drho =
          ghl_tabulated_compute_dP_drho_from_rho(ghl_eos_params.get(), rhoB_in);
    }
#endif
    return dpress_cold_drho;
  }

  static inline double dedp__P_cold(const double rhoB_in,
                                    const double rho_energy_in,
                                    const double P_in) {
    double dpdrho = dpress_cold_drho__rho(rhoB_in);

    auto const rhoh = rho_energy_in + P_in;

    return rhoh / (dpdrho * rhoB_in);
  }

  static inline double P_cold_from_rho(const double rhoB_in) {
    double P;
    if constexpr (eos == margherita_pwp) {
      double rho = rhoB_in;
      double eps;
      Kadath::Margherita::Cold_PWPoly::error_t error;
      P = Kadath::Margherita::Cold_PWPoly::press_cold_eps_cold__rho(eps, rho,
                                                                    error);
    } else if constexpr (eos == margherita_1d) {
      double rho = rhoB_in;
      double eps;
      Kadath::Margherita::Cold_Table::error_t error;
      P = Kadath::Margherita::Cold_Table::press_cold_eps_cold__rho(eps, rho,
                                                                   error);
    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      ghl_hybrid_compute_P_cold(ghl_eos_params.get(), rhoB_in, &P);
    } else if constexpr (eos == ghl_eos_tabulated) {
      P = ghl_tabulated_compute_P_from_rho(ghl_eos_params.get(), rhoB_in);
    }
#endif
    return P;
  }

  // Wrapper to make API consistent without modifying Margherita
  static inline double press_cold_eps_cold__rho(double& eps_cold,
                                                double& rhoB_in) {
    double P;
    if constexpr (eos == margherita_pwp) {
      Kadath::Margherita::Cold_PWPoly::error_t error;
      P = Kadath::Margherita::Cold_PWPoly::press_cold_eps_cold__rho(eps_cold,
                                                                    rhoB_in,
                                                                    error);
    } else if constexpr (eos == margherita_1d) {
      Kadath::Margherita::Cold_Table::error_t error;
      P = Kadath::Margherita::Cold_Table::press_cold_eps_cold__rho(eps_cold,
                                                                   rhoB_in,
                                                                   error);
    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos_params.get(), rhoB_in, &P,
                                             &eps_cold);
    } else if constexpr (eos == ghl_eos_tabulated) {
      P = ghl_tabulated_compute_P_from_rho(ghl_eos_params.get(), rhoB_in);
      eps_cold =
          ghl_tabulated_compute_eps_from_rho(ghl_eos_params.get(), rhoB_in);
    }
#endif
    return P;
  }

  static inline double rho__h_cold(double& h_in) {
    double rho;
    if constexpr (eos == margherita_pwp) {
      Kadath::Margherita::Cold_PWPoly::error_t error;
      rho = Kadath::Margherita::Cold_PWPoly::rho__h_cold(h_in, error);
    } else if constexpr (eos == margherita_1d) {
      Kadath::Margherita::Cold_Table::error_t error;
      rho = Kadath::Margherita::Cold_Table::rho__h_cold(h_in, error);
    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      rho = ghl_hybrid_compute_rho_cold_from_h(ghl_eos_params.get(), h_in);
    } else if constexpr (eos == ghl_eos_tabulated) {
      // GRHayL will throw an error rather than enforce table bounds
      // so we enforce them here instead.
      size_t const N = ghl_eos_params->N_rho;
      const double h_min = std::exp(ghl_eos_params->lh_of_lr[0]);
      const double h_max = std::exp(ghl_eos_params->lh_of_lr[N-1]);
      if(h_in < h_min) {
        h_in = h_min;
        return std::exp(ghl_eos_params->table_logrho[0]);
      } else if(h_in > h_max) {
        h_in = h_max;
        return std::exp(ghl_eos_params->table_logrho[N-1]);
      }

      // It's only safe to linear interpolate in logp vs logrho space
      // Instead we rootfind the correct specific enthalpy rather than
      // attempting a table look-up using the same method as used
      // in Margherita::Cold_Table
      const auto func = [&](const double &lrho) {
        double const rho = std::exp(lrho);
        double const eps = ghl_tabulated_compute_eps_from_rho(ghl_eos_params.get(), rho);
        double const press = ghl_tabulated_compute_P_from_rho(ghl_eos_params.get(), rho);

        return h_in - ( 1. + eps + press/rho);
      };

      // Note the root bounds are quite sensitive given GRHaYL will simply error out
      // rather than enforcing table bounds.
      auto lrho = zero_brent<>(log(ghl_eos_params->rho_min * 1.001), log(ghl_eos_params->rho_max), 1.0e-13, func);
      return std::exp(lrho);
    }
#endif
    return rho;
  }

  static inline double rho_energy_dedp__P_cold(double& rhoE,
                                               double& dedp,
                                               double& P_in) {
    double rho = rho__press_cold(P_in);
    rhoE = rho_energy__rho_cold(rho);
    auto const rhoh = rhoE + P_in;
    auto const dpdrho = dpress_cold_drho__rho(rho);
    dedp = rhoh / (dpdrho * rho);
    return rho;
  }
};
#ifdef WITH_GRHAYL_EOS
template <class FUKA_eos_t, FUKA_eos_t eos>
std::unique_ptr<ghl_eos_parameters>
    FUKA_EOS_Wrapper<FUKA_eos_t, eos>::ghl_eos_params = nullptr;
#endif
}  // namespace FUKA_EOS
}  // namespace Kadath