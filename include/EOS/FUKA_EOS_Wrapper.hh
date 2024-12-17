#include <memory>
#include "standalone/cold_pwpoly.hh"
#include "standalone/cold_pwpoly_implementation.hh"
#include "standalone/cold_table.hh"
#include "standalone/cold_table_implementation.hh"
#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#endif

namespace Kadath {
namespace FUKA_EOS {

typedef enum {margherita_pwp, margherita_1d} margherita_eos_t;

template<class FUKA_eos_t, FUKA_eos_t eos>
struct FUKA_EOS_Wrapper {
  static std::unique_ptr<ghl_eos_parameters> ghl_eos_params;

  static inline void rho_cold__P_cold(
    const double rho_in,
    const double P_in,
    double *restrict rhoB_cold_ptr) {

    if constexpr(eos == margherita_pwp) {
        Cold_PWPoly::error_t error;
        Cold_PWPoly::rho__press_cold()
    } else if constexpr(eos == margherita_1d) {

    }
#ifdef WITH_GRHAYL_EOS
    else if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      // Retrieve K and Gamma from GRHayL
      double K, Gamma, eps, P;
      ghl_hybrid_get_K_and_Gamma(ghl_eos_params.get(), rho_in, &K, &Gamma);

      *rhoB_cold_ptr = pow(P_in / K, 1.0 / Gamma);
    } else if constexpr (eos == ghl_eos_tabulated) {
      *rhoB_cold_ptr = ghl_tabulated_compute_rho_from_P(ghl_eos_params.get(), P_in);
    }
#endif
  }

  static inline void rho_energy__rho_cold(
    const double rhoB_in,
    double *restrict rho_energy_ptr) {
    double eps;

    if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      double P;
      ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos_params.get(), rhoB_in,
                                              &P, &eps);
    } else if constexpr (eos == ghl_eos_tabulated) {
      eps = ghl_tabulated_compute_eps_from_rho(ghl_eos_params.get(), rhoB_in);
    }

    *rho_energy_ptr = rhoB_in * (1.0 + eps);
  }

  static inline void dpress_cold_drho__rho(
    const double rho_in,
    double* dpress_cold_drho_ptr) {

    if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      double K, Gamma, eps, P;
      ghl_hybrid_get_K_and_Gamma(ghl_eos_params.get(), rho_in, &K, &Gamma);

      auto press_cold = K * pow(rho_in, Gamma);
      *dpress_cold_drho_ptr = Gamma * press_cold / rho_in;
    } else if constexpr (eos == ghl_eos_tabulated) {
      *dpress_cold_drho_ptr = ghl_tabulated_compute_dP_drho_from_rho(ghl_eos_params.get(), rho_in);
    }
  }

  static inline void dedp__P_cold(
    const double rhoB_in,
    const double rho_energy_in,
    const double P_in,
    double* dedp){

    double dpdrho;
    dpress_cold_drho__rho(rhoB_in, &dpdrho);

    auto const rhoh = rho_energy_in + P_in;

    *dedp = rhoh/(dpdrho*rhoB_in);
  }

  static inline void P_cold_from_rho(
    const double rho_in,
    double *restrict P_cold_ptr) {

    if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      ghl_hybrid_compute_P_cold(ghl_eos_params.get(), rho_in, P_cold_ptr);
    } else if constexpr (eos == ghl_eos_tabulated) {
      *P_cold_ptr = ghl_tabulated_compute_P_from_rho(ghl_eos_params.get(), rho_in);
    }
  }

  // Wrapper to make API consistent without modifying Margherita
  static inline double press_cold_eps_cold__rho(double & eps_cold, double& rhoB_in) {
    double P;
    if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos_params.get(), rhoB_in,
                                              &P, &eps_cold);
    } else if constexpr (eos == ghl_eos_tabulated) {
      P = ghl_tabulated_compute_P_from_rho(ghl_eos_params.get(), rhoB_in);
      eps_cold = ghl_tabulated_compute_eps_from_rho(ghl_eos_params.get(), rhoB_in);
    }
    return P;
  }

  static inline double rho__h_cold(const double h_in) {
    double rho;
    if constexpr (eos == ghl_eos_simple || eos == ghl_eos_hybrid) {
      rho = ghl_hybrid_compute_rho_cold_from_h(ghl_eos_params.get(), h_in);
    } else if constexpr (eos == ghl_eos_tabulated) {
      rho = ghl_tabulated_compute_rho_cold_from_h(ghl_eos_params.get(), h_in);
    }
    return rho;
  }
};

template<class FUKA_eos_t, FUKA_eos_t eos>
std::unique_ptr<ghl_eos_parameters> FUKA_EOS_Wrapper<FUKA_eos_t, eos>::ghl_eos_params = nullptr;
}
}