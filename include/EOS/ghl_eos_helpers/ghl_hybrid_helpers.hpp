
#pragma once
#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#include <array>
#include <cmath>
#include <string>
#include "EOS/standalone/EOS_parfile_parser.hpp"
#include "EOS/standalone/brent.hh"

namespace Kadath {
namespace GHL_EOS {

inline auto populate_ghl_polytrope(std::string polytrope_file) {
  using namespace Kadath::FUKA_EOS;
  std::unique_ptr<ghl_eos_parameters> eos_params;
  ghl_eos_parameters eos;

  parse_polytrope_file parser;
  parser(polytrope_file);
  const double gamma_th = 0.0;
  ghl_initialize_hybrid_eos_functions_and_params(parser.rhomin, parser.rhomin,
                                                 parser.rhomax,
                                                 parser.num_pieces,
                                                 parser.rho_tab.data(),
                                                 parser.gamma_tab.data(),
                                                 parser.K0, gamma_th, &eos);
  eos.press_atm = parser.Pmin;
  eos.press_min = parser.Pmin;
  eos.p_ppoly[0] = parser.Pmin;
  eos_params = std::make_unique<ghl_eos_parameters>(eos);
  return eos_params;
}

inline auto ghl_setup_polytrope(std::string polytrope_file) {
  auto eos_params = populate_ghl_polytrope(polytrope_file);

  return eos_params;
}

inline double ghl_hybrid_rho__h_cold(
  std::unique_ptr<ghl_eos_parameters>& ghl_eos_params,
  double & h_cold) {

  double epsmin, epsmax, press_min, press_max;
  const double rhomin = ghl_eos_params->rho_min;
  const double rhomax = ghl_eos_params->rho_max;

  ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos_params.get(), rhomin, &press_min, &epsmin);
  ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos_params.get(), rhomax, &press_max, &epsmax);

  const auto h_max = 1. + epsmax + press_max/rhomax;
  const auto h_min = 1. + epsmin + press_min/rhomin;

  // Range check
  if (h_cold < h_min) {
    h_cold = h_min;
    return rhomin;
  }
  if (h_cold > h_max) {
    h_cold = h_max;
    return rhomax;
  }

  const auto func = [&](const double &lrho) {
   double eps_cold;
   double rho = std::exp(lrho);
   const double press_cold = press_cold_eps_cold__rho(eps_cold, rho, error);
   return h_cold - ( 1. + eps_cold + press_cold / exp(lrho) );
  };
  auto lrho = zero_brent<>(log(rhomin/10.), 0.999*log(rhomax), 1.0e-13, func);

  double rho = std::exp(lrho);
  return rho;
}
}  // namespace GHL_EOS
}  // namespace Kadath
#endif