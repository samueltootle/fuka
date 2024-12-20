
#pragma once
#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#include <cmath>
#include <string>
#include <array>
#include "EOS/standalone/polytrope_file_parser.hpp"

namespace Kadath {
namespace GHL_EOS {

inline auto populate_ghl_polytrope(std::string polytrope_file) {
  using namespace Kadath::FUKA_EOS;
  std::unique_ptr<ghl_eos_parameters> eos_params;
  ghl_eos_parameters eos;

  parse_polytrope_file parser;
  parser(polytrope_file);
  const double rho_atm = 1e-15;
  const double gamma_th = 0.0;
  ghl_initialize_hybrid_eos_functions_and_params(
    rho_atm,
    parser.rhomin,
    parser.rhomax,
    parser.num_pieces,
    parser.rho_tab.data(),
    parser.gamma_tab.data(),
    parser.K0,
    gamma_th,
    &eos
  );
  eos_params = std::make_unique<ghl_eos_parameters>(eos);
  return eos_params;
}

inline auto ghl_setup_polytrope(std::string polytrope_file) {
  auto eos_params = populate_ghl_polytrope(polytrope_file);

  return eos_params;
}
}  // namespace GHL_HELPERS
}  // namespace Kadath
#endif