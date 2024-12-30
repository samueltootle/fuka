
#pragma once
#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#include <array>
#include <cmath>
#include <string>
#include "EOS/standalone/EOS_parfile_parser.hpp"

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

inline auto populate_ghl_tabulated(std::string table_par_file) {
  using namespace Kadath::FUKA_EOS;
  std::unique_ptr<ghl_eos_parameters> eos_params;
  ghl_eos_parameters eos;

  parse_3d_EOS_par_file parser;
  parser(table_par_file);
  ghl_initialize_tabulated_eos_functions_and_params(
      parser.table_abspath.c_str(), parser.rho_atm, parser.rho_min,
      parser.rho_max, parser.Ye_atm, parser.Ye_min, parser.Ye_max, parser.T_atm,
      parser.T_min, parser.T_max, &eos);
  ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(1.1e-2, &eos);
  eos_params = std::make_unique<ghl_eos_parameters>(eos);
  return eos_params;
}

inline auto ghl_setup_table(std::string table_par_file) {
  auto eos_params = populate_ghl_tabulated(table_par_file);

  return eos_params;
}
}  // namespace GHL_EOS
}  // namespace Kadath
#endif