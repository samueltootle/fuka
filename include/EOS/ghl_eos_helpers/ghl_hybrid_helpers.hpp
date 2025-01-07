
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


/**
 * Update GRHaYL beta equilibrium table to be more useful for elliptic solvers.
 */
inline void ghl_modify_beta_equ_table(std::unique_ptr<ghl_eos_parameters>& eos) {
  // Find the bounded slice of interest
  std::vector<double> lp_of_lr_v;
  std::vector<double> le_of_lr_v;
  std::vector<double> lh_of_lr_v;
  std::vector<double> logrho_v;
  for(int ir = 0; ir < eos->N_rho; ir++) {
    const double rho = exp(eos->table_logrho[ir]);
    const double h = exp(eos->lh_of_lr[ir]);
    // if(rho < eos->rho_min || rho > eos->rho_max || h < 1) {
    //   continue;
    // }
    lp_of_lr_v.push_back(eos->lp_of_lr[ir]);
    le_of_lr_v.push_back(eos->le_of_lr[ir]);
    lh_of_lr_v.push_back(eos->lh_of_lr[ir]);
    logrho_v.push_back(eos->table_logrho[ir]);
  }

  // Update 1D table slice size
  eos->N_rho = lp_of_lr_v.size();

  // Free old arrays
  free(eos->lp_of_lr);
  free(eos->le_of_lr);
  free(eos->lh_of_lr);
  free(eos->table_logrho);

  // Allocate new arrays
  eos->lp_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  eos->le_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  eos->lh_of_lr = (double *)malloc(sizeof(double) * eos->N_rho);
  eos->table_logrho = (double *)malloc(sizeof(double) * eos->N_rho);

  // Copy data
  for(int ir = 0; ir < eos->N_rho; ir++) {
    eos->lp_of_lr[ir] = lp_of_lr_v[ir];
    eos->le_of_lr[ir] = le_of_lr_v[ir];
    eos->lh_of_lr[ir] = lh_of_lr_v[ir];
    eos->table_logrho[ir] = logrho_v[ir];
  }
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

  ghl_modify_beta_equ_table(eos_params);

  // for(int ir = 0; ir < eos_params->N_rho; ir++) {
  //   std::cout << std::exp(eos_params->table_logrho[ir]) << " "
  //             << std::exp(eos_params->lh_of_lr[ir]) << " "
  //             << std::exp(eos_params->lp_of_lr[ir]) << " "
  //             << std::exp(eos_params->le_of_lr[ir]) << std::endl;
  // }
  // std::cout << "##################################################################################\n";
  // for(int ir = 0; ir < eos_params->N_rho; ir++) {
  //   std::cout << eos_params->table_logrho[ir] << " "
  //             << eos_params->lh_of_lr[ir] << " "
  //             << eos_params->lp_of_lr[ir] << " "
  //             << eos_params->le_of_lr[ir] << std::endl;
  // }
  return eos_params;
}

inline auto ghl_setup_table(std::string table_par_file) {
  auto eos_params = populate_ghl_tabulated(table_par_file);

  return eos_params;
}
}  // namespace GHL_EOS
}  // namespace Kadath
#endif