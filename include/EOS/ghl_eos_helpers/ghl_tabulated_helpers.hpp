#pragma once
#ifdef WITH_GRHAYL_EOS
#include <ghl/ghl.h>
#include <array>
#include <cmath>
#include <string>
#include "EOS/standalone/EOS_parfile_parser.hpp"
#include "EOS/standalone/brent.hh"

namespace Kadath {
namespace GHL_EOS {

/**
 * Update GRHaYL beta equilibrium table to be more useful for elliptic solvers.
 */
inline void ghl_modify_beta_equ_table(
    std::unique_ptr<ghl_eos_parameters>& eos) {
  // This needs more testing regarding its utility
  return;

  // Find the bounded slice of interest
  std::vector<double> lp_of_lr_v;
  std::vector<double> le_of_lr_v;
  std::vector<double> lh_of_lr_v;
  std::vector<double> logrho_v;
  for (int ir = 0; ir < eos->N_rho; ir++) {
    const double rho = exp(eos->table_logrho[ir]);
    const double h = exp(eos->lh_of_lr[ir]);

    if (rho < eos->rho_min || rho > eos->rho_max || h < 1) {
      continue;
    }
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
  eos->lp_of_lr = (double*)malloc(sizeof(double) * eos->N_rho);
  eos->le_of_lr = (double*)malloc(sizeof(double) * eos->N_rho);
  eos->lh_of_lr = (double*)malloc(sizeof(double) * eos->N_rho);
  eos->table_logrho = (double*)malloc(sizeof(double) * eos->N_rho);

  // Copy data
  for (int ir = 0; ir < eos->N_rho; ir++) {
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
  ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(parser.T_beta, &eos);
  eos_params = std::make_unique<ghl_eos_parameters>(eos);

  ghl_modify_beta_equ_table(eos_params);

  return eos_params;
}

inline auto ghl_setup_table(std::string table_par_file) {
  auto eos_params = populate_ghl_tabulated(table_par_file);

  // Immediately free unneeded arrays
  free(eos_params->table_all);
  free(eos_params->table_logT);
  free(eos_params->table_Y_e);
  free(eos_params->table_eps);

  return eos_params;
}

inline double ghl_tabulated_rho__h_cold(
    std::unique_ptr<ghl_eos_parameters>& ghl_eos_params,
    double& h_in) {
  // GRHayL will throw an error rather than enforce table bounds
  // so we enforce them here instead.
  size_t const N = ghl_eos_params->N_rho;
  const double h_min = std::exp(ghl_eos_params->lh_of_lr[0]);
  const double h_max = std::exp(ghl_eos_params->lh_of_lr[N - 1]);
  if (h_in < h_min) {
    h_in = h_min;
    return std::exp(ghl_eos_params->table_logrho[0]);
  } else if (h_in > h_max) {
    h_in = h_max;
    return std::exp(ghl_eos_params->table_logrho[N - 1]);
  }

  // It's only safe to linear interpolate in logp vs logrho space
  // We rootfind the correct specific enthalpy rather than
  // attempting a table look-up using the same method as used
  // in Margherita::Cold_Table
  const auto func = [&](const double& lrho) {
    double const rho = std::exp(lrho);
    double const eps =
        ghl_tabulated_compute_eps_from_rho(ghl_eos_params.get(), rho);
    double const press =
        ghl_tabulated_compute_P_from_rho(ghl_eos_params.get(), rho);

    return h_in - (1. + eps + press / std::exp(lrho));
  };

  // Note the root bounds are quite sensitive given GRHaYL will simply error out
  // rather than enforcing table bounds.
  auto lrho = zero_brent<>(log(ghl_eos_params->rho_min * 1.001),
                           log(ghl_eos_params->rho_max), 1.0e-13, func);
  return std::exp(lrho);
}

// GRHaYL will error out rather than enforce table bounds.
// Therefore we check them here to avoid errors.
inline double ghl_tabulated_check_rho(
    std::unique_ptr<ghl_eos_parameters>& ghl_eos_params,
    const double rho_in) {

  // Enforce lower bound
  double lrho = std::max(std::log(rho_in), ghl_eos_params->table_logrho[0]);

  const auto N = ghl_eos_params->N_rho;
  // Enforce upper bound
  lrho = std::min(lrho, ghl_eos_params->table_logrho[N - 1]);

  return std::exp(lrho);
}

// GRHaYL will error out rather than enforce table bounds.
// Therefore we check them here to avoid errors.
inline double ghl_tabulated_check_press(
    std::unique_ptr<ghl_eos_parameters>& ghl_eos_params,
    const double P_in) {

  // Enforce lower bound
  double lpress = std::max(std::log(P_in), ghl_eos_params->lp_of_lr[0]);

  const auto N = ghl_eos_params->N_rho;
  // Enforce upper bound
  lpress = std::min(lpress, ghl_eos_params->lp_of_lr[N - 1]);

  return std::exp(lpress);
}

}  // namespace GHL_EOS
}  // namespace Kadath
#endif
