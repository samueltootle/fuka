//
// This file is part of Margherita, the light-weight EOS framework
//
//  Copyright (C) 2017, Elias Roland Most
//                      <emost@th.physik.uni-frankfurt.de>
//  Copyright (C) 2019, Ludwig Jens Papenfort
//                      <papenfort@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#define COLDTABLE_SETUP
#include "cold_table.hh"
#include "cold_table_implementation.hh"
#include "cold_pwpoly.hh"
#include "cold_pwpoly_implementation.hh"

#include "./Margherita_EOS.h"
#include "./linear_interp.hh"
#include "./spline.hh"

#include "lorene_io.hh"
#include <fstream>
namespace Kadath {
namespace Margherita {

using Cold_Table_spline_t = Cold_Table_t<0,cubic_spline_t>;

template<bool setup_spline = false>
void setup_Cold_Table(std::string cold_table_name, int cold_lintp_points, double h_cut = 0.0) {
  auto vectors = Lorene_Table(cold_table_name, h_cut);

  assert(vectors[0].front() > 0);  // rho
  assert(vectors[2].front() > 0);  // press

  Cold_Table::rhomin = vectors[0].front();
  Cold_Table::rhomax = vectors[0].back();

  Cold_Table::press_min = vectors[2].front();
  Cold_Table::press_max = vectors[2].back();

  double eps_min = vectors[1].front();
  double eps_max = vectors[1].back();

  double h_min = 1. + eps_min + Cold_Table::press_min / Cold_Table::rhomin;
  double h_max = 1. + eps_max + Cold_Table::press_max / Cold_Table::rhomax;

  Cold_Table::hmin = h_min;
  Cold_Table::hmax = h_max;

  __attribute__((unused)) double K_ext = (h_min - 1) / (2. * Cold_Table::rhomin);

  // Shift eps to ensure positivity
  if (vectors[1].front() < 0) {
    Cold_Table::energy_shift = -2. * vectors[1].front();
    for (auto &v : vectors[1]) v += Cold_Table::energy_shift;
  }

  // Log the table
  for (auto &v : vectors)
    for (auto &w : v) w = log(w);

  // Vector contains rho,eps and press vectors
  // These need to be moved into spline object
  std::unique_ptr<double[]> rho_ptr{new double[vectors[0].size()]};
  std::unique_ptr<double[]> eps_ptr{new double[vectors[1].size()]};
  std::unique_ptr<double[]> press_ptr{new double[vectors[2].size()]};

  for (size_t i = 0; i < vectors[0].size(); ++i) {
    rho_ptr[i] = vectors[0][i];
    eps_ptr[i] = vectors[1][i];
    press_ptr[i] = vectors[2][i];

  }

  // Interpolate table using cubic spline to highly resolved table, suitable for
  // linear  interpolation

  auto spline =
      cubic_spline_t<double, 2>(vectors[0].size(), std::move(rho_ptr),
                                std::move(eps_ptr), std::move(press_ptr));

  std::unique_ptr<double[]> rho_lin_ptr{new double[cold_lintp_points]};
  std::unique_ptr<double[]> eps_lin_ptr{new double[cold_lintp_points]};
  std::unique_ptr<double[]> press_lin_ptr{new double[cold_lintp_points]};

  const auto delta_rho = (log(Cold_Table::rhomax) - log(Cold_Table::rhomin)) /
                         (cold_lintp_points - 1);

  for (int nn = 0; nn < cold_lintp_points; ++nn) {
    const auto rhoL = log(Cold_Table::rhomin) + delta_rho * nn;
    rho_lin_ptr[nn] = rhoL;
    auto res = spline.interpolate<0, 1>(rhoL);
    eps_lin_ptr[nn] = res[0];
    press_lin_ptr[nn] = res[1];
  }

  Cold_Table::lintp = linear_interp_t<double, 2>(
      cold_lintp_points, std::move(rho_lin_ptr), std::move(eps_lin_ptr),
      std::move(press_lin_ptr));

  if(setup_spline){
    Cold_Table_spline_t::rhomin = Cold_Table::rhomin;
    Cold_Table_spline_t::rhomax = Cold_Table::rhomax;

    Cold_Table_spline_t::press_min = Cold_Table::press_min;
    Cold_Table_spline_t::press_max = Cold_Table::press_max;

    Cold_Table_spline_t::energy_shift = Cold_Table::energy_shift;

    Cold_Table_spline_t::lintp = std::move(spline);
  }
}
}
}
