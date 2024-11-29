/*
 * Copyright 2022
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
#include <iostream>
#include <string>
#include <utility>
#include "Solvers/bco_solver_utils.hpp"
#include "Solvers/bh_3d_xcts/bh_3d_xcts_solver.hpp"
#include "bco_utilities.hpp"
#include "mpi.h"

/**
 * \addtogroup BBH_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
namespace bco_u = ::Kadath::bco_utils;

template <typename config_t, typename space_t>
bbh_xcts_solver<config_t, space_t>::bbh_xcts_solver(config_t& config_in,
                                                    space_t& space_in,
                                                    Base_tensor& base_in,
                                                    Scalar& conf_in,
                                                    Scalar& lapse_in,
                                                    Vector& shift_in)
    : XCTS_Solver<config_t, space_t>(config_in, space_in, base_in),
      conf(conf_in),
      lapse(lapse_in),
      shift(shift_in),
      fmet(Metric_flat(space_in, base_in)),
      xc1(bco_u::get_center(space_in, space.BH1)),
      xc2(bco_u::get_center(space_in, space.BH2)),
      excluded_doms({space.BH1, space.BH1 + 1, space.BH2, space.BH2 + 1})

{
  // initialize coordinate vector fields
  coord_vectors = default_binary_vector_ary(space);

  update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);
}

// standardized filename for each converged dataset at the end of each stage.
template <typename config_t, typename space_t>
std::string bbh_xcts_solver<config_t, space_t>::converged_filename(
    const std::string stage) const {
  auto res = space.get_domain(0)->get_nbr_points()(0);
  auto M1 = bconfig(MCH, BCO1);
  auto M2 = bconfig(MCH, BCO2);
  if (M2 > M1)
    std::swap(M1, M2);
  bconfig.set(Q) = M2 / M1;
  auto Mtot = M1 + M2;
  std::stringstream ss;
  ss << "BBH";
  if (stage != "")
    ss << "_" << stage << ".";
  else
    ss << ".";
  ss << bconfig(DIST) << "." << bconfig(CHI, BCO1) << "." << bconfig(CHI, BCO2)
     << "." << Mtot << ".q" << bconfig(Q) << "." << std::setfill('0')
     << std::setw(1) << bconfig(NSHELLS, BCO1) << "." << std::setfill('0')
     << std::setw(1) << bconfig(NSHELLS, BCO2) << "." << std::setfill('0')
     << std::setw(2) << res;
  return ss.str();
}

template <typename config_t, typename space_t>
int bbh_xcts_solver<config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;

  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

  // check if we have solved this before til the last stage
  if (solution_exists(last_stage)) {
    if (rank == 0)
      std::cout << "Solved previously: " << bconfig.config_filename()
                << std::endl;
    return EXIT_SUCCESS;
  }

  if (stage_enabled[TOTAL_BC]) {
    this->solver_stage = TOTAL_BC;
    exit_status = solve_stage("TOTAL_BC");

    if (bconfig.control(FIXED_GOMEGA)) {
      bconfig.control(FIXED_GOMEGA) = false;
      exit_status = solve_stage("TOTAL_BC");
    }
  }

  if (stage_enabled[ECC_RED]) {
    this->solver_stage = ECC_RED;
    exit_status = solve_stage("ECC_RED");
  }

  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);

  return exit_status;
}

template <typename config_t, typename space_t>
void bbh_xcts_solver<config_t, space_t>::syst_init(System_of_eqs& syst) {
  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");

  // define numerical constants
  syst.add_cst("4piG", bconfig(QPIG));
  syst.add_cst("PI", M_PI);

  // BH1 fixing parameters
  syst.add_cst("Mm", bconfig(MIRR, BCO1));
  syst.add_cst("chim", bconfig(CHI, BCO1));
  syst.add_cst("CMm", bconfig(MCH, BCO1));

  // BH2 fixing parameters
  syst.add_cst("Mp", bconfig(MIRR, BCO2));
  syst.add_cst("chip", bconfig(CHI, BCO2));
  syst.add_cst("CMp", bconfig(MCH, BCO2));

  // constant vector fields
  syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

  syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
  syst.add_cst("sm", *coord_vectors[S_BCO1]);
  syst.add_cst("mp", *coord_vectors[BCO2_ROT]);
  syst.add_cst("sp", *coord_vectors[S_BCO2]);

  syst.add_cst("ex", *coord_vectors[EX]);
  syst.add_cst("ey", *coord_vectors[EY]);
  syst.add_cst("ez", *coord_vectors[EZ]);
  syst.add_cst("einf", *coord_vectors[S_INF]);

  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P", conf);
  syst.add_var("N", lapse);
  syst.add_var("bet", shift);

  // setup parameters accordingly for a corotating binary
  if (bconfig.control(COROT_BIN)) {
    bconfig.set(OMEGA, BCO1) = 0.;
    bconfig.set(OMEGA, BCO2) = 0.;
  } else {
    // otherwise we setup the system for arbitrarily rotating BHs
    // based on either fixed rotation frequency or variable
    // rotation frequency fixed by a dimensionless
    // spin parameter CHI.
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
      syst.add_cst("omesm", bconfig(FIXED_BCOMEGA, BCO1));
      bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
    } else
      syst.add_var("omesm", bconfig(OMEGA, BCO1));

    syst.add_def(space.BH1 + 2, "ExOme^i = omesm * mm^i");

    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
      syst.add_cst("omesp", bconfig(FIXED_BCOMEGA, BCO2));
      bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
    } else
      syst.add_var("omesp", bconfig(OMEGA, BCO2));
    syst.add_def(space.BH2 + 2, "ExOme^i = omesp * mp^i");
  }

  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");

  // definition describing the binary rotation
  syst.add_def("Morb^i  = mg^i + xaxis * ey^i + yaxis * ex^i");

  syst.add_def(
      "A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
      "Ntilde");

  // ADM linear momenta
  syst.add_def("intPx   = A_ij * ex^j * einf^i / 8 / PI");
  syst.add_def("intPy   = A_ij * ey^j * einf^i / 8 / PI");
  syst.add_def("intPz   = A_ij * ez^j * einf^i / 8 / PI");

  // spin definition on the apparent horizons
  syst.add_def("intSm   = A_ij * mm^i * sm^j   / 8 / PI");
  syst.add_def("intSp   = A_ij * mp^i * sp^j   / 8 / PI");

  // diagnostic from https://arxiv.org/abs/1506.01689
  syst.add_def(ndom - 1, "COMx  = -3 * P^4 * einf^i * ex_i / 8. / PI");
  syst.add_def(ndom - 1, "COMy  = 3 * P^4 * einf^i * ey_i / 8. / PI");
  syst.add_def(ndom - 1, "COMz  = 3 * P^4 * einf^i * ez_i / 8. / PI");

  // BH irreducible mass
  syst.add_def("intMsq  = P^4 / 16. / PI");

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
}

template <typename config_t, typename space_t>
void bbh_xcts_solver<config_t, space_t>::print_diagnostics(
    const System_of_eqs& syst,
    const int ite,
    const double conv) const {
  std::ios_base::fmtflags f(std::cout.flags());
  int ndom = space.get_nbr_domains();

  // local spin angular momentum BH1
  Val_domain integSm(syst.give_val_def("intSm")()(space.BH1 + 2));
  double Sm = space.get_domain(space.BH1 + 2)->integ(integSm, INNER_BC);

  // local spin angular momentum BH2
  Val_domain integSp(syst.give_val_def("intSp")()(space.BH2 + 2));
  double Sp = space.get_domain(space.BH2 + 2)->integ(integSp, INNER_BC);

  // irreducible mass BH1
  Val_domain integMmsq(syst.give_val_def("intMsq")()(space.BH1 + 2));
  double Mirrmsq = space.get_domain(space.BH1 + 2)->integ(integMmsq, INNER_BC);
  double Mirrm = std::sqrt(Mirrmsq);
  // MCH
  double Mchm = std::sqrt(Mirrmsq + Sm * Sm / 4. / Mirrmsq);

  // irreducible mass BH2
  Val_domain integMpsq(syst.give_val_def("intMsq")()(space.BH2 + 2));
  double Mirrpsq = space.get_domain(space.BH2 + 2)->integ(integMpsq, INNER_BC);
  double Mirrp = std::sqrt(Mirrpsq);
  // MCH
  double Mchp = std::sqrt(Mirrpsq + Sp * Sp / 4. / Mirrpsq);

  // ADM Linear Momenta
  Val_domain integPx(syst.give_val_def("intPx")()(ndom - 1));
  double Px = space.get_domain(ndom - 1)->integ(integPx, OUTER_BC);

  Val_domain integPy(syst.give_val_def("intPy")()(ndom - 1));
  double Py = space.get_domain(ndom - 1)->integ(integPy, OUTER_BC);

  Val_domain integPz(syst.give_val_def("intPz")()(ndom - 1));
  double Pz = space.get_domain(ndom - 1)->integ(integPz, OUTER_BC);

  cout << std::string(42, '=') << endl;
  std::cout << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << std::endl;
  std::cout << FORMAT << "Omega: " << bconfig(GOMEGA) << std::endl;
  if (!std::isnan(bconfig.set(ADOT)))
    std::cout << FORMAT << "Adot: " << bconfig(ADOT) << std::endl;
  std::cout << FORMAT << "Axis: " << bconfig(COM) << std::endl;
  std::cout << FORMAT << "P: " << "[" << Px << ", " << Py << ", " << Pz
            << "]\n\n";

  std::cout << FORMAT << "BH1-Mirr: " << Mirrm << std::endl;
  std::cout << FORMAT << "BH1-Mch: " << Mchm << std::endl;
  std::cout << FORMAT << "BH1-R: "
            << bco_u::get_radius(space.get_domain(space.BH1 + 1), EQUI)
            << std::endl;
  std::cout << FORMAT << "BH1-S: " << Sm << std::endl;
  std::cout << FORMAT << "BH1-Chi: " << Sm / Mchm / Mchm << std::endl;
  std::cout << FORMAT << "BH1-Omega: " << bconfig(OMEGA, BCO1) << "\n\n";

  std::cout << FORMAT << "BH2-Mirr: " << Mirrp << std::endl;
  std::cout << FORMAT << "BH2-Mch: " << Mchp << std::endl;
  std::cout << FORMAT << "BH2-R: "
            << bco_u::get_radius(space.get_domain(space.BH2 + 1), EQUI)
            << std::endl;
  std::cout << FORMAT << "BH2-S: " << Sp << std::endl;
  std::cout << FORMAT << "BH2-Chi: " << Sp / Mchp / Mchp << std::endl;
  std::cout << FORMAT << "BH2-Omega: " << bconfig(OMEGA, BCO2) << "\n\n";
  std::cout.flags(f);

  cout << std::string(42, '=') << endl;
}
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath