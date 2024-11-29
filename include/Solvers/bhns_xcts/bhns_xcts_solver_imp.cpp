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
#include "Solvers/bco_solver_utils.hpp"
#include "Solvers/bh_3d_xcts/bh_3d_xcts_solver.hpp"
#include "Solvers/ns_3d_xcts/ns_3d_xcts_solver.hpp"
#include "bco_utilities.hpp"
#include "kadath.hpp"
#include "mpi.h"

/**
 * \addtogroup BHNS_XCTS
 * \ingroup FUKA
 * @{*/

namespace Kadath {
namespace FUKA_Solvers {
using namespace ::Kadath::Margherita;
using namespace ::Kadath::bco_utils;

template <class eos_t, typename config_t, typename space_t>
bhns_xcts_solver<eos_t, config_t, space_t>::bhns_xcts_solver(
    config_t& config_in,
    space_t& space_in,
    Base_tensor& base_in,
    Scalar& conf_in,
    Scalar& lapse_in,
    Vector& shift_in,
    Scalar& logh_in,
    Scalar& phi_in)
    : XCTS_Solver<config_t, space_t>(config_in, space_in, base_in),
      conf(conf_in),
      lapse(lapse_in),
      shift(shift_in),
      logh(logh_in),
      phi(phi_in),
      fmet(Metric_flat(space_in, base_in)),
      xc1(get_center(space_in, space.NS)),
      xc2(get_center(space_in, space.BH)),
      xo(get_center(space, ndom - 1)),
      excluded_doms({space.BH, space.BH + 1}) {
  // initialize coordinate vector fields
  coord_vectors = default_binary_vector_ary(space);

  update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);
}

// standardized filename for each converged dataset at the end of each stage.
template <class eos_t, typename config_t, typename space_t>
std::string bhns_xcts_solver<eos_t, config_t, space_t>::converged_filename(
    const std::string stage) const {
  const std::string eosname{extract_eos_name(BCO1)};
  std::stringstream ss;
  auto res = space.get_domain(0)->get_nbr_points()(0);
  auto M1 = bconfig(MADM, BCO1);
  auto M2 = bconfig(MCH, BCO2);
  if (M2 > M1)
    std::swap(M1, M2);
  bconfig.set(Q) = M2 / M1;
  auto Mtot = M1 + M2;
  ss << "BHNS";
  if (stage != "")
    ss << "_" << stage << ".";
  else
    ss << ".";
  ss << eosname << "." << bconfig(DIST) << "." << bconfig(CHI, BCO1) << "."
     << bconfig(CHI, BCO2) << "." << Mtot << ".q" << bconfig(Q) << "."
     << std::setfill('0') << std::setw(1) << bconfig(NSHELLS, BCO1) << "."
     << std::setfill('0') << std::setw(1) << bconfig(NSHELLS, BCO2) << "."
     << std::setfill('0') << std::setw(2) << res;
  return ss.str();
}

template <class eos_t, typename config_t, typename space_t>
int bhns_xcts_solver<eos_t, config_t, space_t>::solve() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int exit_status = EXIT_SUCCESS;

  auto stage_enabled = bconfig.return_stages();
  auto [last_stage, last_stage_idx] = get_last_enabled(MSTAGE, stage_enabled);

  // check if we have solved this before
  if (solution_exists(last_stage)) {
    if (rank == 0)
      std::cout << "Solved previously: " << bconfig.config_filename_abs()
                << std::endl;
    return EXIT_SUCCESS;
  }

  // overview output
  if (rank == 0) {
    std::cout << "=================================" << endl;
    std::cout << "BHNS grav input" << endl;
    std::cout << "Distance: " << bconfig(DIST) << endl;
    std::cout << "Omega guess: " << bconfig(GOMEGA) << endl;
    std::cout << "Units: " << bconfig(QPIG) << endl;
    std::cout << "=================================" << endl;
  }

  if (stage_enabled[TOTAL]) {
    this->solver_stage = TOTAL;
    exit_status = hydrostatic_equilibrium_stage("TOTAL");
    if (exit_status != EXIT_SUCCESS)
      std::_Exit(EXIT_FAILURE);
  }

  if (stage_enabled[TOTAL_BC]) {
    this->solver_stage = TOTAL_BC;
    if (bconfig.control(FIXED_GOMEGA)) {
      exit_status = hydro_rescaling_stages("TOTAL_BC_FIXED_OMEGA");

      // In the event we are creating a new BHNS sequence, FIXED_OMEGA
      // is used only for the initial solution before resolving the hydro
      // consistently.  The constrols are then disabled since we only need to
      // do this once
      if (exit_status == EXIT_SUCCESS && bconfig.control(SEQUENCES)) {
        bconfig.control(SEQUENCES) = false;
        bconfig.control(FIXED_GOMEGA) = false;

        exit_status = hydrostatic_equilibrium_stage("TOTAL_BC");
      }
    } else
      exit_status = hydrostatic_equilibrium_stage("TOTAL_BC");

    if (exit_status != EXIT_SUCCESS)
      std::_Exit(EXIT_FAILURE);
  }
  if (stage_enabled[ECC_RED]) {
    this->solver_stage = ECC_RED;
    exit_status = hydro_rescaling_stages("ECC_RED");
    if (exit_status != EXIT_SUCCESS)
      std::_Exit(EXIT_FAILURE);
  }

  // Barrier needed in case we need to read from the previous output
  MPI_Barrier(MPI_COMM_WORLD);

  return exit_status;
}

template <class eos_t, typename config_t, typename space_t>
void bhns_xcts_solver<eos_t, config_t, space_t>::syst_init(
    System_of_eqs& syst) {
  using namespace ::Kadath::Margherita;

  const int ndom = space.get_nbr_domains();
  // call the (flat) conformal metric "f"
  fmet.set_system(syst, "f");

  // define the EOS operators
  Param p;
  syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
  syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

  // define numerical constants
  syst.add_cst("4piG", bconfig(QPIG));

  // baryonic mass and dimensionless spin are fixed input parameters
  // along with the ADM of each NS at infinite separation
  syst.add_cst("Madm1", bconfig(MADM, BCO1));
  syst.add_cst("Mb1", bconfig(MB, BCO1));
  syst.add_cst("chi1", bconfig(CHI, BCO1));

  syst.add_cst("Mch", bconfig(MCH, BCO2));
  syst.add_cst("Mirr", bconfig(MIRR, BCO2));
  syst.add_cst("chi2", bconfig(CHI, BCO2));

  // coordinate dependent (rotational) vector fields
  syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);
  syst.add_cst("mm", *coord_vectors[BCO1_ROT]);
  syst.add_cst("mp", *coord_vectors[BCO2_ROT]);

  // cartesianbasis vector fields
  syst.add_cst("ex", *coord_vectors[EX]);
  syst.add_cst("ey", *coord_vectors[EY]);
  syst.add_cst("ez", *coord_vectors[EZ]);

  // flat spac surface normals
  syst.add_cst("sm", *coord_vectors[S_BCO1]);
  syst.add_cst("sp", *coord_vectors[S_BCO2]);
  syst.add_cst("einf", *coord_vectors[S_INF]);

  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("P", conf);
  syst.add_var("N", lapse);
  syst.add_var("bet", shift);

  // quasi-local measurement of the component ADM masses
  syst.add_var("qlMadm1", bconfig(QLMADM, BCO1));

  // check if corotation is considered and adjust
  // rotational velocity components accordingly
  if (bconfig.control(COROT_BIN)) {
    // for a corotating binary, there is no angular frequency parameter
    // since stellar rotation is fixed by the orbital frequency
    bconfig.set(OMEGA, BCO1) = 0.;
    syst.add_cst("omes1", bconfig(OMEGA, BCO1));

    bconfig.set(OMEGA, BCO2) = 0.;
    syst.add_cst("omes2", bconfig(OMEGA, BCO2));
  } else {
    // in the arbitrary spinning / irrotational case
    // the irrotational part is defined by the velocity potential
    syst.add_var("phi", phi);

    // check if stellar omega is fixed or has to be solved
    // for a specific dimensionless spin
    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
      syst.add_cst("omes1", bconfig(FIXED_BCOMEGA, BCO1));
      bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA, BCO1);
    } else {
      syst.add_var("omes1", bconfig(OMEGA, BCO1));
    }

    if (!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
      syst.add_cst("omes2", bconfig(FIXED_BCOMEGA, BCO2));
      bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA, BCO2);
    } else {
      syst.add_var("omes2", bconfig(OMEGA, BCO2));
    }

    // for mixed spins, we have additional definitions
    // that are required for both objects that describe
    // the local rotation field
    for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
      syst.add_def(d, "s^i  = omes1 * mm^i");
      syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
    }
    syst.add_def(space.ADAPTEDBH + 1, "s^i = omes2 * mp^i");
  }

  // define common combinations of conformal factor and lapse
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
  syst.add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
  syst.add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
  syst.add_def("intMsq  = P^4 / 4. / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable
  // in this system
  syst.add_def("h = exp(H)");

  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");
  syst.add_def("dHdlnrho = dHdlnrho(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

  // the conformal extrinsic curvature
  syst.add_def(
      "A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
      "Ntilde");

  // ADM linear momentum, surface integrant at infinity
  syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i");
  syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i");
  syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i");

  // quasi-local spin, surface integral outside the stellar matter distribution
  syst.add_def(space.ADAPTEDNS + 1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
  syst.add_def(space.ADAPTEDBH + 1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");
}

// runtime diagnostics specific for rotating solutions
template <class eos_t, typename config_t, typename space_t>
void bhns_xcts_solver<eos_t, config_t, space_t>::print_diagnostics(
    System_of_eqs const& syst,
    const int ite,
    const double conv) const {
  std::ios_base::fmtflags f(std::cout.flags());
  using namespace ::Kadath::bco_utils;

  double baryonic_massNS = 0.;
  double ql_massNS = 0.;
  for (int d = space.NS; d <= space.ADAPTEDNS; ++d) {
    baryonic_massNS += syst.give_val_def("intMb")()(d).integ_volume();
    ql_massNS += syst.give_val_def("intM")()(d).integ_volume();
  }

  auto [rmin, rmax] = get_rmin_rmax(space, space.ADAPTEDNS);
  double r_bh = get_radius(space.get_domain(space.ADAPTEDBH), EQUI);

  double mirrsq =
      space.get_domain(space.BH + 2)
          ->integ(syst.give_val_def("intMsq")()(space.BH + 2), INNER_BC);
  double Mirr = std::sqrt(mirrsq);

  Val_domain integS1(syst.give_val_def("intS1")()(space.ADAPTEDNS + 1));
  double S1 = space.get_domain(space.ADAPTEDNS + 1)->integ(integS1, OUTER_BC);
  double chi1 = S1 / bconfig(MADM, BCO1) / bconfig(MADM, BCO1);
  double chiql1 = S1 / ql_massNS / ql_massNS;

  Val_domain integS2(syst.give_val_def("intS2")()(space.ADAPTEDBH + 1));
  double S2 = space.get_domain(space.ADAPTEDBH + 1)->integ(integS2, OUTER_BC);
  double Mch = std::sqrt(mirrsq + S2 * S2 / 4. / mirrsq);
  double chi2 = S2 / Mch / Mch;

  Val_domain integPx(syst.give_val_def("intPx")()(ndom - 1));
  double Px = space.get_domain(ndom - 1)->integ(integPx, OUTER_BC);

  Val_domain integPy(syst.give_val_def("intPy")()(ndom - 1));
  double Py = space.get_domain(ndom - 1)->integ(integPy, OUTER_BC);

  Val_domain integPz(syst.give_val_def("intPz")()(ndom - 1));
  double Pz = space.get_domain(ndom - 1)->integ(integPz, OUTER_BC);

  std::cout << "=======================================" << endl;
  std::cout << FORMAT << "Iter: " << ite << std::endl
            << FORMAT << "Error: " << conv << "\n";
  std::cout << FORMAT << "Omega: " << bconfig(GOMEGA) << endl;
  std::cout << FORMAT << "Axis: " << bconfig(COM) << endl;
  std::cout << FORMAT << "Padm: " << "[" << Px << ", " << Py << ", " << Pz
            << "]\n\n";

  std::cout << FORMAT << "NS-Mb: " << baryonic_massNS << "["
            << bconfig(MB, BCO1) << "]\n";
  std::cout << FORMAT << "NS-Madm_ql: " << ql_massNS << "["
            << bconfig(MADM, BCO1) << "]\n";
  std::cout << FORMAT << "NS-R: " << rmin << " " << rmax << std::endl
            << FORMAT << "NS-S: " << S1 << std::endl;
  std::cout << FORMAT << "NS-Chi: " << chi1 << "[" << bconfig(CHI, BCO1)
            << "]\n";
  std::cout << FORMAT << "NS-qlChi: " << chiql1 << std::endl;
  std::cout << FORMAT << "NS-Omega: " << bconfig(OMEGA, BCO1) << "\n\n";

  std::cout << FORMAT << "BH-Mirr: " << Mirr << "[" << bconfig(MIRR, BCO2)
            << "]\n";
  std::cout << FORMAT << "BH-Mch: " << Mch << "[" << bconfig(MCH, BCO2)
            << "]\n";
  std::cout << FORMAT << "BH-R: " << r_bh << std::endl;
  std::cout << FORMAT << "BH-S: " << S2 << std::endl;
  std::cout << FORMAT << "BH-Chi: " << chi2 << "[" << bconfig(CHI, BCO2)
            << "]\n";
  std::cout << FORMAT << "BH-Omega: " << bconfig(OMEGA, BCO2) << std::endl;
  std::cout.flags(f);
  std::cout << "=======================================" << endl;
}  // end print_diagnostics
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath