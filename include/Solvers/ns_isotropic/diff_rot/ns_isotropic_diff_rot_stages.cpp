#include "mpi.h"
#include "bco_utilities.hpp"

/**
 * \addtogroup Stages
 * \ingroup NS_XCTS
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_diff_rot_solver<eos_t, config_t, space_t>::differential_rot_stage() {
  int exit_status = EXIT_SUCCESS;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int max_iter = bconfig.seq_setting(MAX_ITER);
  
  // logarithm of the central enthalpy, a variable in the system of equations 
  double loghc = std::log(bconfig(BCO_PARAMS::HC));
  std::string stagename = "DIFF_ROT";

  // We use `config_filename()` vs `config_filename_abs()` since
  // `solution_exists` will probe the HOME_KADATH/COs directory
  // auto const current = bconfig.config_filename();
  // if(!bconfig.control(RESOLVE) && solution_exists(stagename)) {    
  //   if(rank == 0)
  //     std::cout << "Solved previously: " \
  //               << bconfig.config_filename_abs() << std::endl;
  //   return (current == bconfig.config_filename()) ? \
  //     EXIT_SUCCESS : RELOAD_FILE;
  // }

  if (rank == 0)
    std::cout << "############################" << std::endl
              << "Differentially Rotating NS Solver" << std::endl
              << "Omega_c: " << bconfig(BCO_PARAMS::OMEGA) <<std::endl
              << "############################" << std::endl;

  // Sad tool to make system of equations work with constants
  Scalar one(space);
  one = 1;
  one.std_base();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  syst_init(syst);

  std::string central_fixing_definition{"h - hc"};
  if(seq) {
    central_fixing_definition = ::Kadath::FUKA_Syst_tools::set_ns_mass_fixing(syst, bconfig, seq);
  } else {
    syst.add_cst("hc" , bconfig(BCO_PARAMS::HC));
    syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
    syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
  }

  syst.add_cst("one", one);

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:

      // sources
      syst.add_def(d, "E = Wsq * press * h - press * delta");
      syst.add_def(d, "Srrtt = press * delta");
      syst.add_def(d, "pphi = B * (E + Srrtt) * U");
      syst.add_def(d, "Spp = delta * press * (1 + Usq) + E * Usq");
      syst.add_def(d, "S = 2 * Srrtt + Spp");
 
      // constraint equations
      syst.add_def(d, "eqnu  = delta * lap(nu) + delta * scal(grad(nu), grad(nu + log(B))) "
                            "- delta * multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w)) "
                            "- 4piG * A^2 * (E + S)");
      syst.add_def(d, "eqAterm = delta * lap2(lapAterm) + delta * scal(grad(nu), grad(nu))"
                      "- 3 * delta * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))"
                      "- 2 * 4piG * A^2 * Spp");
      syst.add_def(d, "eqBterm = delta * lap2(lapBterm) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
      syst.add_def(d, "eqwrsint = delta * lap(wrsint) - delta * multrsint(scal(grad(w), grad(nu - 3 * log(B))))"
                          "+ 4 * 4piG * N * A^2 / B * pphi");
 
      // definition for the baryonic mass integral
      syst.add_def(d, "intMb = W * rho * A^2 * B * 4piG / 2");

      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      // This MUST be populated by the differential rotation law, ie KEH_Law(syst)
      // syst.add_def(d, "firstint = H + log(N) - 0.5 * log(Wsq)");
 
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                      "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))"
                "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqBterm = lap2(lapBterm)");
      syst.add_def(d, "eqwrsint = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
      break;
    }
  }
  KEH_law(syst);
 
  // add the constraint equations and demand continuity their normal derivative across domain boundaries
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqAterm=0", "lapAterm", "dn(lapAterm)");
  space.add_eq(syst, "eqBterm=0", "lapBterm", "dn(lapBterm)");
  space.add_eq(syst, "eqwrsint=0", "wrsint", "dn(wrsint)");
  
  // boundary conditions at infinity
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapAterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapBterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "wrsint=0");

  // Fix surface based on vanishing log specific enthalpy
  syst.add_eq_bc(1, OUTER_BC, "H = 0");
 
  // first integral in the innermost domains with non-zero matter content
  // and condition on the central value, either fixed directly or by the
  // integral below
  syst.add_eq_first_integral(0, 1, "firstint", central_fixing_definition.c_str());
 
  // if surface is resolved, fix the central enthalpy by one of these integrals
  space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
  space.add_eq_int_inf(syst, "integ(intMadm) = Madm");
 
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;
    syst.sec_member();
 
  // solve until convergence is achieved
  while (!endloop) {  
    // do exactly one newton step, given the system above
    endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);
 
    update_config_quantities(logh);
    // output files at this iteration and print diagnostics
    std::stringstream ss;
    ss << "diff_rot_ckpt_";
    ss << "_" << ite - 1;
    bconfig.set_filename(ss.str());
    if (rank == 0) {
      print_diagnostics(syst, ite, conv);
      std::cout << std::endl;
      if(bconfig.control(CHECKPOINT))
        checkpoint();
    }

    ite++;
    check_max_iter_exceeded(rank, ite, conv);
  }
  
  update_config_quantities(logh);
  bconfig.set_filename(converged_filename(stagename));
  if (rank == 0) {
    checkpoint();
  }
  return exit_status;
}

/** @}*/
}}