#include "mpi.h"
#include "bco_utilities.hpp"
#include "Solvers/fuka_syst/fuka_syst_setup.hpp"

/**
 * \addtogroup Stages
 * \ingroup NS_XCTS
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_norot_solver<eos_t, config_t, space_t>::norot_stage(bool fixed) {
  int exit_status = EXIT_SUCCESS;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::string stagename = (fixed) ? "NOROT_FIXED" : "NOROT_BC";

  // We use `config_filename()` vs `config_filename_abs()` since
  // `solution_exists` will probe the HOME_KADATH/COs directory
  /*auto const current = bconfig.config_filename();
  if(!bconfig.control(RESOLVE) && solution_exists(stagename)) {    
    if(rank == 0)
      std::cout << "Solved previously: " \
                << bconfig.config_filename_abs() << std::endl;
    return (current == bconfig.config_filename()) ? \
      EXIT_SUCCESS : RELOAD_FILE;
  }*/

  Scalar one(space);
  one = 1.;
  one.std_base();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  syst_init(syst);

  // in case of "fixed" domain radii
  // syst.add_cst("lev", level); 
  if(fixed && bconfig.control(CONTROLS::USE_FIXED_R)) {
  syst.add_cst("one", one);
  syst.add_cst("fixedR", bconfig(BCO_PARAMS::RMID));
  syst.add_def("lev = multr(multr(one)) - fixedR * fixedR");
  }
 
  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      // sources
      syst.add_def(d, "E = press * h - press * delta");
      syst.add_def(d, "Srrtt = press * delta");
      syst.add_def(d, "Spp = press * delta");
      syst.add_def(d, "S = 2 * Srrtt + Spp");
 
      // constraint equations
      syst.add_def(d, "eqnu  = delta * lap(nu) + delta * scal(grad(nu), grad(nu + log(B))) "
                            "- 4piG * A^2 * (E + S)");
      syst.add_def(d, "eqBterm = delta * lap2(lapBterm) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
      syst.add_def(d, "eqAterm = delta * lap2(lapAterm) + delta * scal(grad(nu), grad(nu))"
                              "- 2 * 4piG * A^2 * Spp");
 
      // definition for the baryonic mass integral
      syst.add_def(d, "intMb = rho * A^2 * B * 4piG / 2");

      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, "firstint = H + log(N)");
 
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");
 
      syst.add_def(d, "eqnu    = lap(nu) + scal(grad(nu), grad(nu + log(B)))");
      syst.add_def(d, "eqBterm = lap2(lapBterm)");
      syst.add_def(d, "eqAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))");
      break;
    }
  }

  std::string central_fixing_definition{"h - hc"};
  std::string output_str{};
  if(seq && !fixed) {
    central_fixing_definition = ::Kadath::FUKA_Syst_tools::set_ns_mass_fixing(syst, bconfig, seq);
    output_str = ::Kadath::FUKA_Syst_tools::get_ns_mass_fixing_output(bconfig, seq);
  } else {
    syst.add_cst("hc", bconfig(BCO_PARAMS::HC));
    std::stringstream output;
    output << "Mass fixed using central enthalpy (hc) = " << bconfig(BCO_PARAMS::HC);
    output_str = output.str();
  }

  if (fixed  && bconfig.control(CONTROLS::USE_FIXED_R)) {
    if (rank == 0)
      std::cout << "############################" << std::endl
                << "TOV with a fixed radius" << std::endl
                << output_str << std::endl
                << "############################" << std::endl;
  } else {
    if (rank == 0) {
      std::cout << "############################" << std::endl
                << "TOV with a resolved surface" << std::endl
                << output_str << std::endl
                << "############################" << std::endl;
    }
  }
 
  // add the constraint equations and demand continuity their normal derivative across domain boundaries
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqAterm=0", "lapAterm", "dn(lapAterm)");
  space.add_eq(syst, "eqBterm=0", "lapBterm", "dn(lapBterm)");
  
  // boundary conditions at infinity
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapAterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapBterm=0");

  // if the radius of the stellar surface domain is fixed
  // use the helper construction, i.e. a level function with a root defining the radius
  if(fixed && bconfig.control(CONTROLS::USE_FIXED_R)){
    syst.add_eq_bc(1, OUTER_BC, "lev = 0");
  }
  // if the surface is resolved, define it to be where the matter vanishes
  else{
    syst.add_eq_bc(1, OUTER_BC, "H = 0");
  }
 
  // first integral in the innermost domains with non-zero matter content
  // and condition on the central value, either fixed directly or by the
  // integral below
  syst.add_eq_first_integral(0, 1, "firstint", central_fixing_definition.c_str());
 
  // if surface is resolved, add relevant equations
  if(!fixed && seq) {
    auto idx{seq->mass_idx()};
    switch(idx) {
      case BCO_PARAMS::MADM:
        space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
        space.add_eq_int_inf(syst, "integ(intMadm) = Madm");
        break;
      case BCO_PARAMS::MB:
        syst.add_var("hc", bconfig(BCO_PARAMS::HC));
        syst.add_cst("Mb"  , bconfig(BCO_PARAMS::MB));
        space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
        break;
      default:
        break;
    }
  }
 
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;
 
  // solve until convergence is achieved
  while (!endloop) {  
    // do exactly one newton step, given the system above
    endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);
    if(!fixed)
      update_config_quantities(syst);
 
    // output files at this iteration and print diagnostics
    std::stringstream ss;
    ss << "norot_2d_";
    if(fixed) {
      ss << "fixed";
    }else {
      ss << "norot_bc";
    }
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
  if(!fixed)
    update_config_quantities(syst);
  bconfig.set_filename(converged_filename(stagename));
  bconfig.control(CONTROLS::SEQUENCES) = false;
  if (rank == 0) {
    checkpoint();
  }
  return exit_status;
}
/** @}*/
}}