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
  const int max_iter = bconfig.seq_setting(MAX_ITER);
  
  std::string stagename = (fixed) ? "NOROT_FIXED" : "NOROT_BC";

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

  if (fixed) {
    if (rank == 0)
      std::cout << "############################" << std::endl
                << "TOV with a fixed radius" << std::endl
                << "############################" << std::endl;
  } else {
    if (rank == 0) {
      std::cout << "############################" << std::endl
                << "TOV with a resolved surface" << std::endl;                
      std::cout << "############################" << std::endl;
    }
  }

  Scalar one(space);
  one = 1.;
  one.std_base();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  syst_init(syst);

  // in case of "fixed" domain radii
  // syst.add_cst("lev", level); 
  if(fixed) {
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
      syst.add_def(d, "S = delta * 3 * press");
      syst.add_def(d, "Spp = press * delta");
 
      // constraint equations
      syst.add_def(d, "eqnu = delta * ( lap(nu) + scal(grad(nu), grad(nulogA)) ) - 4piG * A^2 * (E + S)") ;
      syst.add_def(d, "eqnulogA = delta * ( lap2(nulogA) + scal(grad(nu), grad(nu)) ) - 2 * 4piG * A^2 * Spp") ;
 
      // definition for the baryonic mass integral
      // syst.add_def(d, "intMb = P^6 * rho");

      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, "firstint = H + log(N)");
      syst.add_def(d, "intMb = rho * A^3 * 4piG / 2");
 
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");
 
      syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(nulogA))") ;
      syst.add_def(d, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))") ;
      break;
    }
  }

  std::string central_fixing_definition{"h - hc"};
  if(seq && !fixed) {
    central_fixing_definition = ::Kadath::FUKA_Syst_tools::set_ns_mass_fixing(syst, bconfig, seq);
  } else {
    syst.add_cst("hc", bconfig(BCO_PARAMS::HC));
  }
 
  // add the constraint equations and demand continuity their normal derivative across domain boundaries
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqnulogA=0", "nulogA", "dn(nulogA)");
  
  // boundary conditions at infinity
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nulogA=0");

  // if the radius of the stellar surface domain is fixed
  // use the helper construction, i.e. a level function with a root defining the radius
  if(fixed){
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
 
  // if surface is resolved, fix the central enthalpy by one of these integrals
  if(!fixed && seq) {
    space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
    space.add_eq_int_inf(syst, "integ(intMadm) = Madm");
  }
 
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;
 
  // solve until convergence is achieved
  while (!endloop) {  
    // do exactly one newton step, given the system above
    endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);
 
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
  update_config_quantities(logh);
  bconfig.set_filename(converged_filename(stagename));
  bconfig.control(CONTROLS::SEQUENCES) = false;
  if (rank == 0) {
    checkpoint();
  }
  return exit_status;
}

/** @}*/
}}