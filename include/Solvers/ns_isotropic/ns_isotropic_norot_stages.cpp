#include "mpi.h"
#include "bco_utilities.hpp"

/**
 * \addtogroup Stages
 * \ingroup NS_XCTS
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {

template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_norot_solver<eos_t, config_t, space_t>::norot_stage (bool fixed) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence threshold
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy, a variable in the system of equations 
  double loghc = std::log(bconfig(BCO_PARAMS::HC));
  std::string stagename = (fixed) ? "NOROT_FIXED" : "NOROT";

  // The 2D solver is very sensitive so a 'fixed' stage is always needed
  // in order to stabalize convergence based on a fixed radius
  Scalar level(space);
  for (int d = 0; d < ndom - 1; d++)
    level.set_domain(d) = space.get_domain(d)->get_radius() * space.get_domain(d)->get_radius() 
                         -bconfig(RMID) * bconfig(RMID);
  level.set_domain(ndom - 1) = 1;
  level.std_base();
  
  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);

  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  
  syst.add_cst("Hc", loghc);  

  // the basic fields, conformal factor, lapse and (log) enthalpy
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("nulogA", nulogA);

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(nulogA - nu)");
  syst.add_cst("lev", level); 

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  syst.add_def(ndom - 1, "intMadm = -dr(A) / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      // sources
      syst.add_def(d, "E = press * (1 + eps)");
      syst.add_def(d, "S = delta * 3 * press");
      syst.add_def(d, "Spp = press * delta");
      
      // constraint equations
      syst.add_def(d, "eqnu = delta * ( lap(nu) + scal(grad(nu), grad(nulogA)) ) - 4piG * A^2 * (E + S)") ;
      syst.add_def(d, "eqnulogA = delta * ( lap2(nulogA) + scal(grad(nu), grad(nu)) ) - 2 * 4piG * A^2 * Spp") ;
      // Extra...
      // syst.add_def(d, "eqNA = dr(drNA) + 3 * divr(drNA) - 4 * 4piG * NA * A^2 * press") ;


      // // definition for the baryonic mass integral
      // syst.add_def(d, "intMb = P^6 * rho");
      
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, "firstint = H + log(N)");

      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(nulogA))") ;
      syst.add_def(d, "eqnulogA = lap2(nulogA) + scal(grad(nu), grad(nu))") ;
      break;
    }
  }
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqnulogA=0", "nulogA", "dn(nulogA)");
  
  
  if(fixed)
  syst.add_eq_bc(1, OUTER_BC, "lev=0");
  else
  syst.add_eq_bc(1, OUTER_BC, "H=0");
  
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");

  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nulogA=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // solve until convergence is achieved
  while (!endloop) {  
    // do exactly one newton step, given the system above
    endloop = syst.do_newton(conv_thres, conv);
 
    update_config_quantities(loghc);
    // output files at this iteration and print diagnostics
    std::stringstream ss;
    ss << "norot_3d_";
    if(fixed) {
      ss << "fixed";
    }else {
      ss << "norot_bc";
    }
    ss << "_" << ite - 1;
    bconfig.set(QLMADM) = bconfig(MADM) ;
    bconfig.set_filename(ss.str());
    if (rank == 0) {
      print_diagnostics_norot(syst, ite, conv);
      std::cout << std::endl;
      if(bconfig.control(CHECKPOINT))
        checkpoint();
    }

    ite++;
    check_max_iter_exceeded(rank, ite, conv);
  }
 
  bconfig.set_filename(converged_filename(stagename));
  if (rank == 0) {
    checkpoint();
  }
  bconfig.set(BCO_PARAMS::MADM) = 
    space.get_domain(ndom-1)->integ(syst.give_val_def("intMadm")()(ndom-1), OUTER_BC);
  
  update_config(bconfig, logh);
  bconfig.set_filename(converged_filename(stage_name, bconfig));
  bconfig.control(CONTROLS::SEQUENCES) = false;
  Scalar bet(space);  
    {  
    Scalar B(exp(nulogA - nu));
    Scalar N(exp(nu));
    Scalar tmp(N * B - 1);
    bet = Scalar(tmp.mult_sin_theta().mult_r());
    }
  Scalar wrsint(space);
  wrsint.annule_hard();
  wrsint.std_base();
  if(rank == 0)
    bco_utils::save_to_file(space, bconfig, nulogA, nu, logh, bet, wrsint);
  MPI_Barrier(MPI_COMM_WORLD);
  return EXIT_SUCCESS;
}


template<class eos_t, typename config_t, typename space_t>
int ns_isotropic_norot_solver<eos_t, config_t, space_t>::norot_stage(bool fixed) {
  int exit_status = EXIT_SUCCESS;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int max_iter = bconfig.seq_setting(MAX_ITER);
  double loghc = std::log(bconfig(HC));
  std::string stagename = (fixed) ? "NOROT_FIXED" : "NOROT_BC";

  // We use `config_filename()` vs `config_filename_abs()` since
  // `solution_exists` will probe the HOME_KADATH/COs directory
  auto const current = bconfig.config_filename();
  if(!bconfig.control(RESOLVE) && solution_exists(stagename)) {    
    if(rank == 0)
      std::cout << "Solved previously: " \
                << bconfig.config_filename_abs() << std::endl;
    return (current == bconfig.config_filename()) ? \
      EXIT_SUCCESS : RELOAD_FILE;
  }

  if (fixed) {
    if (rank == 0)
      std::cout << "############################" << std::endl
                << "TOV with a fixed radius" << std::endl
                << "############################" << std::endl;
  } else {
    if (rank == 0) {
      std::cout << "############################" << std::endl
                << "TOV with a resolved surface" << std::endl;
      if(bconfig.control(MB_FIXING))
        std::cout << "with Baryonic Mass fixing\n";
      else
        std::cout << "with ADM Mass fixing\n";
                
      std::cout << "############################" << std::endl;
    }
  }
  
  // set up radius and leve field in case of "fixed"
  scalar_ary_t coord_scalars;
  coord_scalars[R_BCO1] = Scalar(space);

  update_fields_co(cfields, coord_vectors, coord_scalars, 0.);
  
  // a level function, defining a root at a given fixed radius
  // helper construction to force the system to attain a fixed radius
  // instead resolving the correct surface
  Scalar level(space);
  level = (*coord_scalars[R_BCO1]) * (*coord_scalars[R_BCO1]) -  bconfig(RMID) * bconfig(RMID);
  level.std_base();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  syst.add_var("H"   , logh);
  syst_init(syst);

  if(bconfig.control(MB_FIXING)) {
    syst.add_cst("Mb"  , bconfig(MB));
    syst.add_var("Madm", bconfig(MADM));
  }
  else {
    syst.add_var("Mb"  , bconfig(MB));
    syst.add_cst("Madm", bconfig(MADM));
  }

  // in case of a fixed radius solve the TOV with the given fixed central enthalpy
  // in case of a resolved surface, solve for the central enthalpy
  if(fixed){
    syst.add_cst("Hc", loghc);
  }else {
    syst.add_var("Hc", loghc);
  }
 
  // in case of "fixed" domain radii
  syst.add_cst("lev" , level);
 
  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      // sources
      syst.add_def(d, "Etilde = press * h - press * delta") ;
      syst.add_def(d, "Stilde = 3 * press * delta") ;
 
      // constraint equations
      syst.add_def(d, "eqP    = delta * D^i D_i P + 4piG / 2. * P^5 * Etilde") ;
      syst.add_def(d, "eqNP   = delta * D^i D_i NP - 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
 
      // definition for the baryonic mass integral
      syst.add_def(d, "intMb = P^6 * rho");
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst.add_def(d, "firstint = H + log(N)");
 
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");
 
      syst.add_def(d, "eqP = D^i D_i P");
      syst.add_def(d, "eqNP = D^i D_i NP");
      break;
    }
  }
 
  // add the constraint equations and demand continuity their normal derivative across domain boundaries
  space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
  space.add_eq(syst, "eqP = 0", "P", "dn(P)");
  
  // boundary conditions at infinity
  syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");

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
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");
 
  // if surface is resolved, fix the central enthalpy by one of these integrals
  if(!fixed) {
    space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
    space.add_eq_int_inf(syst, "integ(intMadm) = Madm");
  }
 
  // print the variation of the surface radius over the whole star
  if(rank == 0) {
    auto rs = Kadath::bco_utils::get_rmin_rmax(space, 1); 
    std::cout << "[Rmin, Rmax] : [" << rs[0] << ", " << rs[1] << "]\n";
  }
 
  double xo = 0.;
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;
 
  // solve until convergence is achieved
  while (!endloop) {  
    // do exactly one newton step, given the system above
    endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);
 
    update_config_quantities(loghc);
    // output files at this iteration and print diagnostics
    std::stringstream ss;
    ss << "norot_3d_";
    if(fixed) {
      ss << "fixed";
    }else {
      ss << "norot_bc";
    }
    ss << "_" << ite - 1;
    bconfig.set(QLMADM) = bconfig(MADM) ;
    bconfig.set_filename(ss.str());
    if (rank == 0) {
      print_diagnostics_norot(syst, ite, conv);
      std::cout << std::endl;
      if(bconfig.control(CHECKPOINT))
        checkpoint();
    }
 
    // update all coordinate fields, in case the domain extents have changed
    update_fields_co(cfields, coord_vectors, coord_scalars, 0.);

    ite++;
    check_max_iter_exceeded(rank, ite, conv);
  }
 
  bconfig.set_filename(converged_filename(stagename));
  if (rank == 0) {
    checkpoint();
  }
  return exit_status;
}

/** @}*/
}}