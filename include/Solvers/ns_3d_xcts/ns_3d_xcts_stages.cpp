#include "mpi.h"
#include "bco_utilities.hpp"

/**
 * \addtogroup Stages
 * \ingroup NS_XCTS
 * @{*/
namespace Kadath {
namespace FUKA_Solvers {
  
template<class eos_t, typename config_t, typename space_t>
int ns_3d_xcts_solver<eos_t, config_t, space_t>::norot_stage(bool fixed) {
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
                << "FIXED stage is deprecated" << std::endl
                << "############################" << std::endl;
    std::__throw_runtime_error("Fixed stage is deprecated.\n");
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

  update_fields_co(cfields, coord_vectors, {}, 0.);

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);
  syst.add_var("H"   , logh);
  syst_init(syst);
 
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

  // if the surface is resolved, define it to be where the matter vanishes
  syst.add_eq_bc(1, OUTER_BC, "H = 0");
  
  std::string central_fixing_definition{"H - Hc"};
  if(seq && seq->is_set()) {
    auto idx{std::get<0>(seq->get_indices())};
    switch(idx) {
      case BCO_PARAMS::HC:
        syst.add_cst("Hc", loghc);
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        break;
      case BCO_PARAMS::NC:
        syst.add_cst("Nc", bconfig(BCO_PARAMS::NC));
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        central_fixing_definition = "rho - Nc";
        break;
      default:
        syst.add_var("Hc", loghc);
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_cst("Madm", bconfig(BCO_PARAMS::MADM));
        break;
    }
  } else {
    syst.add_var("Hc", loghc);
    /// Future deprecate
    if(bconfig.control(MB_FIXING)) {
      syst.add_cst("Mb"  , bconfig(BCO_PARAMS::MB));
      syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
    }
    else {
      syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
      syst.add_cst("Madm", bconfig(BCO_PARAMS::MADM));
    }
  }
  // first integral in the innermost domains with non-zero matter content
  // and condition on the central value, either fixed directly or by the
  // integral below
  syst.add_eq_first_integral(0, 1, "firstint", central_fixing_definition.c_str());
 
  // constrain stellar mass by these integrals and the central log enthalpy
  space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");
  space.add_eq_int_inf(syst, "integ(intMadm) = Madm");
 
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
 
    update_config_quantities(bco_utils::get_boundary_val(0, logh, INNER_BC));
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
    update_fields_co(cfields, coord_vectors, {}, 0.);

    ite++;
    check_max_iter_exceeded(rank, ite, conv);
  }
 
  bconfig.set_filename(converged_filename(stagename));
  if (rank == 0) {
    checkpoint();
  }
  return exit_status;
}

template<class eos_t, typename config_t, typename space_t>
int ns_3d_xcts_solver<eos_t, config_t, space_t>::uniform_rot_stage() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // We use `config_filename()` vs `config_filename_abs()` since
  // `solution_exists` will probe the HOME_KADATH/COs directory
  auto const current = bconfig.config_filename();
  if(!bconfig.control(RESOLVE) && !seq && solution_exists("TOTAL_BC")) {
    if(rank == 0)
      std::cout << "Solved previously: " \
                << bconfig.config_filename_abs() << std::endl;
    return (current == bconfig.config_filename()) ? \
      EXIT_SUCCESS : RELOAD_FILE;
  }
  
  const int max_iter = bconfig.seq_setting(MAX_ITER);

  double loghc = std::log(bconfig(HC));
  double xo = 0.0;
  update_fields_co(cfields, coord_vectors, {}, xo);
  
  if (rank == 0 && bconfig.control(MB_FIXING))
    std::cout << "###################################" << std::endl
              << "Rotating - with fixed Baryonic Mass" << std::endl
              << "###################################" << std::endl;
  else if (rank == 0) 
    std::cout << "###################################" << std::endl
              << "Rotating - with fixed ADM Mass"      << std::endl
              << "###################################" << std::endl;

  System_of_eqs syst(space, 0, ndom - 1);
  syst.add_var("H"   , logh);
  syst_init(syst);
  syst.add_var("bet" , shift);
  
  syst.add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / "
               "2. / Ntilde");

  // ADM Angular momentum
  syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
  // Quasi-local spin angular momentum
  syst.add_def(2,"intS = A_ij * mg^i * sm^j / 2. / 4piG") ;
  
  std::string central_fixing_definition{"H - Hc"};
  std::string spin_fixing_definition{"integ(intJ) - chi * Madm * Madm = 0"};
  
  if(seq && seq->is_set()) {
    auto idx{std::get<0>(seq->get_indices())};
    switch(idx) {
      case BCO_PARAMS::HC:
        syst.add_cst("Hc", loghc);
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        syst.add_cst("chi" , bconfig(CHI));
        syst.add_var("ome" , bconfig(OMEGA));
        break;
      case BCO_PARAMS::NC:
        syst.add_cst("Nc", bconfig(BCO_PARAMS::NC));
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        syst.add_cst("chi" , bconfig(CHI));
        syst.add_var("ome" , bconfig(OMEGA));
        central_fixing_definition = "rho - Nc";
        break;
      case BCO_PARAMS::OMEGA:
        syst.add_cst("Nc", bconfig(BCO_PARAMS::NC));
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        syst.add_var("chi" , bconfig(CHI));
        syst.add_cst("ome" , bconfig(OMEGA));
        central_fixing_definition = "rho - Nc";
        break;
      case BCO_PARAMS::JADM:
        spin_fixing_definition = "integ(intJ) - Jadm = 0";
        // Fixed MADM and CHI
        syst.add_var("Hc", loghc);
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_cst("Madm", bconfig(BCO_PARAMS::MADM));
        syst.add_cst("Jadm", bconfig(BCO_PARAMS::JADM));
        syst.add_var("ome" , bconfig(OMEGA));
        bconfig.set(BCO_PARAMS::CHI) = bconfig(BCO_PARAMS::JADM) / bconfig(BCO_PARAMS::MADM) / bconfig(BCO_PARAMS::MADM);
        break;
      default:
        // Fixed MADM and CHI
        syst.add_var("Hc", loghc);
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_cst("Madm", bconfig(BCO_PARAMS::MADM));
        syst.add_cst("chi" , bconfig(CHI));
        syst.add_var("ome" , bconfig(OMEGA));
        break;
    }
  } else {
    syst.add_var("Hc", loghc);
    syst.add_cst("chi" , bconfig(CHI));
    syst.add_var("ome" , bconfig(OMEGA));
    /// Future deprecate
    if(bconfig.control(MB_FIXING)) {
      syst.add_cst("Mb"  , bconfig(BCO_PARAMS::MB));
      syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
    }
    else {
      syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
      syst.add_cst("Madm", bconfig(BCO_PARAMS::MADM));
    }
  }

  syst.add_def("omega^i = bet^i + ome * mg^i");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    case 0:
    case 1:
      syst.add_def(d, "U^i = omega^i / N");
      syst.add_def(d, "Usquare = P^4 * U_i * U^i");
      syst.add_def(d, "Wsquare = 1. / (1. - Usquare)");
      syst.add_def(d, "W = sqrt(Wsquare)");

      syst.add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
      syst.add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare") ;
      syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

      syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde") ;
      syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                             "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
      syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                             "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

      syst.add_def(d, "intMb = P^6 * rho(h) * W");
      syst.add_def(d, "firstint = H + log(N) - log(W)");

      break;
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
      syst.add_def(d, "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
      syst.add_def(d, "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                      "A^ij * D_j Ntilde");
      break;
    }
  }

  space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
  space.add_eq(syst, "eqP= 0", "P", "dn(P)");
  space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

  syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

  syst.add_eq_bc(1, OUTER_BC, "H = 0");

  syst.add_eq_first_integral(0, 1, "firstint", central_fixing_definition.c_str());
  space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");

  space.add_eq_int_inf(syst, spin_fixing_definition.c_str());
  space.add_eq_int_inf(syst, "integ(intMadm) = Madm");

  if (rank == 0)
      print_diagnostics(syst, 0, 0);
  bool endloop = false;
  int ite = 1;
  double conv;
  while (!endloop) {
    endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

    update_config_quantities(bco_utils::get_boundary_val(0, logh, INNER_BC));
    std::stringstream ss;
    ss << "rot_3d_total" << ite - 1 ;
    bconfig.set(QLMADM) = bconfig(MADM);
    bconfig.set_filename(ss.str());
    if (rank == 0) {
      print_diagnostics(syst, ite, conv);
      if(bconfig.control(CHECKPOINT))
        checkpoint();
    }
    update_fields_co(cfields, coord_vectors, {}, xo, &syst);
    ite++;
    check_max_iter_exceeded(rank, ite, conv);
  }
  
  bconfig.set_filename(converged_filename("TOTAL_BC"));
  if (rank == 0) {
    checkpoint();
  }
  return EXIT_SUCCESS;
}

template<class eos_t, typename config_t, typename space_t>
int ns_3d_xcts_solver<eos_t, config_t, space_t>::binary_boost_stage(
  kadath_config_boost<BIN_INFO>& binconfig, const size_t bco) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int max_iter = bconfig.seq_setting(MAX_ITER);
  
  // generate filename string unique to this binary setup
  std::stringstream stage_ss;
  stage_ss << "BIN_BOOST" << "_" << binconfig(DIST) << "_" << binconfig(GOMEGA);
  auto const boost_converged_filename{stage_ss.str()};
  if(!bconfig.control(RESOLVE) && solution_exists(boost_converged_filename)) {
    if(rank == 0)
      std::cout << "Solved previously: " 
                << bconfig.config_filename_abs() << std::endl;
    return EXIT_SUCCESS;
  }

  // Since we assume a comoving frame in the COM of the Binary
  // we need to using the same system of equation as in the Binary
  // this includes splitting the fluid velocity into an irrotational
  // term plus a rotation term
  Scalar phi(space);
  phi.annule_hard();
  phi.std_base();
  
  // flag needs to be set in order to be used during import into
  // a BNS or BHNS setup
  bconfig.set_field(PHI) = true;

  double H_scale = 0;
  Scalar logh_const(logh);
  logh_const.std_base();

  double xo = 0.0;
  
	double loghc1 = Kadath::bco_utils::get_boundary_val(0, logh, INNER_BC);
  
  if(rank == 0) std::cout << "############################" << std::endl
                          << "Binary boosted NS" << std::endl
                          << "############################" << std::endl;

  update_fields_co(cfields, coord_vectors, {}, xo);
  System_of_eqs syst(space, 0, ndom - 1);
  
  // irrotational part of the fluid velocity
  syst.add_var("phi" , phi);
  syst.add_var("qlMadm", bconfig(QLMADM));
  syst.add_cst("Hconst", logh_const);
  syst.add_var("Hscale" , H_scale);
  syst.add_cst("ex", *coord_vectors[EX]);
  syst.add_cst("ey", *coord_vectors[EY]);
  for (int d=0 ; d<ndom ; d++)
    // the enthalpy is equal to the constant part everywhere
    // outside of the stars, i.e. zero
    if(d > 1)
      syst.add_def(d, "H  = Hconst");
    else
      syst.add_def(d, "H  = Hconst * (1. + Hscale)");
  
  syst_init(syst);

  syst.add_cst("chi" , bconfig(CHI));
  syst.add_var("ome" , bconfig(OMEGA));
  
  // MB and MADM are now fixed as they would be in a binary
  syst.add_cst("Mb"  , bconfig(MB));
  syst.add_cst("Madm", bconfig(MADM));
  
  syst.add_var("bet" , shift);
  
  // boost based on binary config
  syst.add_cst("omega_boost", binconfig(GOMEGA));
  syst.add_def("B^i = bet^i + omega_boost * (mg^i)");

  syst.add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / "
               "2. / Ntilde");

  syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
  syst.add_def("intS = A_ij * mg^i * sm^j / 2 / 4piG") ;

  // define the irrotational + spinning parts of the fluid velocity
  for(int d = 0; d <= 1; ++d){
    syst.add_def(d, "s^i  = ome * mg^i");
    syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
  }

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    case 0:
    case 1:
      // definitions for the fluid 3-velocity
      // and its Lorentz factor
      syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
      syst.add_def(d, "W      = sqrt(Wsquare)");
      syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
      syst.add_def(d, "Usquare= P^4 * U_i * U^i");
      syst.add_def(d, "V^i    = N * U^i - B^i");
      syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
      syst.add_def(d, "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * W * V^i)");

      // rescaled sources and constraint equations        
      syst.add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
      syst.add_def(d, "Stilde = 4 * press * delta + (Etilde + press * delta) * Usquare") ;
      syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

      syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde") ;
      syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                             "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
      syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                             "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

      // integrant of the baryonic mass integral
      syst.add_def(d, "intMb = P^6 * rho(h) * W");

      break;
    default:
      // outside the star the matter is absent and the sources are zero
      //syst.add_eq_full(d, "H = 0");
      syst.add_eq_full(d, "phi = 0");

      syst.add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
      syst.add_def(d, "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
      syst.add_def(d, "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                      "A^ij * D_j Ntilde");
      break;
    }
  }

  space.add_eq(syst, "eqNP= 0", "N", "dn(N)");
  space.add_eq(syst, "eqP= 0", "P", "dn(P)");
  space.add_eq(syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

  syst.add_eq_bc(ndom - 1, OUTER_BC, "N=1");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "P=1");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

  syst.add_eq_bc(1, OUTER_BC, "H = 0");
  syst.add_eq_bc(1, OUTER_BC, "V^i * D_i H = 0");
  // in case of the stellar domains
  syst.add_eq_vel_pot(0, 2, "eqphi = 0", "phi=0");
  syst.add_eq_matching(0, OUTER_BC, "phi");
  syst.add_eq_matching(0, OUTER_BC, "dn(phi)");
  syst.add_eq_vel_pot(1, 2, "eqphi = 0", "phi=0");

  space.add_eq_int(syst, 2, OUTER_BC, "integ(intS) - chi * Madm * Madm = 0");
  space.add_eq_int_inf(syst, "integ(intMadmalt) = qlMadm");
  
  space.add_eq_int_volume(syst, 2, "integvolume(intMb) = Mb");

  if (rank == 0)
      print_diagnostics(syst, 0, 0);
  bool endloop = false;
  int ite = 1;
  double conv;
  while (!endloop) {
    endloop = syst.do_newton(bconfig.seq_setting(PREC), conv);

    std::stringstream ss;
    ss << "ns_bin_boost_" << ite - 1 ;
    bconfig.set_filename(ss.str());
    
    // overwrite logh with scaled version for later use and output
    logh = syst.give_val_def("H");
    
    if (rank == 0) {
      print_diagnostics(syst, ite, conv);
      if(bconfig.control(CHECKPOINT)) {
        // Manual save here since we need to save PHI
        Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
      }
    }
    update_fields_co(cfields, coord_vectors, {}, xo, &syst);
    ite++;
    check_max_iter_exceeded(rank, ite, conv);
  }

  bconfig.set_filename(converged_filename(boost_converged_filename));
  if (rank == 0) {
    Kadath::bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
  }
  return EXIT_SUCCESS;
}
/** @}*/
}}