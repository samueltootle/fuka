/*
 * Copyright 2021
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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
// FUKA includes
#include "Solvers/solver_startup.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"

#include "kadath_adapted.hpp"
#include "kadath_adapted_polar.hpp"

#include "mpi.h"
#include <sstream>
#define O(X) std::cout << "var " #X ": " << X << std::endl;
using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Solvers;

// // forward declarations
template<typename space_t, typename syst_t, typename config_t>
void print_diagnostics_norot(space_t const & space, syst_t const & syst, const config_t & bconfig, int iter, double err);
template<class eos_t, class config_t>
void update_config(config_t& bconfig, System_of_eqs& syst) ;
// template<typename space_t, typename syst_t, typename config_t>
// void print_diagnostics_rot(space_t const & space, syst_t const & syst, const config_t & bconfig, int const ite, double const conv);

template<typename config_t>
std::string converged_filename(const std::string&  stage, config_t bconfig);

template<class eos_t, typename config_t>
int driver(config_t& bconfig, std::string outputdir);

template<class eos_t, typename config_t>
int NS_solver_2d_norot (config_t& bconfig, bool fixed = false);

template<class eos_t, typename config_t>
int NS_solver_2d_uniform_rot (config_t& bconfig);

template<class eos_t, typename config_t>
int NS_solver_2d_differential_rot (config_t& bconfig);
// end forward declarations

int main(int argc, char **argv) {
  int rc = MPI_Init(&argc, &argv) ;
  if (rc!=MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl ;
    MPI_Abort(MPI_COMM_WORLD, rc) ;
  }
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

  using config_t = kadath_config_boost<BCO_ISO_NS_INFO>;
  using InitSolver = Initialize_Solver<config_t>;
  
  // Initialize static member variables
  InitSolver::input_configname = "initial_ns.info";
  InitSolver::bconfig;
  InitSolver::rank = rank;

  // Run initialize routine based on CLI arguments
  InitSolver::init_solver(argc, argv);
  config_t bconfig = InitSolver::bconfig;

  // specify output directory  
  std::string outputdir = InitSolver::outputdir;

  if(!InitSolver::example_setup) {
    // We now have to assume bconfig is a minimal config
    // that contains sequences _init/_final
    bconfig.open_config();
  
    // load and setup the EOS
    const double h_cut = bconfig.eos<double>(HCUT);
    const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      driver<eos_t>(bconfig, outputdir);
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      driver<eos_t>(bconfig, outputdir);
    } else { 
      std::cerr << "Unknown EOSTYPE." << endl;
      std::_Exit(EXIT_FAILURE);
    }
  }
  MPI_Finalize();

  return EXIT_SUCCESS;
} // end main()

template<class eos_t, typename config_t>
int driver(config_t& bconfig, std::string outputdir) {
  int rank = 0 ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
  int exit_status = EXIT_SUCCESS;

  bconfig.set_outputdir(outputdir);

  // specify which stages are enabled, pulled from the configuration file
  std::array<bool, NUM_STAGES> stage_enabled = bconfig.return_stages();
  
  // last stage of the enabled stages
  auto [ last_stage, last_stage_idx ] 
    = get_last_enabled(MSTAGE, stage_enabled);

  if(rank == 0){
      std::cout << "Last Stage Enabled: " << last_stage << std::endl;
  }

  if(stage_enabled[STAGES::NOROT_BC]) {
    if(bconfig.control(CONTROLS::SEQUENCES))
      exit_status = NS_solver_2d_norot<eos_t>(bconfig, true);
    exit_status = NS_solver_2d_norot<eos_t>(bconfig);
  }
  if(stage_enabled[STAGES::UNIFORM_ROT])
    exit_status = NS_solver_2d_uniform_rot<eos_t>(bconfig);
  if(stage_enabled[STAGES::DIFF_ROT])
    exit_status = NS_solver_2d_differential_rot<eos_t>(bconfig);
  return exit_status;
}

template<class eos_t, typename config_t>
int NS_solver_2d_norot (config_t& bconfig, bool fixed) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    std::cout << "###################################" << std::endl
              << "Non-rotating TOV Solver"      << std::endl
              << "###################################" << "\n\n";

  // convergence threshold
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy
  double loghc = std::log(bconfig(HC));
  std::string stage_name{"NOROT"};
  if(fixed) stage_name += "_FIXED";

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // load the fields defined on the space
	Scalar lapAterm   (space, ff1) ;
	Scalar nu  (space, ff1) ;
  Scalar logh   (space, ff1) ;
	fclose(ff1) ;

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

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
  
  // Fix stellar masss from central log specific enthalpy
  syst.add_cst("Hc", loghc);  

  // the variable fields to solve for: 
  // (log) specific enthalpy
  // log(lapse)
  // lapAterm = log(lapse) + log(A)
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("lapAterm", lapAterm);

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_cst("lev", level); 

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  // eq. 4.21 arxiv.org/abs/1003.5015v2
  // For non-rotating solutions, B = A and 4.21 reduces to
  // the follwing.
  syst.add_def(ndom - 1, "intMadm = -dr(A) / 4piG ");
  // eq. 4.15 arxiv.org/abs/1003.5015v2
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  // specific enthalpy from the logarithmic enthalpy, 
  // the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and 
  // pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  // See arXiv:2103.09911 
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      // source terms 3.36, 3.38, arxiv.org/abs/1003.5015v2
      // after rescaling by P/ rho
      syst.add_def(d, "E = press * h - press * delta");
      syst.add_def(d, "S = delta * 3 * press");
      syst.add_def(d, "Spp = press * delta");
      
      // constraint equations - the full system is rescaled by P/ rho
      // Here the source terms are not rescaled since they have already
      // been rescaled in their definition above
      // eq. 3.50 arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqnu = delta * ( lap(nu) + scal(grad(nu), grad(lapAterm)) ) - 4piG * A^2 * (E + S)") ;
      // eq. 3.52
      syst.add_def(d, "eqlapAterm = delta * ( lap2(lapAterm) + scal(grad(nu), grad(nu)) ) - 2 * 4piG * A^2 * Spp") ;
      
      // eq.3.51 this overdetermines the system, but it comes from the
      // constraint on B - see 3.16 for the full equation.
      // syst.add_def(d, "eqNA = dr(drNA) + 3 * divr(drNA) - 4 * 4piG * NA * A^2 * press") ;

      // definition for the baryonic mass integral
      // eq. 4.5 arxiv.org/abs/1003.5015v2, where we set mb = 1.
      syst.add_def(d, "intMb = A^3 * rho");
      
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      // eq. 3.98 arxiv.org/abs/1003.5015v2, lorentz factor W (Gamma in the paper) = 1
      syst.add_def(d, "firstint = H + log(N)");

      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(lapAterm))") ;
      syst.add_def(d, "eqlapAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))") ;
      break;
    }
  }
  // Add constraint equation so the System of equation ensuring
  // a continuous solution for each grid function and it's normal derivative
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqlapAterm=0", "lapAterm", "dn(lapAterm)");  
  
  if(fixed) {
    //Surface is set to a fixed radius
    syst.add_eq_bc(1, OUTER_BC, "lev=0");
  }
  else {
    // Surface is defined as vanishing log specific enthalpy
    syst.add_eq_bc(1, OUTER_BC, "H=0");
  }
  
  // Add first integral equation and ensure that the
  // log specific enthalpy at the origin = Hc
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");

  // Set boundary conditions
  // eq. 3.23 arxiv.org/abs/1003.5015v2
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapAterm=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // repeat until convergence is achieved
  while (!endloop) {        
    // do exactly one newton step
    endloop = syst.do_newton(conv_thres, conv);

    // output the data and diagnostics at this particular step
    std::stringstream ss;
    ss << "rot_3d_testing" << ite - 1 ;
    bconfig.set_filename(ss.str());

    if (rank == 0) {
      print_diagnostics_norot(space, syst, bconfig, ite, conv);
    }
    // count the steps
    ite++;
  }
  // Update the gravitational mass in the Config file
  bconfig.set(BCO_PARAMS::MADM) = 
    space.get_domain(ndom-1)->integ(syst.give_val_def("intMadm")()(ndom-1), OUTER_BC);
  // Update fluid quantities in the Config file
  update_config<eos_t>(bconfig, syst);
  // Give the config file a new filename
  bconfig.set_filename(converged_filename(stage_name, bconfig));
  bconfig.control(CONTROLS::SEQUENCES) = false;

  std::array<bool, NUM_STAGES> saved_stages = bconfig.return_stages();
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  
  // Disable all stages except for the one related to this solution
  stage_enabled.fill(false);
  stage_enabled[STAGES::UNIFORM_ROT] = true;
  bconfig.set_filename(converged_filename(stage_name, bconfig));

  if(rank == 0)
    bco_utils::save_to_file(space, bconfig, lapAterm, nu, logh);
  
  // reset stages
  stage_enabled = saved_stages;
  // disable current stage
  saved_stages[STAGES::NOROT_BC] = false;

  MPI_Barrier(MPI_COMM_WORLD);
  return EXIT_SUCCESS;
}

template<class eos_t, typename config_t>
int NS_solver_2d_uniform_rot (config_t& bconfig) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    std::cout << "###################################" << std::endl
              << "Uniformly Rotating models"      << std::endl
              << "Omega: " << bconfig(BCO_PARAMS::OMEGA) <<std::endl
              << "###################################" << "\n\n";

  // convergence threshold
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy
  double loghc = std::log(bconfig(HC));
  std::string stage_name{"ROT"};

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

  // Initialize fields in main scope before
  // checking if they should be loaded from file
  Scalar lapBterm(space);
  Scalar wrsint(space);

  // load the fields defined on the space
  // These are the same regardless of which 2D solution
  // you start from
	Scalar lapAterm(space, ff1) ;
	Scalar nu(space, ff1) ;
  Scalar logh(space, ff1) ;
  
  // Check if additional fields are present for the Bterm and omega term
  // to determine if we read them in from file.
  if(bconfig.set_field(BCO_FIELDS::LAP_BTERM) && bconfig.set_field(BCO_FIELDS::LAP_WTERM)) {
    if(rank == 0)
    std::cout << "Starting from shift and omega fields from file\n";
    lapBterm = Scalar(space, ff1);
    wrsint = Scalar(space, ff1);
  } else {
    if(rank == 0)
      std::cout << "Generating new shift and omega fields\n";
    Scalar w(space);
    w.annule_hard();
    w.std_base();
    wrsint = Scalar(w.mult_r().mult_sin_theta());
    
    Scalar B(exp(lapAterm - nu));
    Scalar N(exp(nu));
    Scalar tmp(N * B - 1);
    lapBterm = Scalar(tmp.mult_sin_theta().mult_r());
  }
  fclose(ff1) ;
  // Need to ensure this parameter is set such that
  // lap(wrsint) in the system of equations will
  // use \tilde{lap_3} as defined in eq. 3.21
  wrsint.affect_parameters();
  wrsint.set_parameters()->set_m_quant() = 1 ;
  wrsint.std_base();

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);

  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  syst.add_cst("Omega", bconfig(BCO_PARAMS::OMEGA));

  // Fix stellar masss from central log specific enthalpy
  syst.add_cst("Hc", loghc);

  // the variable fields to solve for: 
  // (log) specific enthalpy
  // log(lapse)
  // lapAterm = log(lapse) + log(A)
  // lapBterm = (Lapse * B - 1) * r * sin(theta)
  // wrsint = omega * r * sin(theta)
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("lapAterm", lapAterm);
  syst.add_var("lapBterm", lapBterm);
  syst.add_var("wrsint", wrsint);  

  // Useful definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_def("B = (divrsint(lapBterm) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");

  // define quantity to be integrated at infinity
  // two (in this case) equivalent definitions of ADM mass
  // as well as the Komar mass
  // eq. 4.21 arxiv.org/abs/1003.5015v2
  // Note: this equation does not give very good results.  For Madm = Mk,
  // intMad = -dr(B) / 4piG is much more accurate.  
  // FIXME - there must be a reason
  syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
  // eq. 4.15 arxiv.org/abs/1003.5015v2
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and 
  // pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  // See arXiv:2103.09911 
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:

      // Velocity terms 3.32, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "U = multrsint(B / N * (Omega - w))");
      syst.add_def(d, "Usq = U*U");
      // Lorentz factor 3.35, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "Wsq = 1 / (1 - Usq)");
      syst.add_def(d, "W = sqrt(Wsq)");
      
      // source terms 3.36 - 3.38, arxiv.org/abs/1003.5015v2
      // after rescaling by P/ rho
      syst.add_def(d, "E = Wsq * press * h - press * delta");
      syst.add_def(d, "Srrtt = press * delta");
      // Note: in the definition of the phi component of the pressure,
      // eq. 3.37, we ignore Brsin(theta) since it cancels analytically
      // with 1/Brsin(theta) that appears in eqw below
      syst.add_def(d, "pressp = (E + Srrtt) * U");
      syst.add_def(d, "Spp = delta * press * (1 + Usq) + E * Usq");
      syst.add_def(d, "S = 2 * Srrtt + Spp");
      
      // eq. 3.14, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqnu  = delta * lap(nu) + delta * scal(grad(nu), grad(nu + log(B))) "
                            "- delta * multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w)) "
                            "- 4piG * A^2 * (E + S)");
      
      // Note: the source term is different from eq 3.15, arxiv.org/abs/1003.5015v2 
      // by a factor of 1/Brsin(theta) as discussed above
      syst.add_def(d, "eqw = delta * lap(wrsint) - delta * multrsint(scal(grad(w), grad(nu - 3 * log(B))))"
                          "+ 4 * 4piG * N * A^2 / B * pressp");
      
      // eq. 3.16, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqlapBterm = delta * lap2(lapBterm) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
      
      // eq. 3.17, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqlapAterm = delta * lap2(lapAterm) + delta * scal(grad(nu), grad(nu))"
                      "- 3 * delta * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))"
                      "- 2 * 4piG * A^2 * Spp");
      
      // definition for the baryonic mass integral
      // eq. 4.5 arxiv.org/abs/1003.5015v2, where we set mb = 1.
      syst.add_def(d, "intMb = W * rho * A^2 * B * 4piG / 2");
      
      // first integral of the euler equation for a uniformly rotating star
      // eq. 3.98 arxiv.org/abs/1003.5015v2
      syst.add_def(d, "firstint = H + log(N) - 0.5 * log(Wsq)");
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");
      syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                      "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
      syst.add_def(d, "eqlapBterm = lap2(lapBterm)");
      syst.add_def(d, "eqlapAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))"
                "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))");
      break;
    }
  }
  
  // Add constraint equation so the System of equation ensuring
  // a continuous solution for each grid function and it's normal derivative
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqlapAterm=0", "lapAterm", "dn(lapAterm)");
  space.add_eq(syst, "eqlapBterm=0", "lapBterm", "dn(lapBterm)");
  space.add_eq(syst, "eqw=0", "wrsint", "dn(wrsint)");
  
  // Add first integral equation and ensure that the
  // log specific enthalpy at the origin = Hc
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");

  // Surface is defined as vanishing log specific enthalpy
  syst.add_eq_bc(1, OUTER_BC, "H=0");  

  // Set boundary conditions
  // eq. 3.23 arxiv.org/abs/1003.5015v2
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapAterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapBterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "wrsint=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // repeat until convergence is achieved
  while (!endloop) {        
    // do exactly one newton step
    endloop = syst.do_newton(conv_thres, conv);

    // output the data and diagnostics at this particular step
    std::stringstream ss;
    ss << "rot_2k_chkpt" << ite - 1 ;
    bconfig.set_filename(ss.str());

    if (rank == 0) {
      print_diagnostics_norot(space, syst, bconfig, ite, conv);
      // if(bconfig.control(CHECKPOINT))
        // bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
    // count the steps
    ite++;
  }
  // Update fluid quantities in Config file
  update_config<eos_t>(bconfig, syst);
  // Update gravitational mass in Config file
  bconfig.set(BCO_PARAMS::MADM) = 
  space.get_domain(ndom-1)->integ(syst.give_val_def("intMadm")()(ndom-1), OUTER_BC);

  std::array<bool, NUM_STAGES> saved_stages = bconfig.return_stages();
  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  
  // Disable all stages except for the one related to this solution
  stage_enabled.fill(false);
  stage_enabled[STAGES::UNIFORM_ROT] = true;
  bconfig.set_filename(converged_filename(stage_name, bconfig));
  // Ensure the additional fields are set in the Config file
  bconfig.set_field(BCO_FIELDS::LAP_BTERM) = true;
  bconfig.set_field(BCO_FIELDS::LAP_WTERM) = true;
  if(rank == 0) {
    std::cout << "Success!\n";
    bco_utils::save_to_file(space, bconfig, lapAterm, nu, logh, lapBterm, wrsint);
  }  
  // reset stages
  stage_enabled = saved_stages;
  // disable current stage
  saved_stages[STAGES::UNIFORM_ROT] = false;
  
  // Make sure all ranks sync
  MPI_Barrier(MPI_COMM_WORLD);  
  return EXIT_SUCCESS;
}

template<class eos_t, typename config_t>
int NS_solver_2d_differential_rot (config_t& bconfig) {
  // initialize MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // convergence threshold
  double& conv_thres = bconfig.seq_setting(SEQ_SETTINGS::PREC);
  // logarithm of the central enthalpy
  double loghc = std::log(bconfig(HC));
  std::string stage_name{"DIFFROT"};

  // (re)construct the numerical space
  std::string spacein = bconfig.space_filename();

  // load the space (and thus the domain setup)
	FILE* ff1 = fopen (spacein.c_str(), "r") ;
	Space_polar_adapted space (ff1) ;

  // number of domains defining the space
  int ndom = space.get_nbr_domains();

  // Sad tool to make system of equations work with constants
  Scalar one(space);
  one = 1;
  one.std_base();

  // Initialize fields in main scope before
  // checking if they should be loaded from file
  Scalar lapBterm(space);
  Scalar wrsint(space);
  Scalar Omega(space);
  Omega = (std::fabs(bconfig(BCO_PARAMS::OMEGA)) < 1e-7) ? 1e-7 : bconfig(BCO_PARAMS::OMEGA);
  Omega.std_base();

  // load the fields defined on the space - these should be there since NOROT stage
	Scalar lapAterm(space, ff1) ;
	Scalar nu(space, ff1) ;
  Scalar logh(space, ff1) ;
  
  // Check if additional fields are present for the Bterm and omega term
  // to determine if we read them in from file.
  if(bconfig.set_field(BCO_FIELDS::LAP_BTERM) && bconfig.set_field(BCO_FIELDS::LAP_WTERM)) {
    if(rank == 0)
      std::cout << "Starting from shift and omega fields from file\n";
    lapBterm = Scalar(space, ff1);
    wrsint = Scalar(space, ff1);
    // If we have a previous differential rotation solution...
    if(bconfig.set_field(BCO_FIELDS::DIFF_OMEGA)) {
      if(rank == 0)
        std::cout << "Starting from Omega field from file\n";
      Omega = Scalar(space, ff1);
    }
  } else {
    if(rank == 0)
      std::cout << "Generating new shift and omega fields\n";
    Scalar w(space);
    w.annule_hard();
    w.std_base();
    wrsint = Scalar(w.mult_r().mult_sin_theta());
    
    Scalar B(exp(lapAterm - nu));
    Scalar N(exp(nu));
    Scalar tmp(N * B - 1);
    lapBterm = Scalar(tmp.mult_sin_theta().mult_r());
  }
  fclose(ff1) ;
  // Need to ensure this parameter is set such that
  // lap(wrsint) in the system of equations will
  // use \tilde{lap_3} as defined in eq. 3.21
  wrsint.affect_parameters();
  wrsint.set_parameters()->set_m_quant() = 1 ;
  wrsint.std_base();

  auto npts = space.get_domain(1)->get_nbr_points();

  // Index locations for collocation points associated with:
  // - the origin
  // - equitorial point on the stellar surface, P(r=R_eq, theta=pi/2)
  // - polar point on the stellar surface P(r=R_pole, theta=0)
  Index pos_origin (npts);
  Index pos_eq (npts);
  pos_eq.set(0) = npts(0) - 1; /// Set to outer radius
  pos_eq.set(1) = npts(1) - 1; /// Set theta to be on the xy plane.

  Index pos_pole (npts);
  pos_pole.set(0) = npts(0) - 1; /// Set to outer radius

  auto adpt_dom = space.get_domain(1);
  // Initialize Radii variables from current solution
  double R0 = adpt_dom->get_radius()(pos_eq);
  double Rp = adpt_dom->get_radius()(pos_pole);
  
  // Initial Guess
  bconfig.set(BCO_PARAMS::RMID) = R0;

  // Extract Constants
  double diffAratio = bconfig.template diffrot<double>(DIFFROT_PARAMS::DIFF_ARATIO);
  double diffRratio = bconfig.template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO);

  // int q = bconfig.template diffrot<int>(DIFFROT_PARAMS::DIFF_Q);
  
  // Initialize KEH rotation law parameter A
  double diffA = diffAratio * R0;

  std::string firstint{"firstint = (H + log(N) - 0.5 * log(Wsq)) - 0.5 * j^2 / diffA^2"};

  if (rank == 0)
    std::cout << "###################################" << std::endl
              << "Differential Rotating models"      << std::endl
              << "Law: KEH\n"
              << "First Integral: " << firstint << std::endl
              << "Fixed A / R0: " << diffAratio << std::endl
              << "Initial Rp/Re: " << Rp / R0 << " {" << diffRratio << "}\n"
              << "Initial R0: " << R0 <<std::endl
              << "###################################" << "\n\n";

  // setup a system of equations
  System_of_eqs syst(space, 0, ndom - 1);

  // define numerical constants
  syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));
  syst.add_cst("Hc", loghc);

  // Constant parameters of the KEH differential rotation law
  syst.add_cst("diffAratio", diffAratio);
  syst.add_cst("Rratio", diffRratio);
  
  // Variable parameters to fix the differential rotation profile for KEH
  syst.add_var("diffA", diffA);
  syst.add_var("omec", bconfig(BCO_PARAMS::OMEGA));
  syst.add_var("R0", bconfig(BCO_PARAMS::RMID));

  // the variable fields to solve for: 
  // (log) specific enthalpy
  // log(lapse)
  // lapAterm = log(lapse) + log(A)
  // lapBterm = (Lapse * B - 1) * r * sin(theta)
  // wrsint = omega * r * sin(theta)
  // Omega = rotation profile
  syst.add_var("H", logh);
  syst.add_var("nu", nu);
  syst.add_var("lapAterm", lapAterm);
  syst.add_var("lapBterm", lapBterm);
  syst.add_var("wrsint", wrsint);
  syst.add_var("Omega", Omega);

  // Useful Metric definitions
  syst.add_def("N = exp(nu)");
  syst.add_def("A = exp(lapAterm - nu)");
  syst.add_def("B = (divrsint(lapBterm) + 1) / N");
  syst.add_def("w = divrsint(wrsint)");

  // Useful velocity definitions - needed for rotation law
  // as well as fluid definitions
  // Velocity terms 3.32, arxiv.org/abs/1003.5015v2
  syst.add_def("U = multrsint(B / N * (Omega - w))");
  syst.add_def("Usq = U*U");
  // Lorentz factor 3.35, arxiv.org/abs/1003.5015v2
  syst.add_def("Wsq = 1 /( 1 - Usq)");
  syst.add_def("W = sqrt(Wsq)");
  
  // KEH related definitions
  syst.add_cst("one", one);
  syst.add_def("diffAField = one * diffA");
  syst.add_def("r = multr(one)");
  syst.add_def("omeratio = omec / Omega");
  // This converges, but isn't correct
  // syst.add_def("j = Wsq / N * U");
  syst.add_def("j = Wsq / N * U * multrsint(B)");
  syst.add_def("omelaw = omec - j^2 / diffA^2");

  // This expression does not give accurate results when compared to M_komar and the computed
  // ADM Mass from the 3D code.  The deviation from the correct answer gets large with increasing
  // differential rotation profiles, but the source of the error is unknown
  // eq. 4.21 arxiv.org/abs/1003.5015v2
  // FIXME - there must be a reason
  // syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
  // Instead, the following definition gives very accurate comparisons with M_komar
  syst.add_def(ndom - 1, "intMadm = - (dr(B)) / 4piG ");
  syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");
  syst.add_def(ndom - 1, "intJ = -multrsint(multrsint(dr(w))) / 4 / 4piG");

  // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
  syst.add_def("h = exp(H)");

  // define the EOS operators
  Param p;
  syst.add_ope ("eps", &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope ("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope ("rho", &EOS<eos_t,DENSITY>::action, &p);

  // define rest-mass density, internal energy and 
  // pressure through the enthalpy
  syst.add_def("rho = rho(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("press = press(h)");

  // definition to rescale the equations
  // delta = p / rho
  // See arXiv:2103.09911 
  syst.add_def("delta = h - eps - 1.");

  for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      
      // source terms 3.36 - 3.38, arxiv.org/abs/1003.5015v2
      // after rescaling by P/ rho
      syst.add_def(d, "E = Wsq * press * h - press * delta");
      syst.add_def(d, "Srrtt = press * delta");
      // Note: in the definition of the phi component of the pressure,
      // eq. 3.37, we ignore Brsin(theta) since it cancels analytically
      // with 1/Brsin(theta) that appears in eqw below
      syst.add_def(d, "pressp = (E + Srrtt) * U");
      syst.add_def(d, "Spp = delta * press * (1 + Usq) + E * Usq");
      syst.add_def(d, "S = 2 * Srrtt + Spp");
      
      // eq. 3.14, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqnu  = delta * lap(nu) + delta * scal(grad(nu), grad(nu + log(B))) "
                            "- delta * multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w)) "
                            "- 4piG * A^2 * (E + S)");
      // Note: the source term is different from eq 3.15, arxiv.org/abs/1003.5015v2 
      // by a factor of 1/Brsin(theta) as discussed above
      syst.add_def(d, "eqw = delta * lap(wrsint) - delta * multrsint(scal(grad(w), grad(nu - 3 * log(B))))"
                          "+ 4 * 4piG * N * A^2 / B * pressp");
      
      // eq. 3.16, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqlapBterm = delta * lap2(lapBterm) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
      
      // eq. 3.17, arxiv.org/abs/1003.5015v2
      syst.add_def(d, "eqlapAterm = delta * lap2(lapAterm) + delta * scal(grad(nu), grad(nu))"
                      "- 3 * delta * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))"
                      "- 2 * 4piG * A^2 * Spp");
      
      // definition for the baryonic mass integral
      // eq. 4.5 arxiv.org/abs/1003.5015v2, where we set mb = 1.
      syst.add_def(d, "intMb = W * rho * A^2 * B * 4piG / 2");
      
      // first integral of the euler equation for a differentially rotating star
      // eq. 3.101 arxiv.org/abs/1003.5015v2
      syst.add_def(d, firstint.c_str());
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst.add_eq_full(d, "H = 0");

      syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                      "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqlapAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))"
                "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))");
      syst.add_def(d, "eqlapBterm = lap2(lapBterm)");
      syst.add_def(d, "eqw = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
      break;
    }
    if( d <= space.ADAPTED_INNER)
      syst.add_eq_full(d, "Omega - omelaw = 0");
    else
      syst.add_eq_full(d, "Omega = 0");
  }
  // Add constraint equation so the System of equation ensuring
  // a continuous solution for each grid function and it's normal derivative
  space.add_eq(syst, "eqnu=0", "nu", "dn(nu)");
  space.add_eq(syst, "eqlapAterm=0", "lapAterm", "dn(lapAterm)");
  space.add_eq(syst, "eqlapBterm=0", "lapBterm", "dn(lapBterm)");
  space.add_eq(syst, "eqw=0", "wrsint", "dn(wrsint)");

  // Add KEH constraints
  syst.add_eq_val(0, "diffAField/R0 - diffAratio", pos_origin);
  syst.add_eq_val(1, "r/R0 - 1", pos_eq);
  syst.add_eq_val(1, "r/R0 - Rratio", pos_pole);
  
  // Add first integral equation and ensure that the
  // log specific enthalpy at the origin = Hc
  syst.add_eq_first_integral(0, 1, "firstint", "H - Hc");
  
  // Surface is defined as vanishing log specific enthalpy
  syst.add_eq_bc(1, OUTER_BC, "H=0");  

  // Set boundary conditions
  // eq. 3.23 arxiv.org/abs/1003.5015v2
  syst.add_eq_bc(ndom - 1, OUTER_BC, "nu=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapAterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "lapBterm=0");
  syst.add_eq_bc(ndom - 1, OUTER_BC, "wrsint=0");
  
  // parameters for the solver loop
  bool endloop = false;
  int ite = 1;
  double conv;

  // repeat until convergence is achieved
  while (!endloop) {        
    // do exactly one newton step
    endloop = syst.do_newton(conv_thres, conv);

    // output the data and diagnostics at this particular step
    std::stringstream ss;
    ss << "rot_2k_chkpt" << ite - 1 ;
    bconfig.set_filename(ss.str());

    if (rank == 0) {
      print_diagnostics_norot(space, syst, bconfig, ite, conv);
      // if(bconfig.control(CHECKPOINT))
        // bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }
    // count the steps
    ite++;
  }
  // Update fluid quantities in Config file
  update_config<eos_t>(bconfig, syst);
  
  // Update gravitational mass in Config file
  bconfig.set(BCO_PARAMS::MADM) = 
  space.get_domain(ndom-1)->integ(syst.give_val_def("intMadm")()(ndom-1), OUTER_BC);

  std::array<bool, NUM_STAGES>& stage_enabled = bconfig.return_stages();
  stage_enabled.fill(false);
  stage_enabled[STAGES::DIFF_ROT] = true;
  bconfig.set_filename(converged_filename(stage_name, bconfig));
  bconfig.set_field(BCO_FIELDS::LAP_BTERM) = true;
  bconfig.set_field(BCO_FIELDS::LAP_WTERM) = true;
  bconfig.set_field(BCO_FIELDS::DIFF_OMEGA) = true;
  if(rank == 0) {
    std::cout << "Success!\n";
    bco_utils::save_to_file(space, bconfig, lapAterm, nu, logh, lapBterm, wrsint, Omega);
  }  
  MPI_Barrier(MPI_COMM_WORLD);
  return EXIT_SUCCESS;
}

template<class eos_t, class config_t>
void update_config(config_t& bconfig, System_of_eqs& syst) {
  auto const &  space = syst.get_space();

  auto rs = bco_utils::get_rmin_rmax(space, 1);
  auto h = syst.give_val_def("h")();
  auto hc = bco_utils::get_boundary_val(0, h, INNER_BC);

  // bconfig.set(BCO_PARAMS::RMID) = rs[0];
  bconfig.set(BCO_PARAMS::HC) = hc;
  bconfig.set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig(BCO_PARAMS::HC));
}

// diagnostics at runtime
#define FORMAT std::setw(10) << std::left << std::showpos 
template<typename space_t, typename syst_t, typename config_t>
void print_diagnostics_norot(space_t const & space, syst_t const & syst, 
    const config_t & bconfig, int iter, double err) {

  // total number of domains	
  int ndom = space.get_nbr_domains() ;

  // compute the baryonic mass at volume integral from the given integrant
  double baryonic_mass =
      syst.give_val_def("intMb")()(0).integ_volume() +
      syst.give_val_def("intMb")()(1).integ_volume();

  // compute the ADM mass as surface integral at infinity  
  Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
  double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

  // compute the Komar mass as surface integral at infinity
  Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
  double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

  // get the maximum and minimum coordinate radius along the surface,
  // i.e. the adapted domain boundary
  auto rs = bco_utils::get_rmin_rmax(space, 1);

  // output to standard output  
  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            << FORMAT << "Iter: " << iter << std::endl
            << FORMAT << "Error: " << err << std::endl
            << FORMAT << "Mb: " << baryonic_mass << std::endl
            << FORMAT << "Madm: " << Madm << std::endl
            << FORMAT << "Mk: " << Mk << " [" 
            << std::abs(Madm - Mk) / Madm << "]" << std::endl;
  std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n\n";
  std::cout.flags(f);
} // end print diagnostics norot


//standardized filename for each converged dataset at the end of each stage.
template<typename config_t>
std::string converged_filename(const std::string&  stage, config_t bconfig) {
  std::stringstream ss;
  ss << "NS2D";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << bconfig(MADM) << "." 
     << bconfig(CHI)<< "."   
     << std::setfill('0')  << std::setw(2) << bconfig(BCO_RES);
  return ss.str();
}
