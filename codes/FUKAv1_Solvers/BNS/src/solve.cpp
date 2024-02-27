/*
 * Copyright 2021
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author:
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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
 *
 * This implementation is designed to receive a NS to start from - i.e.
 * we are not creating a star from scratch.  This is useful for
 * running sequences of NSs.
*/
#include "kadath_bin_ns.hpp"
#include "EOS/EOS.hh"
#include "mpi.h"
#include "coord_fields.hpp"
#include "Configurator/config_binary.hpp"
#include "bco_utilities.hpp"
#include <sstream>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

using namespace Kadath;
using namespace Kadath::Margherita;
using namespace Kadath::FUKA_Config;

// forward declarations
template<class eos_t, typename config_t>
int BNS_solver (config_t& bconfig, std::string outputdir);

template<typename config_t>
void print_diagnostics(Space_bin_ns const & space, System_of_eqs const & syst,
    int const ite, double const conv, config_t& bconfig);

template<typename config_t>
std::stringstream converged_filename(const std::string&  stage, config_t bconfig);
// end forward declarations

int main(int argc, char **argv) {
  // initialize MPI
  int rc = MPI_Init(&argc, &argv);
  if (rc!=MPI_SUCCESS) {
	  cerr << "Error starting MPI" << endl;
	  MPI_Abort(MPI_COMM_WORLD, rc);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // usage
  if(argc < 2) {
    std::cerr << "Usage: ./solve /<path>/<ID base name>.info "\
                 "./<output_path>/" << endl;
    std::cerr << "e.g. ./solve initbin.info ./out" << endl;
    std::cerr << "Note: <output_path> is optional and defaults to <path>" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  std::string input_filename = argv[1];
  kadath_config_boost<BIN_INFO> bconfig(input_filename);

  std::string outputdir = bconfig.config_outputdir();
  if(argc > 2) outputdir = argv[2];

  // setup eos to before calling solver
  const double h_cut = bconfig.eos<double>(HCUT, BCO1);
  const std::string eos_file = bconfig.eos<std::string>(EOSFILE, BCO1);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);

  int return_status = EXIT_SUCCESS;
  if(eos_type == "Cold_PWPoly") {
    using eos_t = Kadath::Margherita::Cold_PWPoly;

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);

    // call solver
    return_status = BNS_solver<eos_t>(bconfig,outputdir);
  } else if(eos_type == "Cold_Table") {
    using eos_t = Kadath::Margherita::Cold_Table;

    const int interp_pts = (bconfig.eos<int>(INTERP_PTS, BCO1) == 0) ? \
                            2000 : bconfig.eos<int>(INTERP_PTS, BCO1);

    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);

    // call solver
    return_status = BNS_solver<eos_t>(bconfig,outputdir);
  }
  else {
    if(rank == 0)
      std::cerr << eos_type << " is not recognized.\n";
    std::_Exit(EXIT_FAILURE);
  }
  // end eos setup and solver

  MPI_Finalize();

  return return_status;
} // end main

template<class eos_t, typename config_t>
int BNS_solver (config_t& bconfig, std::string outputdir) {
  // get filename and stages of/from configuration file
  std::string input_filename = bconfig.config_filename();
  auto& stage_enabled = bconfig.return_stages();

  // print the configuration on rank0
	int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) std::cout << bconfig;

	// make sure the mass ratio is updated based on the component ADM masses
  // at inf separation, which should be given from the single star models
  bconfig(Q) = bconfig(MADM, BCO2) / bconfig(MADM, BCO1);

  // check if binary data file exists
  std::string kadath_filename = bconfig.space_filename();
  if(!fs::exists(kadath_filename)){
    if(rank == 0) {
      std::cerr << "File: " << kadath_filename << " not found.\n\n";
      std::_Exit(EXIT_FAILURE);
    }
  }

  // read in binary data, i.e. domain decomposistion and fields
	FILE* fin = fopen(kadath_filename.c_str(), "r");
	Space_bin_ns space(fin);
	Scalar conf  (space, fin);
	Scalar lapse (space, fin);
  Vector shift (space, fin);
	Scalar logh  (space, fin);
	Scalar phi   (space, fin);

	fclose(fin);

  // drop shift, if wanted
  if(bconfig.control(DELETE_SHIFT))
    shift.annule_hard();

  // set output directory and get range of stages
  bconfig.set_outputdir(outputdir);
  auto [ last_stage, last_stage_idx ] = get_last_enabled(MSTAGE, stage_enabled);

  if(rank == 0){
      std::cout << "Last Stage Enabled: " << last_stage << std::endl;
  }

  // setup cartesian basis for the given space
  Base_tensor basis (space, CARTESIAN_BASIS);

  // get the number of domains
	int ndom = space.get_nbr_domains();

  // get the center of the stellar domains as well as the whole domain decomposistion
 	double xc1 = bco_utils::get_center(space,space.NS1);
 	double xc2 = bco_utils::get_center(space,space.NS2);
  double xo  = bco_utils::get_center(space,ndom-1);

  // get central values of the logarithmic enthalpy
	double loghc1 = bco_utils::get_boundary_val(space.NS1, logh, INNER_BC);
	double loghc2 = bco_utils::get_boundary_val(space.NS2, logh, INNER_BC);

  // define a flat background metric
	Metric_flat fmet (space, basis);

  // define coordinate dependent fields
	CoordFields<Space_bin_ns> cfields(space);

  // initialize coordinate vector fields
  vec_ary_t coord_vectors = default_binary_vector_ary(space);

  update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

  // overview output
	if (rank==0) {
		std::cout << "=================================" << endl;
		std::cout << "BNS grav input" << endl;
		std::cout << "Distance: " << bconfig(DIST) << endl;
		std::cout << "Omega guess: " << bconfig(GOMEGA) << endl;
		std::cout << "Units: " << bconfig(QPIG) << endl;
		std::cout << "=================================" << endl;
	}

  // first stage, baryonic mass fixing but without explicit corrections to the linear momentum at infinty
  if(stage_enabled[TOTAL]) {
    if(rank == 0) std::cout << "############################" << std::endl
                            << "Total with mass fixing, chi = "
                            << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                            << " and q = " << bconfig(Q) << std::endl
                            << "############################" << std::endl;

    // setup a system of equations
    System_of_eqs syst (space, 0, ndom-1);

    // conformally flat background metric
    fmet.set_system(syst, "f");

    // setup EOS operators
    Param p;
    syst.add_ope ("eps"   , &EOS<eos_t,EPSILON>::action, &p);
    syst.add_ope ("press" , &EOS<eos_t,PRESSURE>::action, &p);
    syst.add_ope ("rho"   , &EOS<eos_t,DENSITY>::action, &p);
    syst.add_ope ("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);

    // constants
    syst.add_cst ("4piG"  , bconfig(QPIG));

    // baryonic mass and dimensionless spin are fixed input parameters
    syst.add_cst ("Mb1"   , bconfig(MB    , BCO1));
    syst.add_cst ("chi1"  , bconfig(CHI   , BCO1));

    syst.add_cst ("Mb2"   , bconfig(MB    , BCO2));
    syst.add_cst ("chi2"  , bconfig(CHI   , BCO2));

    // coordinate dependent (rotational) vector fields
    syst.add_cst ("mg"    , *coord_vectors[GLOBAL_ROT]);
    syst.add_cst ("mm"    , *coord_vectors[BCO1_ROT]);
    syst.add_cst ("mp"    , *coord_vectors[BCO2_ROT]);

    // cartesian basis vector fields
    syst.add_cst ("ex"    , *coord_vectors[EX]);
    syst.add_cst ("ey"    , *coord_vectors[EY]);
    syst.add_cst ("ez"    , *coord_vectors[EZ]);

    // flat space surface normals
    syst.add_cst ("sm"    , *coord_vectors[S_BCO1]);
    syst.add_cst ("sp"    , *coord_vectors[S_BCO2]);
    syst.add_cst ("einf"  , *coord_vectors[S_INF]);

    // ADM masses of each star at infinity are fixed input parameters
    // given by single star solutions (either TOV or rigidly rotating)
    syst.add_cst ("Madm1" , bconfig(MADM  , BCO1));
    syst.add_cst ("Madm2" , bconfig(MADM  , BCO2));

    // quasi-local measurement of the component ADM masses
    syst.add_var ("qlMadm1" , bconfig(QLMADM  , BCO1));
    syst.add_var ("qlMadm2" , bconfig(QLMADM  , BCO2));

    // central (logarithmic) enthalpy is a variable, fixed by the baryonic mass integral
    syst.add_var ("Hc1"   , loghc1);
    syst.add_var ("Hc2"   , loghc2);

    // center of mass on the x-axis connecting both companions
    // and orbital angular frequency parameter
    // in this case, both are fixed by the two central force-balance
    // equations
    syst.add_var ("xaxis" , bconfig(COM));
    syst.add_var ("ome"   , bconfig(GOMEGA));

    // the actual solution fields of the equations, i.e.
    // conformal factor, lapse, shift and (logarithmic) enthalpy
    syst.add_var ("P"     , conf);
    syst.add_var ("N"     , lapse);
    syst.add_var ("bet"   , shift);
    syst.add_var ("H"     , logh);

    // check if corotation is considered and adjust
    // rotational velocity components accordingly
    if(bconfig.control(COROT_BIN)) {
      // for a corotating binary, there is no angular frequency parameter
      // since stellar rotation is fixed by the orbital frequency
      bconfig.set(OMEGA, BCO1) = 0.;
      syst.add_cst ("omes1", bconfig(OMEGA,BCO1));

      bconfig.set(OMEGA, BCO2) = 0.;
      syst.add_cst ("omes2", bconfig(OMEGA,BCO2));
    } else {
      // in the arbitrary spinning / irrotational case
      // the irrotational part is defined by the velocity potential
      syst.add_var ("phi" , phi);

      // check if stellar omega is fixed or has to be solved
      // for a specific dimensionless spin
      if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst ("omes1", bconfig(FIXED_BCOMEGA,BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA,BCO1);
      } else {
        syst.add_var ("omes1", bconfig(OMEGA,BCO1));
      }

      if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst ("omes2", bconfig(FIXED_BCOMEGA,BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA,BCO2);
      } else {
        syst.add_var ("omes2", bconfig(OMEGA,BCO2));
      }

      // define the irrotational + spinning parts of the fluid velocity
      for(int d = space.NS1; d <= space.ADAPTED1; ++d){
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      for(int d = space.NS2; d <= space.ADAPTED2; ++d){
        syst.add_def(d, "s^i  = omes2 * mp^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
    }

    // some repeating definitions of the fluid part
    syst.add_def("h = exp(H)");
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("dHdlnrho = dHdlnrho(h)");
    syst.add_def("delta = h - eps - 1.");

    // some repeating definitions of the gravitational part
    syst.add_def ("NP = P*N");
    syst.add_def ("Ntilde = N / P^6");

    // orbital rotation vector field, corrected by the "center of mass" shift
    syst.add_def ("Morb^i = mg^i + xaxis * ey^i");
    // full shift vector field, incoporating the orbital part
    syst.add_def ("omega^i= bet^i + ome * Morb^i");

    // the conformal extrinsic curvature
    syst.add_def ("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

    // ADM linear momentum, surface integrant at infinity
    syst.add_def (ndom-1, "intPx = A_i^j * ex_j * einf^i");
    syst.add_def (ndom-1, "intPy = A_i^j * ey_j * einf^i");
    syst.add_def (ndom-1, "intPz = A_i^j * ez_j * einf^i");

    // quasi-local spin, surface integral outside the stellar matter distribution
    syst.add_def (space.ADAPTED1+1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
    syst.add_def (space.ADAPTED2+1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");

    // the actual equations, defined differently in the different domains
    for (int d=0; d<ndom; d++) {
      // if outside the stellar domains, without matter sources
      // resort to the source-free constraint equations
      // and set matter (and velocity potential) to zero
      if((d >= space.ADAPTED2+1) || d == space.ADAPTED1+1){
        if(!bconfig.control(COROT_BIN))
          syst.add_eq_full(d, "phi= 0");

        syst.add_eq_full(d, "H  = 0");

        syst.add_def(d,"eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst.add_def(d,"eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,"eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");

      }
      // in case of the domains harboring the two stars
      // define all derived matter related quantities
      // and enforce the transformed contraint equations (i.e. multiplied by p / rho)
      else {
        // 3-velocity and first integral of the Euler equation
        // in case of corotation
        if(bconfig.control(COROT_BIN)) {
          syst.add_def(d, "U^i    = omega^i / N");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "firstint = log(h * N / W)");
        }
        // 3-velocity, (approximate) first integral of the Euler equation
        // and velocity potential equation
        // in case of irrotational or spinning companions
        else {
          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "V^i    = N * U^i - omega^i");
          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
          syst.add_def(d, "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * W * V^i)");
        }
        // transformed source terms
        syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
        syst.add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare");
        syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

        // transformed constraint equations
        syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde");
        syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                               "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                               "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

        // baryonic mass volume integrant
        syst.add_def(d, "intMb  = P^6 * rho * W");
        // quasi-local ADM mass volume integrant
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");

      }
    }

    // add equations to the system and demand continuity
    // along the normal of the domains
    space.add_eq (syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq (syst, "eqP= 0", "P", "dn(P)");
    space.add_eq (syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    syst.add_eq_bc (ndom-1, OUTER_BC, "N=1");
    syst.add_eq_bc (ndom-1, OUTER_BC, "P=1");
    syst.add_eq_bc (ndom-1, OUTER_BC, "bet^i=0");

    // boundary conditions defining the boundary of the adapted domains, i.e. vanishing matter
    syst.add_eq_bc (space.ADAPTED1, OUTER_BC, "H = 0");
    syst.add_eq_bc (space.ADAPTED2, OUTER_BC, "H = 0");

    // in case of irrotational or spinning companions
    // solve also for the velocity potential
    if(!bconfig.control(COROT_BIN)) {
      // boundary conditions at the stellar surface,
      // i.e. where the velocity potential equation becomes
      // degenerate
      syst.add_eq_bc (space.ADAPTED1, OUTER_BC, "V^i * D_i H = 0");
      syst.add_eq_bc (space.ADAPTED2, OUTER_BC, "V^i * D_i H = 0");

      // add the equation and the matchings to the system
      // in case of the stellar domains
      for(int i = space.NS1; i < space.ADAPTED1; ++i) {
        syst.add_eq_vel_pot (i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching (i, OUTER_BC, "phi");
        syst.add_eq_matching (i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot   (space.ADAPTED1, 2, "eqphi = 0", "phi=0");

      for(int i = space.NS2; i < space.ADAPTED2; ++i) {
        syst.add_eq_vel_pot (i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching (i, OUTER_BC, "phi");
        syst.add_eq_matching (i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot   (space.ADAPTED2, 2, "eqphi = 0", "phi=0");

      // in case the frequency parameter of the star is not fixed by hand
      // solve for the given dimensionless spin through the quasi-local
      // spin surface integral
      if(std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_outer_sphere_one (syst, "integ(intS1) / Madm1 / Madm1 = chi1");

      if(std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_outer_sphere_two (syst, "integ(intS2) / Madm2 / Madm2 = chi2");
    }

    // force-balance equations at the center of each star
    // to fix the orbital frequency as well as the "center of mass"
    // (the latter can get very inaccurate in extreme configurations
    // with large residual linear momenta, for which the next stage is given)
    Index posori1 (space.get_domain(space.NS1)->get_nbr_points());
    syst.add_eq_val(space.NS1, "ex^i * D_i H", posori1);

    Index posori2 (space.get_domain(space.NS2)->get_nbr_points());
    syst.add_eq_val(space.NS2, "ex^i * D_i H", posori2);

    // add the first integral to get (approximate) hydrostatic equilibrium
    syst.add_eq_first_integral (space.NS1, space.ADAPTED1, "firstint", "H - Hc1");
    syst.add_eq_first_integral (space.NS2, space.ADAPTED2, "firstint", "H - Hc2");

    // fix the central enthalpy by baryonic mass volume integrals
    space.add_eq_int_volume (syst, space.NS1, space.ADAPTED1, "integvolume(intMb) = Mb1");
    space.add_eq_int_volume (syst, space.NS2, space.ADAPTED2, "integvolume(intMb) = Mb2");

    // compute a quasi-local approximation of the ADM component masses
    space.add_eq_int_volume (syst, space.NS1, space.ADAPTED1, "integvolume(intM) = qlMadm1");
    space.add_eq_int_volume (syst, space.NS2, space.ADAPTED2, "integvolume(intM) = qlMadm2");

    // print initial diagnostics
    if (rank==0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    // iterative solver variables
    bool endloop = false;
	  int ite = 1;
		double conv;

    // loop until desired convergence is achieved
		while (!endloop) {
      // do exactly one Newton step
			endloop = syst.do_newton(1e-8, conv);

      // recompute all coordinate dependent fields
      // to make sure that they are updated correctly along
      // with the changing adapted domains
      update_fields(cfields, coord_vectors, {}, xo, xc1, xc2, &syst);

      // generate output filename for this iteration
      std::stringstream ss;
      ss << "total_" << ite-1;
      bconfig.set_filename(ss.str());

      // print diagnostics and output configuration as well as the binary data
		  if (rank==0) {
        print_diagnostics(space, syst, ite, conv, bconfig);
  	    if(bconfig.control(CHECKPOINT))
          bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
		  }

			ite++;
		}

    // since the ADM mass at infinite separation is not known
    // a priori for corotating stars,
    // we use the quasi-local measurement to update the ADM masses
    // approximately
    if(bconfig.control(COROT_BIN)) {
      bconfig.set(MADM, BCO1) = bconfig(QLMADM, BCO1);
      bconfig.set(MADM, BCO2) = bconfig(QLMADM, BCO2);

      std::stringstream ss = converged_filename("TOTAL_COROT", bconfig);
      bconfig.set_filename(ss.str());
    } else {
      std::stringstream ss = converged_filename("TOTAL", bconfig);
      bconfig.set_filename(ss.str());
    }

    // output final configuration and binary data
    if (rank==0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
	}

  // second stage, baryonic mass fixing with explicit corrections to the linear momentum at infinty
  // but without solving the first integral again
  if(stage_enabled[TOTAL_BC]) {
    if(rank == 0) std::cout << "############################" << std::endl
                            << "Py fixing without first integral and fixed omega, chi = "
                            << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                            << " and q = " << bconfig(Q) << std::endl
                            << "with Py BC fixing" << std::endl
                            << "############################" << std::endl;

    // residual scaling factors for the ethalpy
    // used to correct the baryonic mass
    double H_scale1 = 0;
    double H_scale2 = 0;

    // "background", unscaled logarithmic enthalpy from the previous step
    Scalar logh_const(logh);
    logh_const.std_base();

    // setup a system of equations
    System_of_eqs syst (space, 0, ndom-1);

    // conformally flat background metric
    fmet.set_system(syst, "f");

    // setup EOS operators
    Param p;
    syst.add_ope ("eps"   , &EOS<eos_t,EPSILON>::action, &p);
    syst.add_ope ("press" , &EOS<eos_t,PRESSURE>::action, &p);
    syst.add_ope ("rho"   , &EOS<eos_t,DENSITY>::action, &p);
    syst.add_ope ("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);

    // constants
    syst.add_cst ("4piG"  , bconfig(QPIG));

    // baryonic mass and dimensionless spin are fixed input parameters
    syst.add_cst ("Mb1"   , bconfig(MB    , BCO1));
    syst.add_cst ("chi1"  , bconfig(CHI   , BCO1));

    syst.add_cst ("Mb2"   , bconfig(MB    , BCO2));
    syst.add_cst ("chi2"  , bconfig(CHI   , BCO2));

    // coordinate dependent (rotational) vector fields
    syst.add_cst ("mg"    , *coord_vectors[GLOBAL_ROT]);
    syst.add_cst ("mm"    , *coord_vectors[BCO1_ROT]);
    syst.add_cst ("mp"    , *coord_vectors[BCO2_ROT]);

    // cartesian basis vector fields
    syst.add_cst ("ex"    , *coord_vectors[EX]);
    syst.add_cst ("ey"    , *coord_vectors[EY]);
    syst.add_cst ("ez"    , *coord_vectors[EZ]);

    // flat space surface normals
    syst.add_cst ("sm"    , *coord_vectors[S_BCO1]);
    syst.add_cst ("sp"    , *coord_vectors[S_BCO2]);
    syst.add_cst ("einf"  , *coord_vectors[S_INF]);

    // ADM masses of each star at infinity are fixed input parameters
    // given by single star solutions (either TOV or rigidly rotating)
    syst.add_cst ("Madm1" , bconfig(MADM  , BCO1));
    syst.add_cst ("Madm2" , bconfig(MADM  , BCO2));

    // quasi-local measurement of the component ADM masses
    syst.add_var ("qlMadm1" , bconfig(QLMADM  , BCO1));
    syst.add_var ("qlMadm2" , bconfig(QLMADM  , BCO2));

    // residual scaling factors for the ethalpy
    syst.add_var ("Hscale1" , H_scale1);
    syst.add_var ("Hscale2" , H_scale2);

    // "center of mass" on the x-axis, connecting both stellar centers
    // fixed by the vanishing of the ADM linear momentum at infinity
    syst.add_var ("xaxis" , bconfig(COM));
    // no additional force-balance is computed,
    // the matter distribution is fixed modulo the scaling factos above,
    // therefore the orbital frequency is a constant similar to
    // eccentricity reduced ID
    syst.add_cst ("ome"   , bconfig(GOMEGA));

    // the actual solution fields of the equations, i.e.
    // conformal factor, lapse, shift and (logarithmic) enthalpy
    syst.add_var ("P"     , conf);
    syst.add_var ("N"     , lapse);
    syst.add_var ("bet"   , shift);
    // constant "background" enthalpy
    syst.add_cst ("Hconst", logh_const);

    // check if corotation is considered and adjust
    // rotational velocity components accordingly
    if(bconfig.control(COROT_BIN)) {
      // for a corotating binary, there is no angular frequency parameter
      // since stellar rotation is fixed by the orbital frequency
      bconfig.set(OMEGA, BCO1) = 0.;
      syst.add_cst ("omes1", bconfig(OMEGA,BCO1));

      bconfig.set(OMEGA, BCO2) = 0.;
      syst.add_cst ("omes2", bconfig(OMEGA,BCO2));
    } else {
      // in the arbitrary spinning / irrotational case
      // the irrotational part is defined by the velocity potential
      syst.add_var ("phi" , phi);

      // check if stellar omega is fixed or has to be solved
      // for a specific dimensionless spin
      if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst ("omes1", bconfig(FIXED_BCOMEGA,BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA,BCO1);
      } else {
        syst.add_var ("omes1", bconfig(OMEGA,BCO1));
      }

      if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst ("omes2", bconfig(FIXED_BCOMEGA,BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA,BCO2);
      } else {
        syst.add_var ("omes2", bconfig(OMEGA,BCO2));
      }

      // define the irrotational + spinning parts of the fluid velocity
      for(int d = space.NS1; d <= space.ADAPTED1; ++d){
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      for(int d = space.NS2; d <= space.ADAPTED2; ++d){
        syst.add_def(d, "s^i  = omes2 * mp^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
    }

    for (int d=0 ; d<ndom ; d++)
      // the enthalpy is equal to the constant part everywhere
      // outside of the stars, i.e. zero
      if((d >= space.ADAPTED2+1) || d == space.ADAPTED1+1)
        syst.add_def(d, "H  = Hconst");

    // inside the stars, it is the constant part
    // and a small correction factor to fix the correct
    // baryonic mass, due to slightly changing
    // fluid velocities
    for(int d = space.NS1; d <= space.ADAPTED1; ++d){
      syst.add_def(d, "H  = Hconst * (1. + Hscale1)");
    }
    for(int d = space.NS2; d <= space.ADAPTED2; ++d){
      syst.add_def(d, "H  = Hconst * (1. + Hscale2)");
    }

    // some repeating definitions of the fluid part
    syst.add_def("h = exp(H)");
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("dHdlnrho = dHdlnrho(h)");
    syst.add_def("delta = h - eps - 1.");

    // some repeating definitions of the gravitational part
    syst.add_def ("NP = P*N");
    syst.add_def ("Ntilde = N / P^6");

    // orbital rotation vector field, corrected by the "center of mass" shift
    syst.add_def ("Morb^i = mg^i + xaxis * ey^i");
    // full shift vector field, incoporating the orbital part
    syst.add_def ("omega^i= bet^i + ome * Morb^i");

    // the conformal extrinsic curvature
    syst.add_def ("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

    // ADM linear momentum, surface integrant at infinity
    syst.add_def (ndom-1, "intPx = A_i^j * ex_j * einf^i");
    syst.add_def (ndom-1, "intPy = A_i^j * ey_j * einf^i");
    syst.add_def (ndom-1, "intPz = A_i^j * ez_j * einf^i");

    // quasi-local spin, surface integral outside the stellar matter distribution
    syst.add_def (space.ADAPTED1+1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
    syst.add_def (space.ADAPTED2+1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");

    // the actual equations, defined differently in the different domains
    for (int d=0; d<ndom; d++) {
      // if outside the stellar domains, without matter sources
      // resort to the source-free constraint equations
      // and set matter (and velocity potential) to zero
      if((d >= space.ADAPTED2+1) || d == space.ADAPTED1+1){
        if(!bconfig.control(COROT_BIN))
          syst.add_eq_full(d, "phi= 0");

        syst.add_def(d,"eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst.add_def(d,"eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,"eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");

      }
      // in case of the domains harboring the two stars
      // define all derived matter related quantities
      // and enforce the transformed contraint equations (i.e. multiplied by p / rho)
      else {
        // 3-velocity and first integral of the Euler equation
        // in case of corotation
        if(bconfig.control(COROT_BIN)) {
          syst.add_def(d, "U^i    = omega^i / N");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "firstint = log(h * N / W)");
        }
        // 3-velocity, (approximate) first integral of the Euler equation
        // and velocity potential equation
        // in case of irrotational or spinning companions
        else {
          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "V^i    = N * U^i - omega^i");
          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
          syst.add_def(d, "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * W * V^i)");
        }
        // transformed source terms
        syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
        syst.add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare");
        syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

        // transformed constraint equations
        syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde");
        syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                               "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                               "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

        // baryonic mass volume integrant
        syst.add_def(d, "intMb  = P^6 * rho * W");
        // quasi-local ADM mass volume integrant
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");

      }
    }

    // add equations to the system and demand continuity
    // along the normal of the domains
    space.add_eq (syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq (syst, "eqP= 0", "P", "dn(P)");
    space.add_eq (syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    syst.add_eq_bc (ndom-1, OUTER_BC, "N=1");
    syst.add_eq_bc (ndom-1, OUTER_BC, "P=1");
    syst.add_eq_bc (ndom-1, OUTER_BC, "bet^i=0");

    // boundary conditions defining the boundary of the adapted domains, i.e. vanishing matter
    syst.add_eq_bc (space.ADAPTED1, OUTER_BC, "H = 0");
    syst.add_eq_bc (space.ADAPTED2, OUTER_BC, "H = 0");

    // in case of irrotational or spinning companions
    // solve also for the velocity potential
    if(!bconfig.control(COROT_BIN)) {
      // boundary conditions at the stellar surface,
      // i.e. where the velocity potential equation becomes
      // degenerate
      syst.add_eq_bc (space.ADAPTED1, OUTER_BC, "V^i * D_i H = 0");
      syst.add_eq_bc (space.ADAPTED2, OUTER_BC, "V^i * D_i H = 0");

      // add the equation and the matchings to the system
      // in case of the stellar domains
      for(int i = space.NS1; i < space.ADAPTED1; ++i) {
        syst.add_eq_vel_pot (i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching (i, OUTER_BC, "phi");
        syst.add_eq_matching (i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot   (space.ADAPTED1, 2, "eqphi = 0", "phi=0");

      for(int i = space.NS2; i < space.ADAPTED2; ++i) {
        syst.add_eq_vel_pot (i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching (i, OUTER_BC, "phi");
        syst.add_eq_matching (i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot   (space.ADAPTED2, 2, "eqphi = 0", "phi=0");

      // in case the frequency parameter of the star is not fixed by hand
      // solve for the given dimensionless spin through the quasi-local
      // spin surface integral
      if(std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_outer_sphere_one (syst, "integ(intS1) / Madm1 / Madm1 = chi1");

      if(std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_outer_sphere_two (syst, "integ(intS2) / Madm2 / Madm2 = chi2");
    }

    // volume integral fixing the baryonic mass by the matter scaling factor
    space.add_eq_int_volume (syst, space.NS1, space.ADAPTED1, "integvolume(intMb) = Mb1");
    space.add_eq_int_volume (syst, space.NS2, space.ADAPTED2, "integvolume(intMb) = Mb2");

    // volume integral computing a quasi-local approximation to the component ADM masses
    space.add_eq_int_volume (syst, space.NS1, space.ADAPTED1, "integvolume(intM) = qlMadm1");
    space.add_eq_int_volume (syst, space.NS2, space.ADAPTED2, "integvolume(intM) = qlMadm2");

    // enforcing the vanishing of the y-component of the ADM linear momentum at infinity,
    // fixing the "center of mass" shift on the x-axis
    space.add_eq_int_inf (syst, "integ(intPy) = 0");

    // print initial diagnostics
    if (rank==0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    // iterative solver variables
    bool endloop = false;
	  int ite = 1;
		double conv;

    // loop until desired convergence is achieved
		while (!endloop) {
      // do exactly one Newton step
			endloop = syst.do_newton(1e-8, conv);

      // recompute all coordinate dependent fields
      // to make sure that they are updated correctly along
      // with the changing adapted domains
      update_fields(cfields, coord_vectors, {}, xo, xc1, xc2, &syst);

      // overwrite logh with scaled version for later use and output
      logh = syst.give_val_def("H");

      // generate output filename for this iteration
      std::stringstream ss;
      ss << "total_bc_" << ite-1;
      bconfig.set_filename(ss.str());

      // print diagnostics and output configuration as well as the binary data
      if (rank==0) {
        print_diagnostics(space, syst, ite, conv, bconfig);
        std::cout << "Rescale Factors: " << H_scale1 << " " << H_scale2 << std::endl;
        if(bconfig.control(CHECKPOINT))
              bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
      }

      ite++;
    }

    // since the ADM mass at infinite separation is not known
    // a priori for corotating stars,
    // we use the quasi-local measurement to update the ADM masses
    // approximately
    if(bconfig.control(COROT_BIN)) {
      bconfig.set(MADM, BCO1) = bconfig(QLMADM, BCO1);
      bconfig.set(MADM, BCO2) = bconfig(QLMADM, BCO2);
      std::stringstream ss = converged_filename("TOTAL_BC_COROT", bconfig);
      bconfig.set_filename(ss.str());
    } else {
      std::stringstream ss = converged_filename("TOTAL_BC", bconfig);
      bconfig.set_filename(ss.str());
    }

    // output final configuration and binary data
    if(rank == 0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
    }

  // third stage, baryonic mass fixing with explicit corrections to the linear momentum at infinty
  // but without solving the first integral again the use of iteratively correct orbital
  // frequency and expansion coefficient
  if(stage_enabled[ECC_RED]) {
    // should the orbital frequency and expansion coefficient set by PN estimates?
    // This is useful either when no manual eccentricity iterations are carried out
    // or as a starting point for the latter.
    if(std::isnan(bconfig.set(ADOT)) ||
       std::isnan(bconfig.set(ECC_OMEGA)) ||
       bconfig.control(USE_PN)) {
      bco_utils::KadathPNOrbitalParams(bconfig, \
        bconfig(MADM, BCO1), bconfig(MADM,BCO2));

      if(rank == 0)
        std::cout << "### Using PN estimate for adot and omega! ###" << std::endl;
    }
    // get the eccentricity reduced orbital frequency as an input if given
    bconfig.set(GOMEGA) = bconfig(ECC_OMEGA);
    
    if (rank==0) {
      std::cout << "=================================" << endl;
      std::cout << "BNS grav input" << endl;
      std::cout << "Distance: " << bconfig(DIST) << endl;
      std::cout << "Omega: " << bconfig(GOMEGA) << endl;
      std::cout << "adot: " << bconfig(ADOT) << endl;
      std::cout << "Units: " << bconfig(QPIG) << endl;
      std::cout << "=================================" << endl;
    }

    // get a cartesian coordinate "vector field", i.e. the position vector wrt origin
    Vector CART(space, CON, basis);
    CART = cfields.cart();

    if(rank == 0) std::cout << "############################" << std::endl
                            << "Eccentricity reduction step with fixed omega, chi = "
                            << bconfig(CHI, BCO1) << "," << bconfig(CHI, BCO2)
                            << " and q = " << bconfig(Q) << std::endl
                            << "############################" << std::endl;

    // residual scaling factors for the ethalpy
    // used to correct the baryonic mass
    double H_scale1 = 0;
    double H_scale2 = 0;

    // "background", unscaled logarithmic enthalpy from the previous step
    Scalar logh_const(logh);
    logh_const.std_base();

    // set the shift of the "center of mass" along the y-axis
    // to zero, potentially fixing x-components of the ADM linear
    // momentum at infinity introduced by the eccentricity reduction
    // parameters
    if(std::isnan(bconfig.set(COMY)))
      bconfig.set(COMY) = 0.;

    // setup a system of equations
    System_of_eqs syst (space, 0, ndom-1);

    // conformally flat background metric
    fmet.set_system(syst, "f");

    // setup EOS operators
    Param p;
    syst.add_ope ("eps"   , &EOS<eos_t,EPSILON>::action, &p);
    syst.add_ope ("press" , &EOS<eos_t,PRESSURE>::action, &p);
    syst.add_ope ("rho"   , &EOS<eos_t,DENSITY>::action, &p);
    syst.add_ope ("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);

    // constants
    syst.add_cst ("4piG"  , bconfig(QPIG));

    // baryonic mass and dimensionless spin are fixed input parameters
    syst.add_cst ("Mb1"   , bconfig(MB    , BCO1));
    syst.add_cst ("chi1"  , bconfig(CHI   , BCO1));

    syst.add_cst ("Mb2"   , bconfig(MB    , BCO2));
    syst.add_cst ("chi2"  , bconfig(CHI   , BCO2));

    // coordinate dependent (rotational) vector fields
    syst.add_cst ("mg"    , *coord_vectors[GLOBAL_ROT]);
    syst.add_cst ("mm"    , *coord_vectors[BCO1_ROT]);
    syst.add_cst ("mp"    , *coord_vectors[BCO2_ROT]);

    // cartesian basis vector fields
    syst.add_cst ("ex"    , *coord_vectors[EX]);
    syst.add_cst ("ey"    , *coord_vectors[EY]);
    syst.add_cst ("ez"    , *coord_vectors[EZ]);

    // flat space surface normals
    syst.add_cst ("sm"    , *coord_vectors[S_BCO1]);
    syst.add_cst ("sp"    , *coord_vectors[S_BCO2]);
    syst.add_cst ("einf"  , *coord_vectors[S_INF]);

    // cartesian coordiante position vector wrt the origin
    syst.add_cst ("r"     , CART);

    // ADM masses of each star at infinity are fixed input parameters
    // given by single star solutions (either TOV or rigidly rotating)
    syst.add_cst ("Madm1" , bconfig(MADM  , BCO1));
    syst.add_cst ("Madm2" , bconfig(MADM  , BCO2));

    // quasi-local measurement of the component ADM masses
    syst.add_var ("qlMadm1" , bconfig(QLMADM  , BCO1));
    syst.add_var ("qlMadm2" , bconfig(QLMADM  , BCO2));

    // residual scaling factors for the ethalpy
    syst.add_var ("Hscale1" , H_scale1);
    syst.add_var ("Hscale2" , H_scale2);

    // "center of mass" on the x-axis, connecting both stellar centers
    // fixed by the vanishing of the ADM linear momentum at infinity
    syst.add_var ("xaxis" , bconfig(COM));
    // same on the y-axis in case finite momenta develope by
    // the eccentricity reduction parameters
    syst.add_var ("yaxis" , bconfig(COMY));

    // no additional force-balance is computed,
    // the matter distribution is fixed modulo the scaling factos above.
    // the orbital eccentricity and the expansion coefficient are thus an input parameter
    // (from iterative eccentricity reduction or PN estimates)
    syst.add_cst ("ome"   , bconfig(GOMEGA));
    syst.add_cst ("adot"  , bconfig(ADOT));

    // the actual solution fields of the equations, i.e.
    // conformal factor, lapse, shift and (logarithmic) enthalpy
    syst.add_var ("P"     , conf);
    syst.add_var ("N"     , lapse);
    syst.add_var ("bet"   , shift);
    // constant "background" enthalpy
    syst.add_cst ("Hconst", logh_const);

    // check if corotation is considered and adjust
    // rotational velocity components accordingly
    if(bconfig.control(COROT_BIN)) {
      // for a corotating binary, there is no angular frequency parameter
      // since stellar rotation is fixed by the orbital frequency
      bconfig.set(OMEGA, BCO1) = 0.;
      syst.add_cst ("omes1", bconfig(OMEGA,BCO1));

      bconfig.set(OMEGA, BCO2) = 0.;
      syst.add_cst ("omes2", bconfig(OMEGA,BCO2));
    } else {
      // in the arbitrary spinning / irrotational case
      // the irrotational part is defined by the velocity potential
      syst.add_var ("phi" , phi);

      // check if stellar omega is fixed or has to be solved
      // for a specific dimensionless spin
      if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1))) {
        syst.add_cst ("omes1", bconfig(FIXED_BCOMEGA,BCO1));
        bconfig.set(OMEGA, BCO1) = bconfig(FIXED_BCOMEGA,BCO1);
      } else {
        syst.add_var ("omes1", bconfig(OMEGA,BCO1));
      }

      if(!std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2))) {
        syst.add_cst ("omes2", bconfig(FIXED_BCOMEGA,BCO2));
        bconfig.set(OMEGA, BCO2) = bconfig(FIXED_BCOMEGA,BCO2);
      } else {
        syst.add_var ("omes2", bconfig(OMEGA,BCO2));
      }

      // define the irrotational + spinning parts of the fluid velocity
      for(int d = space.NS1; d <= space.ADAPTED1; ++d){
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      for(int d = space.NS2; d <= space.ADAPTED2; ++d){
        syst.add_def(d, "s^i  = omes2 * mp^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
    }

    for (int d=0 ; d<ndom ; d++)
      // the enthalpy is equal to the constant part everywhere
      // outside of the stars, i.e. zero
      if((d >= space.ADAPTED2+1) || d == space.ADAPTED1+1)
        syst.add_def(d, "H  = Hconst");

    // inside the stars, it is the constant part
    // and a small correction factor to fix the correct
    // baryonic mass, due to slightly changing
    // fluid velocities
    for(int d = space.NS1; d <= space.ADAPTED1; ++d){
      syst.add_def(d, "H  = Hconst * (1. + Hscale1)");
    }
    for(int d = space.NS2; d <= space.ADAPTED2; ++d){
      syst.add_def(d, "H  = Hconst * (1. + Hscale2)");
    }

    // some repeating definitions of the fluid part
    syst.add_def("h = exp(H)");
    syst.add_def("press = press(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("rho = rho(h)");
    syst.add_def("dHdlnrho = dHdlnrho(h)");
    syst.add_def("delta = h - eps - 1.");

    // some repeating definitions of the gravitational part
    syst.add_def ("NP = P*N");
    syst.add_def ("Ntilde = N / P^6");

    // orbital rotation vector field, corrected by the "center of mass" shift
    syst.add_def ("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");
    // full shift vector field, incoporating the orbital part
    syst.add_def ("rcom^i = r^i - xaxis * ex^i + yaxis * ey^i");
    syst.add_def ("omega^i= bet^i + ome * Morb^i + adot * rcom^i");

    // the conformal extrinsic curvature
    syst.add_def ("A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

    // ADM linear momentum, surface integrant at infinity
    syst.add_def (ndom-1, "intPx = A_i^j * ex_j * einf^i");
    syst.add_def (ndom-1, "intPy = A_i^j * ey_j * einf^i");
    syst.add_def (ndom-1, "intPz = A_i^j * ez_j * einf^i");

    // quasi-local spin, surface integral outside the stellar matter distribution
    syst.add_def (space.ADAPTED1+1, "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
    syst.add_def (space.ADAPTED2+1, "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");

    // the actual equations, defined differently in the different domains
    for (int d=0; d<ndom; d++) {
      // if outside the stellar domains, without matter sources
      // resort to the source-free constraint equations
      // and set matter (and velocity potential) to zero
      if((d >= space.ADAPTED2+1) || d == space.ADAPTED1+1){
        if(!bconfig.control(COROT_BIN))
          syst.add_eq_full(d, "phi= 0");

        syst.add_def(d,"eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst.add_def(d,"eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst.add_def(d,"eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");

      }
      // in case of the domains harboring the two stars
      // define all derived matter related quantities
      // and enforce the transformed contraint equations (i.e. multiplied by p / rho)
      else {
        // 3-velocity and first integral of the Euler equation
        // in case of corotation
        if(bconfig.control(COROT_BIN)) {
          syst.add_def(d, "U^i    = omega^i / N");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "Wsquare= 1 / (1 - Usquare)");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "firstint = log(h * N / W)");
        }
        // 3-velocity, (approximate) first integral of the Euler equation
        // and velocity potential equation
        // in case of irrotational or spinning companions
        else {
          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");
          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");
          syst.add_def(d, "V^i    = N * U^i - omega^i");
          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");
          syst.add_def(d, "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * W * V^i)");
        }
        // transformed source terms
        syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
        syst.add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare");
        syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

        // transformed constraint equations
        syst.add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde");
        syst.add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                               "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst.add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                               "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

        // baryonic mass volume integrant
        syst.add_def(d, "intMb  = P^6 * rho * W");
        // quasi-local ADM mass volume integrant
        syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");

      }
    }

    // add equations to the system and demand continuity
    // along the normal of the domains
    space.add_eq (syst, "eqNP= 0", "N", "dn(N)");
    space.add_eq (syst, "eqP= 0", "P", "dn(P)");
    space.add_eq (syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");

    // boundary conditions at infinity
    syst.add_eq_bc (ndom-1, OUTER_BC, "N=1");
    syst.add_eq_bc (ndom-1, OUTER_BC, "P=1");
    syst.add_eq_bc (ndom-1, OUTER_BC, "bet^i=0");

    // boundary conditions defining the boundary of the adapted domains, i.e. vanishing matter
    syst.add_eq_bc (space.ADAPTED1, OUTER_BC, "H = 0");
    syst.add_eq_bc (space.ADAPTED2, OUTER_BC, "H = 0");

    // in case of irrotational or spinning companions
    // solve also for the velocity potential
    if(!bconfig.control(COROT_BIN)) {
      // boundary conditions at the stellar surface,
      // i.e. where the velocity potential equation becomes
      // degenerate
      syst.add_eq_bc (space.ADAPTED1, OUTER_BC, "V^i * D_i H = 0");
      syst.add_eq_bc (space.ADAPTED2, OUTER_BC, "V^i * D_i H = 0");

      // add the equation and the matchings to the system
      // in case of the stellar domains
      for(int i = space.NS1; i < space.ADAPTED1; ++i) {
        syst.add_eq_vel_pot (i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching (i, OUTER_BC, "phi");
        syst.add_eq_matching (i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot   (space.ADAPTED1, 2, "eqphi = 0", "phi=0");

      for(int i = space.NS2; i < space.ADAPTED2; ++i) {
        syst.add_eq_vel_pot (i, 2, "eqphi = 0", "phi=0");
        syst.add_eq_matching (i, OUTER_BC, "phi");
        syst.add_eq_matching (i, OUTER_BC, "dn(phi)");
      }
      syst.add_eq_vel_pot   (space.ADAPTED2, 2, "eqphi = 0", "phi=0");

      // in case the frequency parameter of the star is not fixed by hand
      // solve for the given dimensionless spin through the quasi-local
      // spin surface integral
      if(std::isnan(bconfig.set(FIXED_BCOMEGA, BCO1)))
        space.add_eq_int_outer_sphere_one (syst, "integ(intS1) / Madm1 / Madm1 = chi1");

      if(std::isnan(bconfig.set(FIXED_BCOMEGA, BCO2)))
        space.add_eq_int_outer_sphere_two (syst, "integ(intS2) / Madm2 / Madm2 = chi2");
    }

    // volume integral fixing the baryonic mass by the matter scaling factor
    space.add_eq_int_volume (syst, space.NS1, space.ADAPTED1, "integvolume(intMb) = Mb1");
    space.add_eq_int_volume (syst, space.NS2, space.ADAPTED2, "integvolume(intMb) = Mb2");

    // volume integral computing a quasi-local approximation to the component ADM masses
    space.add_eq_int_volume (syst, space.NS1, space.ADAPTED1, "integvolume(intM) = qlMadm1");
    space.add_eq_int_volume (syst, space.NS2, space.ADAPTED2, "integvolume(intM) = qlMadm2");

    // enforcing the vanishing of the x- and y-component of the ADM linear momentum at infinity,
    // fixing the "center of mass" shift on the x,y-axis
    space.add_eq_int_inf (syst, "integ(intPx) = 0");
    space.add_eq_int_inf (syst, "integ(intPy) = 0");

    // print initial diagnostics
    if (rank==0)
      print_diagnostics(space, syst, 0, 0, bconfig);

    // iterative solver variables
    bool endloop = false;
	  int ite = 1;
		double conv;

    // loop until desired convergence is achieved
		while (!endloop) {
      // do exactly one Newton step
			endloop = syst.do_newton(1e-8, conv);

      // recompute all coordinate dependent fields
      // to make sure that they are updated correctly along
      // with the changing adapted domains
      update_fields(cfields, coord_vectors, {}, xo, xc1, xc2, &syst);

      // overwrite logh with scaled version for later use and output
      logh = syst.give_val_def("H");

      // generate output filename for this iteration
      std::stringstream ss;
      ss << "total_ecc_" << ite-1;
      bconfig.set_filename(ss.str());

      // print diagnostics and output configuration as well as the binary data
      if (rank==0) {
        print_diagnostics(space, syst, ite, conv, bconfig);

        std::cout << "Rescale Factors: " << H_scale1 << " " << H_scale2 << std::endl;

        if(bconfig.control(CHECKPOINT))
              bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
      }

      ite++;
    }

    // since the ADM mass at infinite separation is not known
    // a priori for corotating stars,
    // we use the quasi-local measurement to update the ADM masses
    // approximately
    if(bconfig.control(COROT_BIN)) {
      bconfig.set(MADM, BCO1) = bconfig(QLMADM, BCO1);
      bconfig.set(MADM, BCO2) = bconfig(QLMADM, BCO2);

      std::stringstream ss = converged_filename("ECC_RED_COROT", bconfig);
      bconfig.set_filename(ss.str());
    } else {
      std::stringstream ss = converged_filename("ECC_RED", bconfig);
      bconfig.set_filename(ss.str());
    }

    // output final configuration and binary data
    if (rank==0)
      bco_utils::save_to_file(space, bconfig, conf, lapse, shift, logh, phi);

  }
	return EXIT_SUCCESS;
} // end BNS_Solver

// runtime diagnostics
template<typename config_t>
void print_diagnostics(Space_bin_ns const & space, System_of_eqs const & syst,
    int const ite, double const conv, config_t& bconfig) {

  // get number of domains
	int ndom = space.get_nbr_domains();

  // volume integration of baryonic mass and quasi-local ADM mass for both stars  
  double baryonic_mass1 = 0.;
  double ql_mass1 = 0.;
  for(int d = space.NS1; d <= space.ADAPTED1; ++d){
    baryonic_mass1 += syst.give_val_def("intMb")()(d).integ_volume();
    ql_mass1 += syst.give_val_def("intM")()(d).integ_volume();
  }

  double baryonic_mass2 = 0.;
  double ql_mass2 = 0.;
  for(int d = space.NS2; d <= space.ADAPTED2; ++d){
    baryonic_mass2 += syst.give_val_def("intMb")()(d).integ_volume();
    ql_mass2 += syst.give_val_def("intM")()(d).integ_volume();
  }

  // surface integration of the ADM linear momentum at infinity
  Val_domain integPx(syst.give_val_def("intPx")()(ndom - 1));
  double Px = space.get_domain(ndom - 1)->integ(integPx, OUTER_BC);

  Val_domain integPy(syst.give_val_def("intPy")()(ndom - 1));
  double Py = space.get_domain(ndom - 1)->integ(integPy, OUTER_BC);

  Val_domain integPz(syst.give_val_def("intPz")()(ndom - 1));
  double Pz = space.get_domain(ndom - 1)->integ(integPz, OUTER_BC);

  // surface integration of the quasi-local (dimensionless) spin angular momentum of both stars
  Val_domain integS1(syst.give_val_def("intS1")()(space.ADAPTED1+1));
  double S1 = space.get_domain(space.ADAPTED1+1)->integ(integS1, OUTER_BC);
  double chi1 = S1 / bconfig(MADM, BCO1) / bconfig(MADM, BCO1);
  double chiql1 = S1 / ql_mass1 / ql_mass1;

  Val_domain integS2(syst.give_val_def("intS2")()(space.ADAPTED2+1));
  double S2 = space.get_domain(space.ADAPTED2+1)->integ(integS2, OUTER_BC);
  double chi2 = S2 / bconfig(MADM, BCO2) / bconfig(MADM, BCO2);
  double chiql2 = S2 / ql_mass2 / ql_mass2;

  #define FORMAT std::setw(13) << std::left << std::showpos

  // print mininmal and maximal radius of the surface adapted domains
  auto print_rs = [&](int dom) {
    auto rs = bco_utils::get_rmin_rmax(space, dom);
    std::cout << FORMAT << "NS-Rs: " << rs[0] << " " << rs[1] << std::endl;
  };

  std::ios_base::fmtflags f( std::cout.flags() );
  std::cout << "=======================================" << std::endl
            // iteration of the iterative solver
            << FORMAT << "Iter: " << ite << std::endl
            // remaining residual of all the equations
            << FORMAT << "Error: " << conv << std::endl
            // orbital angular frequency parameter
            << FORMAT << "Omega: " << bconfig(GOMEGA) << std::endl
            // "center of mass" shift
            << FORMAT << "Axis: " << bconfig(COM) << std::endl
            // ADM linear momentum at infinity
            // i.e. residual momentum of the spacetime
            << FORMAT << "P: " << "[" << Px << ", " << Py << ", " << Pz << "]\n\n"

            // baryonic mass
            << FORMAT << "NS1-Mb: " << baryonic_mass1 << std::endl
            // quasi-local ADM mass
            << FORMAT << "NS1-Madm_ql: " << ql_mass1 << std::endl
            // quasi-local spin angular momentum
            << FORMAT << "NS1-S: " << S1 << std::endl
            // quasi-local dimensionless spin
            << FORMAT << "NS1-Chi: " << chi1 << std::endl
            // quasi-local dimensionless spin, normalized by the quasi-local ADM mass
            << FORMAT << "NS1-Chi_ql: " << chiql1 << std::endl
            // angular frequency paramter of the stellar rotation
            << FORMAT << "NS1-OmegaS: " << bconfig(OMEGA, BCO1) << "\n";
  print_rs(space.ADAPTED1);

  std::cout << "\n"
            << FORMAT << "NS2-Mb: " << baryonic_mass2 << std::endl
            << FORMAT << "NS2-Madm_ql: " << ql_mass2 << std::endl
            << FORMAT << "NS2-S: " << S2 << std::endl
            << FORMAT << "NS2-Chi: " << chi2 << std::endl
            << FORMAT << "NS2-Chi_ql: " << chiql2 << std::endl
            << FORMAT << "NS2-OmegaS: " << bconfig(OMEGA, BCO2) << std::endl;
  print_rs(space.ADAPTED2);
  std::cout.flags(f);
  std::cout << "=======================================" << "\n\n";
} // end print_diagnostics


// standardized filename for each converged dataset at the end of each stage.
template<typename config_t>
std::stringstream converged_filename(const std::string&  stage, config_t bconfig) {
  std::stringstream ss;
  ss << "converged_BNS";
  if(stage != "") ss  << "_" << stage << ".";
  else ss << ".";
  ss << bconfig(DIST)      << "."
     << bconfig(CHI, BCO1) << "."
     << bconfig(CHI, BCO2) << "."
     << bconfig(MADM, BCO1)+bconfig(MADM, BCO2) << ".q"
     << bconfig(Q)         << "."
     << std::setfill('0')  << std::setw(2) << bconfig(BIN_RES);
  return ss;
}
