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
#include "kadath_adapted_polar.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "EOS/FUKA_EOS_Utilities.hh"
#include "coord_fields.hpp"
#include "mpi.h"

using namespace Kadath;
using namespace Kadath::FUKA_Config;
using namespace Kadath::FUKA_EOS;

// conversion from solar mass to km
constexpr double M2km = 1.4769994423016508;

template<class eos_t>
struct reader_3d {
  template<class config_t>
  void operator()(config_t& bconfig) {
    std::string spacein = bconfig.space_filename();
    // load domain decomposition and fields from binary file
    FILE* ff1 = fopen(spacein.c_str(), "r") ;
    Space_spheric_adapted space (ff1) ;
    Scalar conf   (space, ff1) ;
    Scalar lapse  (space, ff1) ;
    Vector shift  (space, ff1) ;
    Scalar logh   (space, ff1) ;

    Scalar diff_omega(space);
    if(bconfig.field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)){
      diff_omega = Scalar(space, ff1);
    }
    fclose(ff1) ;

    // number of dimensions, a cartesian type of basis and the flat bg metric
    int ndom = space.get_nbr_domains();
    Base_tensor basis(space, CARTESIAN_BASIS);
    Metric_flat fmet(space, basis);

    // fields depending on the coords
    CoordFields<Space_spheric_adapted> cf_generator(space);
    vec_ary_t coord_vectors {default_co_vector_ary(space)};

    // get origin of the system and initialize coordinate fields
    double xo = bco_utils::get_center(space,0);
    update_fields_co(cf_generator, coord_vectors, {}, xo);

    // central values of the matter fields
    double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);
    double hc = std::exp(loghc);
    double nc = EOS<eos_t,DENSITY>::get(hc);
    double pc = EOS<eos_t,PRESSURE>::get(hc);

    // minimal and maximal radius of the adapted surface domain
    auto [ rmin, rmax ] = bco_utils::get_rmin_rmax(space, 1);
    // inner radius of the nucleus
    double rin1 = bco_utils::get_radius(space.get_domain(0), OUTER_BC);

    // setup the system of equations to define derived quantities
    System_of_eqs syst(space, 0, ndom - 1);

    // flat background metric
    fmet.set_system(syst, "f");

    // constants
    syst.add_cst("4piG", bconfig(BCO_QPIG));
    syst.add_cst("PI"  , M_PI);

    // baryonic mass, dimensionless spin, central log enthalpy, angular frequency parameter
    syst.add_cst("Mb"  , bconfig(MB));
    syst.add_cst("chi" , bconfig(CHI));
    syst.add_cst("Hc"  , loghc);
    if(bconfig.field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)){
      syst.add_cst("ome", diff_omega);
    } else {
      syst.add_cst("ome" , bconfig(OMEGA));
    }
    syst.add_cst("Madm", bconfig(MADM));

    // coordinate dependent fields
    syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
    syst.add_cst("ex"  , *coord_vectors[EX]);
    syst.add_cst("ey"  , *coord_vectors[EX]);
    syst.add_cst("ez"  , *coord_vectors[EX]);
    syst.add_cst("einf", *coord_vectors[S_INF]);

    // gravitational fields
    syst.add_cst("P"   , conf);
    syst.add_cst("N"   , lapse);
    syst.add_cst("bet" , shift);

    // matter represented by log enthalpy
    syst.add_cst("H"   , logh);

    // common definitions
    syst.add_def("NP = P*N");
    syst.add_def("Ntilde = N / P^6");

    // beta combined with the rotational part
    if(std::isnan(bconfig.set(BVELY)))
      bconfig.set(BVELY) = 0.;
    syst.add_def("omega^i = bet^i + ome * mg^i");


    // extrinsic curvature
    syst.add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");
    syst.add_def (ndom-1, "intPx = A_i^j * ex_j * einf^i") ;
    syst.add_def (ndom-1, "intPy = A_i^j * ey_j * einf^i") ;
    syst.add_def (ndom-1, "intPz = A_i^j * ez_j * einf^i") ;

    // integrants at infinity
    syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 8. / PI");
    syst.add_def(ndom - 1, "intMadmalt = -dr(P) / 2 / PI");
    syst.add_def(ndom - 1, "intMadm = - einf^i * D_i P / 2 / PI");
    syst.add_def(ndom - 1, "intMk = (einf^i * D_i N - A_ij * einf^i * bet^j) / 4 / PI");

    syst.add_def ("dH = (ex^i * D_i H) / H");

    // setup operators for the EOS
    Param p;
    set_eos_ope_struct<eos_t>()(syst, p);

    // setup derived matter quantities
    for (int d = 0; d < ndom; d++) {
      switch (d) {
      case 0:
      case 1:
        syst.add_def(d, "h = exp(H)");

        syst.add_def(d, "U^i = omega^i / N");
        syst.add_def(d, "Usquare = P^4 * U_i * U^i");
        syst.add_def(d, "Wsquare = 1. / (1. - Usquare)");
        syst.add_def(d, "W = sqrt(Wsquare)");

        syst.add_def(d, "intMb = P^6 * rho(h) * W");
        syst.add_def(d, "firstint = H + log(N) - log(W)");
        syst.add_def(d, "intH  = P^6 * H * W") ;
        break;
      }
    }
    double Px       = space.get_domain(ndom-1)->integ(syst.give_val_def("intPx")()(ndom-1), OUTER_BC);
    double Py       = space.get_domain(ndom-1)->integ(syst.give_val_def("intPy")()(ndom-1), OUTER_BC);
    double Pz       = space.get_domain(ndom-1)->integ(syst.give_val_def("intPz")()(ndom-1), OUTER_BC);

    // ADM angular momentum at infinity
    Val_domain integJ(syst.give_val_def("intJ")()(ndom - 1));
    double J = space.get_domain(ndom - 1)->integ(integJ, OUTER_BC);

    // two (equivalent) computations of the ADM mass at infinity
    Val_domain integMadm(syst.give_val_def("intMadm")()(ndom - 1));
    double Madm = space.get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

    Val_domain integMadmalt(syst.give_val_def("intMadmalt")()(ndom - 1));
    double Madmalt = space.get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

    // Komar mass at infinity
    Val_domain integMk(syst.give_val_def("intMk")()(ndom - 1));
    double Mk = space.get_domain(ndom - 1)->integ(integMk, OUTER_BC);

    // baryonic mass as volume integral over the two stellar domains
    double baryonic_mass = syst.give_val_def("intMb")()(0).integ_volume() +
                          syst.give_val_def("intMb")()(1).integ_volume();
    double H_integral = syst.give_val_def("intH")()(0).integ_volume() +
                        syst.give_val_def("intH")()(1).integ_volume();

    syst.add_def ("intMsq = P^4");
    double A = space.get_domain(2)->integ(syst.give_val_def("intMsq")()(2), INNER_BC);
    double AR = sqrt(A / 4. / acos(-1.));

    auto npts = space.get_domain(1)->get_nbr_points();

    Index pos_eq (npts);
    pos_eq.set(0) = npts(0) - 1; /// Set to outer radius
    pos_eq.set(1) = npts(1) - 1; /// Set theta to be on the xy plane.

    Index pos_pole (npts);
    pos_pole.set(0) = npts(0) - 1; /// Set to outer radius

    // cout << "Cart (Pole): ";
    // for(auto i = 1; i <=3; ++i)
    //   cout << space.get_domain(1)->get_cart(i)(pos_pole) << ", ";
    // cout << endl;

    // cout << "Cart (Equitorial): ";
    // for(auto i = 1; i <=3; ++i)
    //   cout << space.get_domain(1)->get_cart(i)(pos_eq) << ", ";
    // cout << endl;
    double conf_eq = conf(1)(pos_eq);
    double Circumferential_R = conf_eq * conf_eq * space.get_domain(1)->get_radius()(pos_eq);
    double& CR = Circumferential_R;

    Scalar logh_dr(logh.der_r());
    double mass_shedding_parameter = logh_dr(1)(pos_eq) / logh_dr(1)(pos_pole);

    Point P(3);
    double central_dHdx = syst.give_val_def("dH")().val_point(P);
    double central_euler = syst.give_val_def("firstint")().val_point(P);

    #ifdef FORMAT
      #undef FORMAT
    #endif
    #define FORMAT std::setw(25) << std::right << std::setprecision(5) << std::fixed << std::showpos
    auto print_shells = [&](int dom_min, int dom_max)
    {
      int cnt = 1;
      for(int i = dom_min; i < dom_max; ++i) {
        std::string shell{"SHELL"+std::to_string(cnt)+" = "};
        std::cout << FORMAT << shell << bco_utils::get_radius(space.get_domain(i), OUTER_BC) << std::endl;
        cnt++;
      }
    };

    auto res_r = space.get_domain(0)->get_nbr_points()(0);
    auto res_t = space.get_domain(0)->get_nbr_points()(1);
    auto res_p = space.get_domain(0)->get_nbr_points()(2);

    // output to stdout
    std::cout << FORMAT << "RES = "  << "[" << res_r << "," << res_t << "," << res_p << "]\n"
              << FORMAT << "Coord R_IN = "  << rin1 << std::endl
              << FORMAT << "Coord R = "     << "[" << rmin << ", " << rmax << "]\n";
    print_shells(2, ndom-2);
    std::cout << FORMAT << "Coord R_OUT = " << bco_utils::get_radius(space.get_domain(ndom-2), OUTER_BC) << "\n\n";

    std::cout << FORMAT << "Areal R = "    << AR << " [" << AR * M2km << "km]\n"
              << FORMAT << "Circumferential R = " << CR << " [" << CR * M2km << "km]\n"
              << FORMAT << "Mass Shedding = " << mass_shedding_parameter << "\n"
              << FORMAT << "Baryonic Mass = " << baryonic_mass << std::endl
              << FORMAT << "ADM Mass = " << Madm << " [" << bconfig(MADM) << "]\n"
              << FORMAT << "ADM Momentum = " << J << std::endl
              << FORMAT << "Chi = " << J / Madm / Madm << " [" << bconfig(CHI) << "]\n"
              << FORMAT << "Omega = "<< bconfig(OMEGA) << std::endl
              << FORMAT << std::scientific << "Central Density = " << nc  << std::endl
              << FORMAT << std::scientific << "Central log(h) = " << loghc << std::endl
              << FORMAT << std::scientific << "Central Pressure = " << pc << std::endl
              << FORMAT << std::scientific << "Central dlog(h)/dx = " << central_dHdx << std::endl
              << FORMAT << std::scientific << "Central Euler Constant = "<< central_euler << std::endl
              << FORMAT << "Integrated log(h) = "    << H_integral << "\n\n";

    std::cout << FORMAT << "Mk = "   << Mk << std::scientific
              << ", Diff: " << 2. * fabs(Madm-Mk)/(Madm+Mk) << std::endl
              << FORMAT << "Px = "   << Px   << std::endl
              << FORMAT << "Py = "   << Py   << std::endl
              << FORMAT << "Pz = "   << Pz   << std::endl;
  }
};

int main(int argc, char **argv) {
  // initialize MPI
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    cerr << "Error starting MPI" << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  // expecting a configuration file on execution
  if(argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  // load the configuration
  std::string ifilename{argv[1]};
  kadath_config_boost<BCO_NS_INFO> bconfig(ifilename);
  // setup the EOS
  EOS_initialize::init(bconfig);
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);
  EOS_Function_Dispatcher::dispatch<reader_3d>(bconfig, eos_type, bconfig);
  // if(eos_type == "Cold_PWPoly") {
  //   using eos_t = Kadath::Margherita::Cold_PWPoly;

  //   EOS<eos_t,PRESSURE>::init(eos_file, h_cut);

  //   // call reader to output diagnostics
  //   if(bconfig(DIM) == 3) reader_3d<eos_t>(bconfig);
  // } else if(eos_type == "Cold_Table") {
  //   using eos_t = Kadath::Margherita::Cold_Table;

  //   const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
  //                           2000 : bconfig.eos<int>(INTERP_PTS);

  //   EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);

  //   // call reader to output diagnostics
  //   if(bconfig(DIM) == 3) reader_3d<eos_t>(bconfig);
  // } else {
  //   std::cerr << "Unknown EOSTYPE." << endl;
  //   std::_Exit(EXIT_FAILURE);
  // }

  MPI_Finalize();
  return EXIT_SUCCESS;
}


