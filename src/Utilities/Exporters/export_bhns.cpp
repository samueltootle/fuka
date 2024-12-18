#include <kadath.hpp>
#include <EOS/EOS.hh>
#include <coord_fields.hpp>
#include <Configurator/config_binary.hpp>
#include <exporter_utilities.hpp>
#include <cmath>
#include "EOS/FUKA_EOS_Utilities.hh"

using namespace Kadath;
using namespace export_utils;
using namespace Kadath::FUKA_Config;

std::array<std::vector<double>,NUM_OUT> KadathExportBHNS(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn,
                                                        double const interpolation_offset, int const interp_order,
                                                        double const delta_r_rel) {
  std::string filename(fn);

  kadath_config_boost<BIN_INFO> bconfig(filename);

  using namespace Kadath::FUKA_EOS;
  EOS_initialize::init(bconfig, NODES::BCO1);

  // get const EOS information
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);

  /* file containing KADATH fields must have same name as config file
   * with only the extension being different */
  std::string kadath_filename = filename.substr(0,filename.size()-5) + ".dat";

	FILE* fin = fopen(kadath_filename.c_str(), "r") ;
	Space_bhns space(fin) ;
	Scalar conf  (space, fin) ;
	Scalar lapse (space, fin) ;
  Vector shift (space, fin) ;
	Scalar logh  (space, fin) ;
	Scalar phi   (space, fin) ;

	fclose(fin) ;

  std::vector<std::reference_wrapper<const Scalar>> quants;
  for (int i = 0; i < NUM_QUANTS; ++i)
    quants.push_back(std::cref(conf));

  quants[PSI] = std::cref(conf);
  quants[ALP] = std::cref(lapse);

  Base_tensor basis(shift.get_basis());

  double rbh;
  int ndom = space.get_nbr_domains();

 	Index center_pos (space.get_domain(space.NS)->get_nbr_points()) ;
 	double xm = space.get_domain(space.NS)->get_cart(1)(center_pos);
 	double xp = space.get_domain(space.BH)->get_cart(1)(center_pos);
  double xo  = space.get_domain(ndom-1)->get_cart(1)(center_pos);

  /* Setup system of equations for Aij and Matter terms */
  Metric_flat fmet(space, basis);

  CoordFields<Space_bhns> cfields(space);

  Vector global_rot = cfields.rot_z();
  Vector star_rot = cfields.rot_z(xm);
  Vector ez = cfields.e_cart(3);
  Vector ey = cfields.e_cart(2);
  Vector ex = cfields.e_cart(1);
  Vector esurf = cfields.e_rad();

  System_of_eqs syst(space, 0, ndom - 1);

  fmet.set_system(syst, "f");

  Param p;
  EOS_Function_Dispatcher::dispatch<set_eos_ope_struct>(bconfig, eos_type, syst, p);

  syst.add_cst("4piG", bconfig(QPIG));
  syst.add_cst("PI", M_PI);

  syst.add_cst("omes1", bconfig(OMEGA,BCO1));

  syst.add_cst("Mg", global_rot);
  syst.add_cst("ome", bconfig(GOMEGA));
  syst.add_cst("xaxis", bconfig(COM));
  syst.add_cst("yaxis", bconfig(COMY));

  syst.add_cst("ex", ex);
  syst.add_cst("ey", ey);

  syst.add_cst("m1", star_rot);

  syst.add_cst("P", conf);
  syst.add_cst("N", lapse);
  syst.add_cst("bet", shift);
  syst.add_cst("phi", phi);

  syst.add_cst("H", logh);

  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");

  syst.add_def ("Morb^i = Mg^i + xaxis * ey^i + yaxis * ex^i");
  syst.add_def ("B^i= bet^i + ome * Morb^i");

  for(int d = space.NS; d <= space.ADAPTEDNS; ++d){
    syst.add_def  (d, "s^i  = omes1 * m1^i");
  }

  syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");

  syst.add_def("h = exp(H)");

  for(int d = 0; d < ndom; ++d) {
    if(d <= space.ADAPTEDNS)
      syst.add_def(d, "eta_i = D_i phi + P^4 * s_i");
    else
      syst.add_def(d, "eta_i = D_i phi");
  }

  syst.add_def("Wsquare = eta^i * eta_i / h^2 / P^4 + 1.");
  syst.add_def("W = sqrt(Wsquare)");

  syst.add_def("U^i = eta^i / P^4 / h / W");

  Tensor A(syst.give_val_def("A"));
  Index ind(A);

  quants[AXX] = std::cref(A(ind));
  ind.inc();
  quants[AXY] = std::cref(A(ind));
  ind.inc();
  quants[AXZ] = std::cref(A(ind));
  ind.inc();
  ind.inc();
  quants[AYY] = std::cref(A(ind));
  ind.inc();
  quants[AYZ] = std::cref(A(ind));
  ind.inc();
  ind.inc();
  ind.inc();
  quants[AZZ] = std::cref(A(ind));

  quants[H] = std::cref(logh);

  Vector vel_kad(syst.give_val_def("U"));

  quants[UX] = std::cref(vel_kad(1));
  quants[UY] = std::cref(vel_kad(2));
  quants[UZ] = std::cref(vel_kad(3));

  quants[BETX] = std::cref(shift(1));
  quants[BETY] = std::cref(shift(2));
  quants[BETZ] = std::cref(shift(3));
  /* end Aij, shift, and matter setup */

  /* setup for filling in BH with junk */
  Index I2(space.get_domain(space.BH+2)->get_radius().get_conf().get_dimensions());
  rbh = space.get_domain(space.BH+2)->get_radius()(I2);
  /* end BH setup */

  std::array<std::vector<double>,NUM_OUT> out;
  for(auto& v : out)
    v.resize(npoints);

  for (int i = 0; i < npoints; ++i) {
    std::vector<double> quant_vals(NUM_QUANTS);

    // get number of dimensions and construct index
    int ndim = 3;
    double x_shifted = xx[i] - bconfig(COM); //shifting carpet [0,0] to the initial center of mass
    double y_shifted = yy[i] - bconfig(COMY); //shifting carpet [0,0] to the initial center of mass
    double xxp = x_shifted - xp;
    double r2yz = y_shifted * y_shifted + zz[i] * zz[i];

    //Radius measurement centered on the BH
    double r_plus  = std::sqrt(xxp * xxp + r2yz);

    //lambda function for interpolation inside BH.
    auto interp_f = [&](auto& ah_r, auto& extrap_r, auto& bh_ori) {
      double theta = std::acos(zz[i] / extrap_r);
      double Phi = std::atan2(y_shifted, (x_shifted - bh_ori)); // atan2 is needed here

      std::vector<double> r_points(interp_order);
      for (int j = 0; j < interp_order; j++) {
        r_points[j] = (1. + interpolation_offset) * (1. + j * delta_r_rel) * ah_r;
      }

      //std::vector<double> psi6(interp_order);

      for (int k = 0; k < NUM_QUANTS; ++k) {
        std::vector<double> vals(interp_order);

        for (int j = 0; j < interp_order; j++) {
          vals[j] = quants[k].get().val_point(
              point_spherical(r_points[j], theta, Phi, bh_ori));
        }
        if(k == H)
          quant_vals[k] = 0;
        else if(k == UX || k == UY || k == UZ)
          quant_vals[k] = 0;
        else
          quant_vals[k] =
              lagrange_gen_k(interp_order, extrap_r, r_points.data(), vals.data());
      }
    };

    if(r_plus <= (1. + interpolation_offset) * rbh) {
      interp_f(rbh, r_plus, xp);
    } else {
      // construct Kadath point
      Point abs_coords(ndim);
      abs_coords.set(1) = x_shifted;
      abs_coords.set(2) = y_shifted;
      abs_coords.set(3) = zz[i];

      for (int k = 0; k < NUM_QUANTS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
      }
    }

    auto const psi = quant_vals[PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    out[ALPHA][i] = quant_vals[ALP];

    out[BETAX][i] = quant_vals[BETX];
    out[BETAY][i] = quant_vals[BETY];
    out[BETAZ][i] = quant_vals[BETZ];

    double g[3][3];
    g[0][0] = psi4;
    g[0][1] = 0.0;
    g[0][2] = 0.0;
    g[1][1] = psi4;
    g[1][2] = 0.0;
    g[2][2] = psi4;
    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];

    out[GXX][i] = g[0][0];
    out[GXY][i] = g[0][1];
    out[GXZ][i] = g[0][2];
    out[GYY][i] = g[1][1];
    out[GYZ][i] = g[1][2];
    out[GZZ][i] = g[2][2];

    out[KXX][i] = quant_vals[AXX] * psi4;
    out[KXY][i] = quant_vals[AXY] * psi4;
    out[KXZ][i] = quant_vals[AXZ] * psi4;
    out[KYY][i] = quant_vals[AYY] * psi4;
    out[KYZ][i] = quant_vals[AYZ] * psi4;
    out[KZZ][i] = quant_vals[AZZ] * psi4;

    double h = std::exp(quant_vals[H]);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at h=1
    if(h == 1.) {
      out[RHO][i]   = 0.;
      out[EPS][i]   = 0.;
      out[PRESS][i] = 0.;
    }
    else {
      if(eos_type == "Cold_Table") {
        using eos_t = FUKA_EOS_Wrapper<margherita_eos_t, margherita_1d>;

        out[RHO][i]   = EOS<eos_t, DENSITY>::get(h);
        out[EPS][i]   = EOS<eos_t, EPSILON>::get(h);
        out[PRESS][i] = EOS<eos_t, PRESSURE>::get(h);
      }

      if(eos_type == "Cold_PWPoly") {
        using eos_t = FUKA_EOS_Wrapper<margherita_eos_t, margherita_pwp>;

        out[RHO][i] = EOS<eos_t, DENSITY>::get(h);
        out[EPS][i] = EOS<eos_t, EPSILON>::get(h);
        out[PRESS][i] = EOS<eos_t, PRESSURE>::get(h);
      }
    }

    out[VELX][i] = quant_vals[UX];
    out[VELY][i] = quant_vals[UY];
    out[VELZ][i] = quant_vals[UZ];
  } // for i

  return out;
}
