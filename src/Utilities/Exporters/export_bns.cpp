#include <kadath_bin_ns.hpp>
#include <EOS/EOS.hh>
#include <coord_fields.hpp>
#include <Configurator/config_binary.hpp>
#include <exporter_utilities.hpp>
#include <bco_utilities.hpp>

#include <cmath>
#include "EOS/FUKA_EOS_Utilities.hh"

using namespace Kadath;
using namespace export_utils;
using namespace Kadath::FUKA_Config;

std::array<std::vector<double>,NUM_OUT> KadathExportBNS(int const npoints,
                                                        double const * xx, double const * yy, double const * zz,
                                                        char const * fn) {
  std::string filename(fn);

  kadath_config_boost<BIN_INFO> bconfig(filename);

  using namespace Kadath::FUKA_EOS;
  EOS_initialize::init(bconfig, NODES::BCO1);
  // get const EOS information
  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE, BCO1);

  double& units = bconfig(QPIG);
  double& distance = bconfig(DIST);
  double& omega = bconfig(GOMEGA);
  double& ome1 = bconfig(OMEGA, BCO1);
  double& ome2 = bconfig(OMEGA, BCO2);
  double& Madm1 = bconfig(MADM, BCO1);
  double& Madm2 = bconfig(MADM, BCO2);
  double& axis = bconfig(COM);

  /* file containing KADATH fields must have same name as config file
   * with only the extension being different */
  std::string kadath_filename = bconfig.space_filename();

	FILE* fin = fopen(kadath_filename.c_str(), "r") ;
	Space_bin_ns space(fin) ;
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

  int ndom = space.get_nbr_domains();

 	double xc1 = Kadath::bco_utils::get_center(space,space.NS1);
 	double xc2 = Kadath::bco_utils::get_center(space,space.NS2);
  double xo  = Kadath::bco_utils::get_center(space,ndom-1);

  Metric_flat fmet(space, basis);

  // init coordinate fields
  CoordFields<Space_bin_ns> cfields(space);
  vec_ary_t coord_vectors {default_binary_vector_ary(space)};
  update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

  System_of_eqs syst(space, 0, ndom - 1);

  fmet.set_system(syst, "f");

  Param p;
  EOS_Function_Dispatcher::dispatch<set_eos_ope_struct>(bconfig, eos_type, syst, p);

  syst.add_cst("4piG", units);
  syst.add_cst("PI", M_PI);

  syst.add_cst("omes1", ome1);
  syst.add_cst("omes2", ome2);

  syst.add_cst("mg"   , *coord_vectors[GLOBAL_ROT]) ;
  syst.add_cst("mm"   , *coord_vectors[BCO1_ROT]) ;
  syst.add_cst("mp"   , *coord_vectors[BCO2_ROT]) ;

  syst.add_cst("ex"   , *coord_vectors[EX])  ;
  syst.add_cst("ey"   , *coord_vectors[EY])  ;
  syst.add_cst("ez"   , *coord_vectors[EZ])  ;

  syst.add_cst("sm"   , *coord_vectors[S_BCO1])  ;
  syst.add_cst("sp"   , *coord_vectors[S_BCO2])  ;
  syst.add_cst("einf" , *coord_vectors[S_INF])  ;

  syst.add_cst("xaxis", axis);
  syst.add_cst("ome", omega);

  syst.add_cst("P", conf);
  syst.add_cst("N", lapse);
  syst.add_cst("bet", shift);
  syst.add_cst("phi", phi);

  syst.add_cst("H", logh);

  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");

  syst.add_def("Morb^i = mg^i + xaxis * ey^i");
  syst.add_def("omega^i = bet^i + ome * Morb^i");

  for(int d = space.NS1; d <= space.ADAPTED1; ++d){
    syst.add_def  (d, "s^i  = omes1 * mm^i");
  }
  for(int d = space.NS2; d <= space.ADAPTED2; ++d){
    syst.add_def  (d, "s^i  = omes2 * mp^i");
  }

  syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");

  syst.add_def("h = exp(H)");

  for(int d = 0; d < ndom; ++d) {
    if((d <= space.ADAPTED1) || (d <= space.ADAPTED2 && d >= space.NS2))
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

  std::array<std::vector<double>,NUM_OUT> out;
  for(auto& v : out)
    v.resize(npoints);

  for (int i = 0; i < npoints; ++i) {
    std::vector<double> quant_vals(NUM_QUANTS);

    // get number of dimensions and construct index
    int ndim = 3;

    // construct Kadath point, shifted wrt the COM
    Point abs_coords(ndim);
    abs_coords.set(1) = xx[i] - axis;
    abs_coords.set(2) = yy[i];
    abs_coords.set(3) = zz[i];

    for (int k = 0; k < NUM_QUANTS; ++k) {
      quant_vals[k] = quants[k].get().val_point(abs_coords);
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

    // get quantities point-wise, since h is smoothest, and cut data at h=1
    if(h == 1.) {
      out[RHO][i] = 0.;
      out[EPS][i] = 0.;
      out[PRESS][i] = 0.;
    }
    else {
      if(eos_type == "Cold_Table") {
        using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, margherita_1d>;

        out[RHO][i]   = EOS<eos_t, DENSITY>::get(h);
        out[EPS][i]   = EOS<eos_t, EPSILON>::get(h);
        out[PRESS][i] = EOS<eos_t, PRESSURE>::get(h);
      }

      if(eos_type == "Cold_PWPoly") {
        using eos_t = FUKA_EOS_Wrapper<fuka_eos_t, margherita_pwp>;

        out[RHO][i]   = EOS<eos_t, DENSITY>::get(h);
        out[EPS][i]   = EOS<eos_t, EPSILON>::get(h);
        out[PRESS][i] = EOS<eos_t, PRESSURE>::get(h);
      }
    }

    out[VELX][i] = quant_vals[UX];
    out[VELY][i] = quant_vals[UY];
    out[VELZ][i] = quant_vals[UZ];
  } // for i

  return out;
}
