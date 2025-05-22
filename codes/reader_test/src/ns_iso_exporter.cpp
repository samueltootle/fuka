#include "Configurator/config_bco.hpp"
#include "EOS/EOS.hh"
#include "Solvers/ns_isotropic/ns_isotropic_exporter.hpp"
#include "bco_utilities.hpp"
// #include "EOS/standalone/tov.hh"

// Kadath includes
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "kadath_adapted_polar.hpp"

#include "./reader_test_tools.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <thread>

using namespace Kadath;
using namespace Kadath::FUKA_Config;

using config_t = kadath_config_boost<BCO_NS_INFO>;
using reader_t = Kadath::FUKA_Solvers::CFMS_NS_ISO_Exporter;
constexpr unsigned int Npts = 1e4;
constexpr double range = 10;
constexpr double dx = range / Npts;

//template<class eos_t>
//void reader_test(config_t& bconfig) {
//    auto spacein = bconfig.space_filename();
//	FILE* ff1 = fopen (spacein.c_str(), "r") ;
//	Space_polar_adapted space (ff1) ;
//
//  // load the fields defined on the space
//	Scalar lap_Aterm(space, ff1) ;
//	Scalar nu       (space, ff1) ;
//    Scalar logh     (space, ff1) ;
//    Scalar lap_Bterm(space, ff1) ;
//	fclose(ff1) ;
//
//    auto nbrpts(space.get_domain(1)->get_nbr_points());
//    Index i(nbrpts);
//
//    auto x_v = space.get_domain(1)->get_cart(1);
//    auto z_v = space.get_domain(1)->get_cart(2);
//    Scalar x_s(space);
//    for(int d = 0; d < space.get_nbr_domains(); ++d)
//        x_s.set_domain(d) = space.get_domain(d)->get_cart(1);
//    x_s.std_base();
//    x_s.coef();
//    i.set(0) = nbrpts(0) - 1;
//    i.set(1) = nbrpts(1) - 1;
//    auto x = x_v(i);
//    auto z = z_v(i);
//    Point p(2);
//    p.set(1) = x;
//    p.set(2) = z;
//    cout << x_s.val_point(p) << endl;
//    // Point p(2);
//
//}

int main(int argc, char** argv) {
  // expecting a configuration file on execution
  if (argc < 2) {
    std::cerr << "Usage: ./reader /<path>/<ID base name>.info" << std::endl;
    std::cerr << "e.g. ./reader converged.NS.9.info" << endl;
    std::_Exit(EXIT_FAILURE);
  }

  // load the configuration
  std::string ifilename{argv[1]};

  // ID Coords
  std::vector<double> xx(Npts);
  std::vector<double> yy(Npts);
  std::vector<double> zz(Npts);

#pragma omp parallel for
  for (auto i = 0; i < Npts; ++i) {
    xx[i] += i * dx;
    yy[i] += i * dx;
    zz[i] += i * dx;
  }
  // END ID Coords

  config_t bconfig(ifilename);
  reader_t input_reader(ifilename);
  interp_data(input_reader, xx, yy, zz);

  //  kadath_config_boost<BCO_NS_INFO> bconfig(ifilename);
  //
  //  // setup the EOS
  //  const double h_cut = bconfig.eos<double>(HCUT);
  //  const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
  //  const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);
  //
  //  if(eos_type == "Cold_PWPoly") {
  //    using eos_t = Kadath::Margherita::Cold_PWPoly;
  //
  //    EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
  //
  //    // call reader to output diagnostics
  //    reader_test<eos_t>(bconfig);
  //  } else if(eos_type == "Cold_Table") {
  //    using eos_t = Kadath::Margherita::Cold_Table;
  //
  //    const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
//                            2000 : bconfig.eos<int>(INTERP_PTS);
  //
  //    EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
  //
  //    // call reader to output diagnostics
  //    reader_test<eos_t>(bconfig);
  //  } else {
  //    std::cerr << "Unknown EOSTYPE." << endl;
  //    std::_Exit(EXIT_FAILURE);
  //  }

  return EXIT_SUCCESS;
}
