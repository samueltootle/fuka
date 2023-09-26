#pragma once
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"
#include "EOS/EOS.hh"
#include "name_tools.hpp"
#include "exporter_utilities.hpp"
#include "Solvers/fuka_syst/fuka_syst_setup.hpp"
#include "Solvers/fuka_syst/fuka_syst_tools.hpp"
#include "codes_utilities.hpp"
#include "bco_utilities.hpp"
#include "exporter_utilities.hpp"

#include <cstdlib>
#include <string>
#include <filesystem>
#include <mutex>
#ifdef _OPENMP
  #include <omp.h>
#endif
namespace fs = std::filesystem;

namespace Kadath::FUKA_Solvers {
static std::mutex copy_mutex;
template<class space_t>
struct find_dom {
  int dom{-1};
  find_dom(std::unique_ptr<space_t>& space, Point& p) {
    auto ndom = space->get_nbr_domains();
    for(auto d = 2; d < ndom; ++d) {
      if(space->get_domain(d)->is_in(p)) {
        dom = d;
        break;
      }
    }
    if(dom == -1) {
      std::stringstream msg;
      msg << "Point " << p << " not found in the numerical space. ";
      msg << space.get() << ", " << bco_utils::get_radius(space->get_domain(2), INNER_BC) << endl;
      // cout << *(space->get_domain(d)) << endl;
      // throw std::runtime_error(msg.str().c_str());
      cout << msg.str() << endl;
    }
  }
  int operator()() { return dom; }
};

/**
 * @brief The following "Reader" is really a bandage until the FUKA_Solvers
 * can be rewritten.  Essentially many of the utilities her are duplicate to
 * the solvers, however, pointers are used instead of reference thereby allowing
 * the possibility for dynamic allocation for, e.g. multi-threaded importing of
 * the initial data.  It would have no impact on the solving of the system of equations
 * 
 * @tparam config_t 
 * @tparam space_t 
 */
template<class config_t, class space_t>
struct Reader {
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;

  ptr_data_member(space_t, space, unique);
  ptr_data_member(System_of_eqs, syst, unique);
  ptr_data_member(base_config_t, bconfig, unique);
  
  protected:
  int ndom{};

  public:
  Reader() : space(nullptr), syst(nullptr) {}
  Reader(std::string config_filename) : 
    space(nullptr), syst(nullptr), bconfig(nullptr) {
      bconfig.reset(new base_config_t{config_filename});
      bconfig->open_config();
  }
  
  virtual void load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space.reset(new space_t{ff1});
    syst.reset(new System_of_eqs{*space});
    fclose(ff1);
    ndom = space->get_nbr_domains();
  }
};

template<class config_t, class space_t>
struct CFMS_BH_Reader : public Reader<config_t, space_t> {
  enum XCTS_VARS : size_t {  
    XCTS_PSI,
    XCTS_ALPHA,
    XCTS_BETAX,
    XCTS_BETAY,
    XCTS_BETAZ,
    XCTS_AXX,
    XCTS_AXY,
    XCTS_AXZ,
    XCTS_AYY,
    XCTS_AYZ,
    XCTS_AZZ,
    NUM_XCTS_VARS 
  };

  enum OUTPUT_VARS : size_t {
    ALPHA,
    BETAX,
    BETAY,
    BETAZ,
    GXX,
    GXY,
    GXZ,
    GYY,
    GYZ,
    GZZ,
    KXX,
    KXY,
    KXZ,
    KYY,
    KYZ,
    KZZ,
    NUM_OUTPUT_VARS
  };
  using pointwise_ary_t = std::array<double, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  // Types
  using Reader<config_t, space_t>::base_space_t;
  using Reader<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Reader<config_t, space_t>::space;
  using Reader<config_t, space_t>::bconfig;
  using Reader<config_t, space_t>::syst;
  using Reader<config_t, space_t>::ndom;

  // CFMS_BH imported fields from file
  ptr_data_member(Scalar, conformal_factor, shared);
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Vector, shift, shared);

  // Constructed objects
  ptr_data_member(Base_tensor, basis, shared);
  ptr_data_member(Metric_flat, fmet, shared);
  ptr_data_member(Tensor, A, shared);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  bool export_ready{false};
  int const ndim{3};

  void load_solution_from_file() override {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space.reset(new space_t{ff1});
    conformal_factor.reset( new Scalar(*space.get(), ff1)) ;
    lapse.reset( new Scalar(*space.get(), ff1)) ;
    shift.reset( new Vector(*space.get(), ff1)) ;
    fclose(ff1);
    basis.reset(new Base_tensor{shift->get_basis()});
    fmet.reset(new Metric_flat(*space, *basis));
    ndom = space->get_nbr_domains();
  }

  void extract_xcts_grid_functions() {
    if(quants.capacity() != XCTS_VARS::NUM_XCTS_VARS) {
      for (int i = 0; i < XCTS_VARS::NUM_XCTS_VARS; ++i)
        quants.push_back(std::cref(*conformal_factor));
    }
    syst.reset(new System_of_eqs (*space));    
    fmet->set_system(*syst, "f") ;

    // Fields - must be initialized before common setup
    syst->add_cst("N"  , *lapse) ;
    syst->add_cst("bet", *shift) ;

    // syst->add_def("Ts^i = N * bet^i");
    syst->add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    A.reset(new Tensor(syst->give_val_def("A")));

    quants[XCTS_VARS::XCTS_PSI] = std::cref(*conformal_factor);
    quants[XCTS_VARS::XCTS_ALPHA] = std::cref(*lapse);
    quants[XCTS_VARS::XCTS_BETAX] = std::cref((*shift)(1));
    quants[XCTS_VARS::XCTS_BETAY] = std::cref((*shift)(2));
    quants[XCTS_VARS::XCTS_BETAZ] = std::cref((*shift)(3));

    export_utils::add_tensor_refs(quants, {
      XCTS_VARS::XCTS_AXX, 
      XCTS_VARS::XCTS_AXY, 
      XCTS_VARS::XCTS_AXZ, 
      XCTS_VARS::XCTS_AYY, 
      XCTS_VARS::XCTS_AYZ, 
      XCTS_VARS::XCTS_AZZ}, *A);
    export_ready = true;
  }

  public:
  bool is_export_ready() const { return export_ready; }
  CFMS_BH_Reader() : Reader<config_t, space_t>(),
    basis(nullptr), fmet(nullptr),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr) {}
  CFMS_BH_Reader(std::string config_filename) :
    Reader<config_t, space_t>(config_filename),
      basis(nullptr), fmet(nullptr),
        conformal_factor(nullptr), lapse(nullptr), shift(nullptr) {
  
    load_solution_from_file();
    extract_xcts_grid_functions();
  }

  CFMS_BH_Reader(CFMS_BH_Reader const & r) : fmet(nullptr) {
    copy_mutex.lock();
    space.reset(new space_t((*r.get_space())));
    basis.reset(new Base_tensor(*space, r.get_basis()->get_basis(0)));
    fmet.reset(new Metric_flat(*space, *basis));
    bconfig.reset(new config_t(*r.bconfig));
    
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    bconfig->set_filename("test");
    export_ready = false;

    extract_xcts_grid_functions();
    copy_mutex.unlock();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
    // std::cout << "copy\n";
  }

  CFMS_BH_Reader(CFMS_BH_Reader&& b) noexcept = delete;
  CFMS_BH_Reader& operator=(const CFMS_BH_Reader& b)
  {
    if (this == &b) return *this;

    CFMS_BH_Reader tmp(b);
    *this = std::move(tmp);
// std::cout << "assignment\n";
    return *this;
  }

  public:

  std::vector<double> interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
    
    std::vector<double> quant_vals(XCTS_VARS::NUM_XCTS_VARS);

    // Initial guess of the excision radius - needed for filling
    double rbh = bco_utils::get_radius(space->get_domain(2), INNER_BC);

    double r2yz = y * y + z * z;

    double r = std::sqrt(x * x + r2yz);

    // lambda function for filling excised region
    auto interp_f = [&](auto& ah_r, auto& extrap_r, auto bh_ori) {
      // Avoid division by "0"
      extrap_r = (extrap_r <= 1e-12) ? 1e-12 : extrap_r;
      double x_shifted = (x == 0.) ? 1e-14 : x;
      double theta = std::acos(z / extrap_r);
      
      // atan2 is needed here
      double phi = std::atan2(y, (x_shifted - bh_ori)); 

      // Where the filling takes places
      export_utils::spherical_turduck(
        quants, quant_vals, interp_order, delta_r_rel, interpolation_offset, 
        rbh, extrap_r, theta, phi, 2, bh_ori
      );
    };

    if (r <= (1. + interpolation_offset) * rbh) {
      interp_f(rbh, r, 0.);
      quant_vals[XCTS_VARS::XCTS_ALPHA] = -1;
    } else { 
      Point abs_coords(ndim);
      abs_coords.set(1) = x;
      abs_coords.set(2) = y;
      abs_coords.set(3) = z;
      find_dom fd(space, abs_coords);
      quant_vals[XCTS_VARS::XCTS_ALPHA] = fd();
      // for (int k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
      //   quant_vals[k] = quants[k].get().val_point(abs_coords);
      // }
    }
    return quant_vals;
  }

  std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS> export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
    std::array<std::vector<double>,OUTPUT_VARS::NUM_OUTPUT_VARS> out;
    for(auto& v : out)
      v.resize(npoints);
    
    for (int i = 0; i < npoints; ++i) {
      
      auto quant_vals = interpolate_pointwise(xx[i], yy[i], zz[i], interpolation_offset, interp_order, delta_r_rel);
      
      // Fill output vector by storing non-conformal quantities
      auto const psi = quant_vals[XCTS_VARS::XCTS_PSI];
      auto const psi2 = psi * psi;
      auto const psi4 = psi2 * psi2;

      out[OUTPUT_VARS::ALPHA][i] = quant_vals[XCTS_VARS::XCTS_ALPHA];

      out[OUTPUT_VARS::BETAX][i] = quant_vals[XCTS_VARS::XCTS_BETAX];
      out[OUTPUT_VARS::BETAY][i] = quant_vals[XCTS_VARS::XCTS_BETAY];
      out[OUTPUT_VARS::BETAZ][i] = quant_vals[XCTS_VARS::XCTS_BETAZ];

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

      out[OUTPUT_VARS::GXX][i] = g[0][0];
      out[OUTPUT_VARS::GXY][i] = g[0][1];
      out[OUTPUT_VARS::GXZ][i] = g[0][2];
      out[OUTPUT_VARS::GYY][i] = g[1][1];
      out[OUTPUT_VARS::GYZ][i] = g[1][2];
      out[OUTPUT_VARS::GZZ][i] = g[2][2];

      out[OUTPUT_VARS::KXX][i] = quant_vals[XCTS_VARS::XCTS_AXX] * psi4;
      out[OUTPUT_VARS::KXY][i] = quant_vals[XCTS_VARS::XCTS_AXY] * psi4;
      out[OUTPUT_VARS::KXZ][i] = quant_vals[XCTS_VARS::XCTS_AXZ] * psi4;
      out[OUTPUT_VARS::KYY][i] = quant_vals[XCTS_VARS::XCTS_AYY] * psi4;
      out[OUTPUT_VARS::KYZ][i] = quant_vals[XCTS_VARS::XCTS_AYZ] * psi4;
      out[OUTPUT_VARS::KZZ][i] = quant_vals[XCTS_VARS::XCTS_AZZ] * psi4;
    }
  }
  
  std::array<double, OUTPUT_VARS::NUM_OUTPUT_VARS> export_pointwise(
    double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
    
    std::array<double, OUTPUT_VARS::NUM_OUTPUT_VARS> out;
      
    auto quant_vals = interpolate_pointwise(x, y, z, interpolation_offset, interp_order, delta_r_rel);
    
    // Fill output vector by storing non-conformal quantities
    auto const psi = quant_vals[XCTS_VARS::XCTS_PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    out[OUTPUT_VARS::ALPHA] = quant_vals[XCTS_VARS::XCTS_ALPHA];

    out[OUTPUT_VARS::BETAX] = quant_vals[XCTS_VARS::XCTS_BETAX];
    out[OUTPUT_VARS::BETAY] = quant_vals[XCTS_VARS::XCTS_BETAY];
    out[OUTPUT_VARS::BETAZ] = quant_vals[XCTS_VARS::XCTS_BETAZ];

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

    out[OUTPUT_VARS::GXX] = g[0][0];
    out[OUTPUT_VARS::GXY] = g[0][1];
    out[OUTPUT_VARS::GXZ] = g[0][2];
    out[OUTPUT_VARS::GYY] = g[1][1];
    out[OUTPUT_VARS::GYZ] = g[1][2];
    out[OUTPUT_VARS::GZZ] = g[2][2];

    out[OUTPUT_VARS::KXX] = quant_vals[XCTS_VARS::XCTS_AXX] * psi4;
    out[OUTPUT_VARS::KXY] = quant_vals[XCTS_VARS::XCTS_AXY] * psi4;
    out[OUTPUT_VARS::KXZ] = quant_vals[XCTS_VARS::XCTS_AXZ] * psi4;
    out[OUTPUT_VARS::KYY] = quant_vals[XCTS_VARS::XCTS_AYY] * psi4;
    out[OUTPUT_VARS::KYZ] = quant_vals[XCTS_VARS::XCTS_AYZ] * psi4;
    out[OUTPUT_VARS::KZZ] = quant_vals[XCTS_VARS::XCTS_AZZ] * psi4;
    return out;
  }
};
}