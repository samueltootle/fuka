#pragma once
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "codes_utilities.hpp"
#include "bco_utilities.hpp"
#include "adapted_bh.hpp"
#include "exporter_utilities.hpp"
#include "name_tools.hpp"

#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"

#include "EOS/EOS.hh"

#include "Solvers/fuka_syst/fuka_syst_setup.hpp"
#include "Solvers/fuka_syst/fuka_syst_tools.hpp"


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
  ptr_data_member(base_config_t, bconfig, unique);
  
  protected:
  int ndom{};

  public:
  Reader() : space(nullptr) {}
  Reader(std::string config_filename) : 
    space(nullptr), bconfig(nullptr) {
      bconfig.reset(new base_config_t{config_filename});
      bconfig->open_config();
  }
  
  virtual void load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space.reset(new space_t{ff1});
    fclose(ff1);
    ndom = space->get_nbr_domains();
  }
};

struct CFMS_BH_Reader : public Reader<Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_BH_INFO>, Space_adapted_bh> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_BH_INFO>;
  using space_t = Space_adapted_bh;

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

  using interp_ary_t = std::array<double, NUM_XCTS_VARS>; 
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>; 
  using grid_ary_t = std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  // Types
  using Reader<config_t, space_t>::base_space_t;
  using Reader<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Reader<config_t, space_t>::space;
  using Reader<config_t, space_t>::bconfig;
  using Reader<config_t, space_t>::ndom;

  // CFMS_BH imported fields from file
  ptr_data_member(Scalar, conformal_factor, unique);
  ptr_data_member(Scalar, lapse, unique);
  ptr_data_member(Vector, shift, unique);

  // Constructed objects
  ptr_data_member(Tensor, A, unique);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  interp_ary_t quant_vals;
  output_ary_t out_pw;
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

    ndom = space->get_nbr_domains();
  }

  void extract_computed_grid_functions() {
    System_of_eqs syst(*space);
    Base_tensor basis(shift->get_basis());
    Metric_flat fmet(*space, basis);
    fmet.set_system(syst, "f") ;

    // Fields - must be initialized before common setup
    syst.add_cst("N"  , *lapse) ;
    syst.add_cst("bet", *shift) ;

    syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    A.reset(new Tensor(syst.give_val_def("A")));
    A->coef();
  }

  void populate_quants() {
    if(quants.capacity() != XCTS_VARS::NUM_XCTS_VARS) {
      for (size_t i = 0; i < XCTS_VARS::NUM_XCTS_VARS; ++i)
        quants.push_back(std::cref(*conformal_factor));
    }
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
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  bool is_export_ready() const { return export_ready; }  
  const int & get_ndim() const { return ndim; }

  CFMS_BH_Reader() : Reader<config_t, space_t>(),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr) { }
  
  CFMS_BH_Reader(std::string config_filename) :
    Reader<config_t, space_t>(config_filename),
        conformal_factor(nullptr), lapse(nullptr), shift(nullptr) {

    load_solution_from_file();
    extract_computed_grid_functions();
    populate_quants();
    this->export_pointwise(0.5, 0., 0.);
  }

  CFMS_BH_Reader(CFMS_BH_Reader const & r) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;

    space.reset(new space_t((*r.get_space())));
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    A.reset(new Tensor(*space, *r.A.get()));

    bconfig.reset(new config_t(*r.bconfig));
    
    populate_quants();
  }

  CFMS_BH_Reader(CFMS_BH_Reader&& b) noexcept = delete;
  CFMS_BH_Reader& operator=(const CFMS_BH_Reader& b)
  {
    if (this == &b) return *this;

    CFMS_BH_Reader tmp(b);
    *this = std::move(tmp);
    return *this;
  }

  public:

  interp_ary_t interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {

    // Initial guess of the excision radius - needed for filling
    // FIXME make excision generic
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
    } else { 
      Point abs_coords(ndim);
      abs_coords.set(1) = x;
      abs_coords.set(2) = y;
      abs_coords.set(3) = z;

      for (size_t k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
      }
    }
    return quant_vals;
  }
  
  output_ary_t export_pointwise(
    double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
      
    quant_vals = interpolate_pointwise(x, y, z, interpolation_offset, interp_order, delta_r_rel);
    
    // Fill output vector by storing non-conformal quantities
    auto const psi = quant_vals[XCTS_VARS::XCTS_PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    out_pw[OUTPUT_VARS::ALPHA] = quant_vals[XCTS_VARS::XCTS_ALPHA];

    out_pw[OUTPUT_VARS::BETAX] = quant_vals[XCTS_VARS::XCTS_BETAX];
    out_pw[OUTPUT_VARS::BETAY] = quant_vals[XCTS_VARS::XCTS_BETAY];
    out_pw[OUTPUT_VARS::BETAZ] = quant_vals[XCTS_VARS::XCTS_BETAZ];

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

    out_pw[OUTPUT_VARS::GXX] = g[0][0];
    out_pw[OUTPUT_VARS::GXY] = g[0][1];
    out_pw[OUTPUT_VARS::GXZ] = g[0][2];
    out_pw[OUTPUT_VARS::GYY] = g[1][1];
    out_pw[OUTPUT_VARS::GYZ] = g[1][2];
    out_pw[OUTPUT_VARS::GZZ] = g[2][2];

    out_pw[OUTPUT_VARS::KXX] = quant_vals[XCTS_VARS::XCTS_AXX] * psi4;
    out_pw[OUTPUT_VARS::KXY] = quant_vals[XCTS_VARS::XCTS_AXY] * psi4;
    out_pw[OUTPUT_VARS::KXZ] = quant_vals[XCTS_VARS::XCTS_AXZ] * psi4;
    out_pw[OUTPUT_VARS::KYY] = quant_vals[XCTS_VARS::XCTS_AYY] * psi4;
    out_pw[OUTPUT_VARS::KYZ] = quant_vals[XCTS_VARS::XCTS_AYZ] * psi4;
    out_pw[OUTPUT_VARS::KZZ] = quant_vals[XCTS_VARS::XCTS_AZZ] * psi4;
    return out_pw;
  }

  grid_ary_t export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
    
    grid_ary_t out;
    for(auto& v : out) {
      v.resize(npoints);
    }
    
    for (int i = 0; i < npoints; ++i) {
      export_pointwise(xx[i], yy[i], zz[i], interpolation_offset, interp_order, delta_r_rel);
      
      out[OUTPUT_VARS::ALPHA][i] = out_pw[OUTPUT_VARS::ALPHA];

      out[OUTPUT_VARS::BETAX][i] = out_pw[OUTPUT_VARS::BETAX];
      out[OUTPUT_VARS::BETAY][i] = out_pw[OUTPUT_VARS::BETAY];
      out[OUTPUT_VARS::BETAZ][i] = out_pw[OUTPUT_VARS::BETAZ];

      out[OUTPUT_VARS::GXX][i] = out_pw[OUTPUT_VARS::GXX];
      out[OUTPUT_VARS::GXY][i] = out_pw[OUTPUT_VARS::GXY];
      out[OUTPUT_VARS::GXZ][i] = out_pw[OUTPUT_VARS::GXZ];
      out[OUTPUT_VARS::GYY][i] = out_pw[OUTPUT_VARS::GYY];
      out[OUTPUT_VARS::GYZ][i] = out_pw[OUTPUT_VARS::GYZ];
      out[OUTPUT_VARS::GZZ][i] = out_pw[OUTPUT_VARS::GZZ];

      out[OUTPUT_VARS::KXX][i] = out_pw[OUTPUT_VARS::KXX];
      out[OUTPUT_VARS::KXY][i] = out_pw[OUTPUT_VARS::KXY];
      out[OUTPUT_VARS::KXZ][i] = out_pw[OUTPUT_VARS::KXZ];
      out[OUTPUT_VARS::KYY][i] = out_pw[OUTPUT_VARS::KYY];
      out[OUTPUT_VARS::KYZ][i] = out_pw[OUTPUT_VARS::KYZ];
      out[OUTPUT_VARS::KZZ][i] = out_pw[OUTPUT_VARS::KZZ];
    }
    return out;
  }
};
}
