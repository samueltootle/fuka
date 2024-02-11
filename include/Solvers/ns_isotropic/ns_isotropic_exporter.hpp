#include "Solvers/exporter.hpp"
namespace Kadath::FUKA_Solvers {

struct CFMS_NS_ISO_Exporter : public Exporter<Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_ISO_NS_INFO>, Space_polar_adapted> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_ISO_NS_INFO>;
  using space_t = Space_polar_adapted;

  enum OUTPUT_BASIS : size_t {
    CARTESIAN,
    SPHERICAL,
    UNDEFINED,
  };
  
  // Input ID grid functions that we interpolate on
  enum ISO_VARS : size_t {
    ISO_ALPHA,
    ISO_METRIC_A,
    ISO_METRIC_B,
    ISO_METRIC_OMEGA,
    ISO_DOMEGA_DR,
    ISO_DOMEGA_DTHETA,
    ISO_H,
    ISO_U,
    NUM_ISO_VARS 
  };

  enum CART_VARS : size_t {
    CART_ALPHA,
    CART_G11,
    CART_G12,
    CART_G13,
    CART_G22,
    CART_G23,
    CART_G33,
    CART_K11,
    CART_K12,
    CART_K13,
    CART_K22,
    CART_K23,
    CART_K33,
    CART_VEL1,
    CART_VEL2,
    CART_VEL3,
    CART_H,
    NUM_CART_VARS 
  };

  // Output physical grid functions
  enum OUTPUT_VARS : size_t {
    ALPHA,
    BETA1,
    BETA2,
    BETA3,
    G11,
    G12,
    G13,
    G22,
    G23,
    G33,
    K11,
    K12,
    K13,
    K22,
    K23,
    K33,
    RHO,
    EPS,
    PRESS,
    VEL1,
    VEL2,
    VEL3,
    NUM_OUTPUT_VARS
  };

  std::map<std::string, OUTPUT_VARS> output_var_map {
    {"lapse", OUTPUT_VARS::ALPHA},
    {"beta1", OUTPUT_VARS::BETA1},
    {"beta2", OUTPUT_VARS::BETA2},
    {"beta3", OUTPUT_VARS::BETA3},
    {"g11"  , OUTPUT_VARS::G11},
    {"g12"  , OUTPUT_VARS::G12},
    {"g13"  , OUTPUT_VARS::G13},
    {"g22"  , OUTPUT_VARS::G22},
    {"g23"  , OUTPUT_VARS::G23},
    {"g33"  , OUTPUT_VARS::G33},
    {"k11"  , OUTPUT_VARS::K11},
    {"k12"  , OUTPUT_VARS::K12},
    {"k13"  , OUTPUT_VARS::K13},
    {"k22"  , OUTPUT_VARS::K22},
    {"k23"  , OUTPUT_VARS::K23},
    {"k33"  , OUTPUT_VARS::K33},
    {"rho"  , OUTPUT_VARS::RHO},
    {"eps"  , OUTPUT_VARS::EPS},
    {"press", OUTPUT_VARS::PRESS},
    {"vel1" , OUTPUT_VARS::VEL1},
    {"vel2" , OUTPUT_VARS::VEL2},
    {"vel3" , OUTPUT_VARS::VEL3},
  };

  using interp_ary_t = std::array<double, NUM_ISO_VARS>; 
  using interp_cart_ary_t = std::array<double, NUM_CART_VARS>; 
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>; 
  using grid_ary_t = std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  std::vector<ISO_VARS> spacetime_indicies__iso{
    ISO_VARS::ISO_ALPHA,
    ISO_VARS::ISO_METRIC_A,
    ISO_VARS::ISO_METRIC_B,
    ISO_VARS::ISO_METRIC_OMEGA,
    ISO_VARS::ISO_DOMEGA_DR,
    ISO_VARS::ISO_DOMEGA_DTHETA,
  };
  std::vector<ISO_VARS> fluid_indicies__iso{
    ISO_VARS::ISO_H,
    ISO_VARS::ISO_U
  };
  std::vector<CART_VARS> spacetime_indicies__cart{
    CART_VARS::CART_ALPHA,
    CART_VARS::CART_G11,
    CART_VARS::CART_G12,
    CART_VARS::CART_G13,
    CART_VARS::CART_G22,
    CART_VARS::CART_G23,
    CART_VARS::CART_G33,
    CART_VARS::CART_K11,
    CART_VARS::CART_K12,
    CART_VARS::CART_K13,
    CART_VARS::CART_K22,
    CART_VARS::CART_K23,
    CART_VARS::CART_K33,
  };
  std::vector<CART_VARS> fluid_indicies__cart{
    CART_VARS::CART_H,
    CART_VARS::CART_VEL1,
    CART_VARS::CART_VEL2,
    CART_VARS::CART_VEL3,
  };  

  // Types from Base class
  using Exporter<config_t, space_t>::base_space_t;
  using Exporter<config_t, space_t>::base_config_t;

  // Class members from base class to avoid using this->
  using Exporter<config_t, space_t>::space;
  using Exporter<config_t, space_t>::bconfig;
  using Exporter<config_t, space_t>::ndom;

  // CFMS_NS imported fields from file
  // In principle these can be reset to nullptr
  // after the Spherical basis quantities are
  // computed.
  ptr_data_member(Scalar, lap_Aterm, shared);
  ptr_data_member(Scalar, lap_Bterm, shared);
  ptr_data_member(Scalar, Nu, shared);
  ptr_data_member(Scalar, lap_omega_term, shared);
  ptr_data_member(Scalar, logh, shared);
  ptr_data_member(Scalar, omega, shared);

  // Constructed objects - Spherical Basis
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Scalar, metric_A, shared);
  ptr_data_member(Scalar, metric_B, shared);  
  ptr_data_member(Scalar, metric_omega, shared);
  ptr_data_member(Scalar, domega_dr, shared);
  ptr_data_member(Scalar, domega_dt, shared);
  ptr_data_member(Scalar, fluidvel, shared);

  // Constructed objects - Cartesian Basis
  // see basis_transform()
  ptr_data_member(Scalar, gxx, shared);
  ptr_data_member(Scalar, gxy, shared);
  ptr_data_member(Scalar, gxz, shared);
  ptr_data_member(Scalar, gyy, shared);
  ptr_data_member(Scalar, gyz, shared);
  ptr_data_member(Scalar, gzz, shared);
  ptr_data_member(Scalar, Kxx, shared);
  ptr_data_member(Scalar, Kxy, shared);
  ptr_data_member(Scalar, Kxz, shared);
  ptr_data_member(Scalar, Kyy, shared);
  ptr_data_member(Scalar, Kyz, shared);
  ptr_data_member(Scalar, Kzz, shared);
  ptr_data_member(Scalar, velUx, shared);
  ptr_data_member(Scalar, velUy, shared);
  ptr_data_member(Scalar, velUz, shared);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  std::vector<std::reference_wrapper<const Scalar>> output_quants;
  interp_cart_ary_t quant_vals;
  output_ary_t out_pw;
  bool export_ready{false};
  int const ndim{2};
  OUTPUT_BASIS output_base{OUTPUT_BASIS::UNDEFINED};
  
  // EOS Parameters
  double h_cut{0};
  std::string eos_file{};
  std::string eos_type{};
  
  /**
   * @brief We initialize the FUKA EOS module here to avoid needing to handle this by the
   * an interface code on the importer side since the FUKA EOS module is independent of the
   * EOS module used in the evolution code. 
   */
  void initialize_eos();

  /**
   * @brief Here we load the ID solution from file in order to initialize the ID fields
   * and Exporter variables
   * 
   */
  void load_solution_from_file() override;
  
  /**
   * @brief Setup the EOS operators in System_of_eqs
   * 
   * @tparam eos_t C++ Polytrope/Table type
   * @param syst System of equations
   * @param p Parameter container (unused, but required by add_ope)
   */
  template<class eos_t>
  void set_eos_ope(System_of_eqs& syst, Param& p) {
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  }

  /**
   * @brief Compute needed grid functions (gf) within Kadath's System_of_eqs before
   * populating the gf pointers
   * 
   */
  void extract_computed_grid_functions();

  /**
   * @brief Using the populated struct pointers, populate the "quants" array with the
   * associated gf references
   * 
   */
  void populate_quants();

  void initialize_Cartesian_fields();

    // Initialize 
    

  public:
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  bool is_export_ready() const { return export_ready; }  
  const int & get_ndim() const { return ndim; }
  
  CFMS_NS_ISO_Exporter() : Exporter<config_t, space_t>(),
    lap_Aterm(nullptr), lap_Bterm(nullptr), Nu(nullptr), lap_omega_term(nullptr), logh(nullptr), 
      metric_A(nullptr), metric_B(nullptr), lapse(nullptr), metric_omega(nullptr), 
        domega_dr(nullptr), domega_dt(nullptr), fluidvel(nullptr), gxx(nullptr),
          gxy(nullptr), gxz(nullptr), gyy(nullptr), gyz(nullptr), gzz(nullptr), Kxx(nullptr),
            Kxy(nullptr), Kxz(nullptr), Kyy(nullptr), Kyz(nullptr), Kzz(nullptr),
            velUx(nullptr), velUy(nullptr), velUz(nullptr) {}
  
  CFMS_NS_ISO_Exporter(std::string config_filename) :
    Exporter<config_t, space_t>(config_filename),
    lap_Aterm(nullptr), lap_Bterm(nullptr), Nu(nullptr), lap_omega_term(nullptr), logh(nullptr), 
      metric_A(nullptr), metric_B(nullptr), lapse(nullptr), metric_omega(nullptr), 
        domega_dr(nullptr), domega_dt(nullptr), fluidvel(nullptr),  gxx(nullptr),
          gxy(nullptr), gxz(nullptr), gyy(nullptr), gyz(nullptr), gzz(nullptr), Kxx(nullptr),
            Kxy(nullptr), Kxz(nullptr), Kyy(nullptr), Kyz(nullptr), Kzz(nullptr),
            velUx(nullptr), velUy(nullptr), velUz(nullptr) {

    load_solution_from_file();
    initialize_Cartesian_fields();
    initialize_eos();
    extract_computed_grid_functions();
    basis_transform__Spherical_to_Cart();
    populate_quants();
    
    // This is to avoid a "bug" where "something" in kadath is not
    // correctly initialized prior to copying to other threads resulting
    // in undefined behavior.  By running the interpolator once, this
    // bug seems to be avoided.
    this->export_pointwise(0.5, 0., 0.);
  }

  /**
   * @brief QI solutions are computed in a spherical basis.  Here we transform
   * to a Cartesian basis prior to interpolation since most code
   * 
   */
  void basis_transform__Spherical_to_Cart();

  /**
   * @brief Construct a new object - a mutex is used to allow for thread safety in
   * copying.
   *  
   */
  CFMS_NS_ISO_Exporter(CFMS_NS_ISO_Exporter const & r);
  CFMS_NS_ISO_Exporter(CFMS_NS_ISO_Exporter&& b) noexcept = delete;
  CFMS_NS_ISO_Exporter& operator=(const CFMS_NS_ISO_Exporter& b);

  public:

  /**
   * @brief For a given coordinate, interpolate the ID solution
   * 
   * @param x 
   * @param y 
   * @param z 
   * @return interp_ary_t 
   */
  interp_cart_ary_t interpolate_pointwise(double const & x, double const & y, double const & z);

    /**
   * @brief For a given coordinate, interpolate the ID solution in a spherical basis
   * 
   * @param x 
   * @param y 
   * @param z 
   * @return interp_ary_t 
   */
  interp_ary_t interpolate_pointwise__spherical(double const & x, double const & y, double const & z);
  
  /**
   * @brief For a given coordinate, interpolate the ID solution for a subset of the ID variables
   * 
   * @param x 
   * @param y 
   * @param z 
   * @param slice indicies of array subset to interpolate
   * @return interp_ary_t 
   */
  interp_cart_ary_t interpolate_pointwise_subset(double const & x, double const & y, double const & z,
    std::vector<CART_VARS> slice);
  
  /**
   * @brief Export: Load only the spacetime variables in Cartesian coordinates
   * 
   * @param x 
   * @param y 
   * @param z
   * @return output_ary_t 
   */
  output_ary_t export_pointwise_spacetime_vars(double const & x, double const & y, double const & z);

  /**
   * @brief Export: Load only the fluid variables
   * 
   * @tparam eos_t EOS type
   * @param x 
   * @param y 
   * @param z
   * @return output_ary_t 
   */
  template<class eos_t>
  output_ary_t export_pointwise_fluid_vars_imp(double const & x, double const & y, double const & z) {
    
    // Reset to NAN
    quant_vals.fill(NAN);
    quant_vals = interpolate_pointwise_subset(x, y, z, fluid_indicies__cart);

    double const H = quant_vals[CART_VARS::CART_H];
    double h = std::exp(H);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-14) {
      rho = 0.;
      eps = 0.;
      press = 0.;
    }
    else {
      rho = EOS<eos_t, DENSITY>::get(h);
      eps = EOS<eos_t, EPSILON>::get(h);
      press = EOS<eos_t, PRESSURE>::get(h);
    }
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
    out_pw[OUTPUT_VARS::VEL1] = quant_vals[CART_VARS::CART_VEL1];
    out_pw[OUTPUT_VARS::VEL2] = quant_vals[CART_VARS::CART_VEL2];
    out_pw[OUTPUT_VARS::VEL3] = quant_vals[CART_VARS::CART_VEL3];
    output_base = OUTPUT_BASIS::CARTESIAN;
    return out_pw;
  }

  /**
   * @brief Export: Interface for only loading the fluid variables
   * 
   * @tparam eos_t EOS type
   * @param x 
   * @param y 
   * @param z
   * @return output_ary_t 
   */
  output_ary_t export_pointwise_fluid_vars(double const & x, double const & y, double const & z);

  /**
   * @brief Export an array of OUTPUT_VARS for a given point
   * 
   * @tparam eos_t C++ Polytrope/Table type
   * @param x 
   * @param y 
   * @param z 
   * @return output_ary_t Interpolated solution at x,y,z
   */
  template<class eos_t>
  output_ary_t export_pointwise_imp(double const & x, double const & y, double const & z) {        

    // Reset to NAN
    quant_vals.fill(NAN);
    quant_vals = interpolate_pointwise(x, y, z);

    out_pw[OUTPUT_VARS::ALPHA] = quant_vals[CART_VARS::CART_ALPHA];

    out_pw[OUTPUT_VARS::BETA1] = 0.0;
    out_pw[OUTPUT_VARS::BETA2] = 0.0;
    out_pw[OUTPUT_VARS::BETA3] = 0.0;

    out_pw[OUTPUT_VARS::G11] = quant_vals[CART_VARS::CART_G11];
    out_pw[OUTPUT_VARS::G12] = quant_vals[CART_VARS::CART_G12];
    out_pw[OUTPUT_VARS::G13] = quant_vals[CART_VARS::CART_G13];
    out_pw[OUTPUT_VARS::G22] = quant_vals[CART_VARS::CART_G22];
    out_pw[OUTPUT_VARS::G23] = quant_vals[CART_VARS::CART_G23];
    out_pw[OUTPUT_VARS::G33] = quant_vals[CART_VARS::CART_G33];
    out_pw[OUTPUT_VARS::K11] = quant_vals[CART_VARS::CART_K11];
    out_pw[OUTPUT_VARS::K12] = quant_vals[CART_VARS::CART_K12];
    out_pw[OUTPUT_VARS::K13] = quant_vals[CART_VARS::CART_K13];
    out_pw[OUTPUT_VARS::K22] = quant_vals[CART_VARS::CART_K22];
    out_pw[OUTPUT_VARS::K23] = quant_vals[CART_VARS::CART_K23];
    out_pw[OUTPUT_VARS::K33] = quant_vals[CART_VARS::CART_K33];

    double const H = quant_vals[CART_VARS::CART_H];
    double h = std::exp(H);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-14) {
      rho = 0.;
      eps = 0.;
      press = 0.;
      out_pw[OUTPUT_VARS::VEL1] = 0.;
      out_pw[OUTPUT_VARS::VEL2] = 0.;
      out_pw[OUTPUT_VARS::VEL3] = 0.;
    }
    else {
      rho = EOS<eos_t, DENSITY>::get(h);
      eps = EOS<eos_t, EPSILON>::get(h);
      press = EOS<eos_t, PRESSURE>::get(h);
    }
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
    out_pw[OUTPUT_VARS::VEL1] = quant_vals[CART_VARS::CART_VEL1];
    out_pw[OUTPUT_VARS::VEL2] = quant_vals[CART_VARS::CART_VEL2];
    out_pw[OUTPUT_VARS::VEL3] = quant_vals[CART_VARS::CART_VEL3];
    output_base = OUTPUT_BASIS::CARTESIAN;
    return out_pw;
  }

  /**
   * @brief Interface to export an array of OUTPUT_VARS for a given point.  The logic
   * for determining the EOS type is here and adds a bit of overhead as compared to
   * a templated class where eos_t is known at compile time.  Here we error on
   * convenience rather than speed since the cost is very small.
   * 
   * @tparam eos_t C++ Polytrope/Table type
   * @param x 
   * @param y 
   * @param z 
   * @return output_ary_t Interpolated solution at x,y,z
   */
  output_ary_t export_pointwise(double const & x, double const & y, double const & z);

  /**
   * @brief Export an array of OUTPUT_VARS for an array of npoints
   * 
   * @param npoints Number of coordinate points
   * @param xx 
   * @param yy 
   * @param zz 
   * @return grid_ary_t 
   */
  grid_ary_t export_coordinate_array(int const npoints, double const * xx, double const * yy, double const * zz);

  /**
   * @brief Export an array of OUTPUT_VARS for a given point in a Spherical basis
   * 
   * @tparam eos_t C++ Polytrope/Table type
   * @param x 
   * @param y 
   * @param z 
   * @return output_ary_t Interpolated solution at x,y,z
   */
  template<class eos_t>
  output_ary_t export_pointwise__spherical_imp(double const x, double const y, double const z) {        

    // Reset to NAN
    auto quant_vals__sph = interpolate_pointwise__spherical(x, y, z);

    const double r2 = x * x + y * y + z * z;
    const double r  = std::sqrt(r2);

    auto get_theta = [&]() {
      double theta;
      if(r2 < 1e-12) {
        theta = acos(z / 1e-10); //theta (angle from Z to xy plane)
      } else {
        theta = acos(z / r); //theta (angle from Z to xy plane)
      }
      return theta;
    };
    const double theta = get_theta();
    const double sint  = std::sin(theta);

    const double domega_dr = quant_vals__sph[ISO_VARS::ISO_DOMEGA_DR];
    const double domega_dt = quant_vals__sph[ISO_VARS::ISO_DOMEGA_DTHETA];

    const double N = quant_vals__sph[ISO_VARS::ISO_ALPHA];
    const double A = quant_vals__sph[ISO_VARS::ISO_METRIC_A];
    const double B = quant_vals__sph[ISO_VARS::ISO_METRIC_B];

    out_pw[OUTPUT_VARS::ALPHA] = N;

    out_pw[OUTPUT_VARS::BETA1] = 0.0;
    out_pw[OUTPUT_VARS::BETA2] = 0.0;
    out_pw[OUTPUT_VARS::BETA3] = -quant_vals__sph[ISO_VARS::ISO_METRIC_OMEGA];

    out_pw[OUTPUT_VARS::G11] = A * A;
    out_pw[OUTPUT_VARS::G12] = 0.0;
    out_pw[OUTPUT_VARS::G13] = 0.0;
    out_pw[OUTPUT_VARS::G22] = A * A * r2;
    out_pw[OUTPUT_VARS::G23] = 0.0;
    out_pw[OUTPUT_VARS::G33] = B * B * r2 * sint * sint;
    out_pw[OUTPUT_VARS::K11] = 0.0;
    out_pw[OUTPUT_VARS::K12] = 0.0;
    out_pw[OUTPUT_VARS::K13] = -out_pw[OUTPUT_VARS::G33] / 2.0 / N * domega_dr;
    out_pw[OUTPUT_VARS::K22] = 0.0;
    out_pw[OUTPUT_VARS::K23] = -out_pw[OUTPUT_VARS::G33] / 2.0 / N * domega_dt;
    out_pw[OUTPUT_VARS::K33] = 0.0;

    double const H = quant_vals__sph[ISO_VARS::ISO_H];
    double h = std::exp(H);
    double rho, eps, press, vphiU;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-14) {
      rho = 0.;
      eps = 0.;
      press = 0.;
      vphiU = 0.;
    }
    else {
      rho = EOS<eos_t, DENSITY>::get(h);
      eps = EOS<eos_t, EPSILON>::get(h);
      press = EOS<eos_t, PRESSURE>::get(h);
      vphiU = quant_vals__sph[ISO_VARS::ISO_U];
    }
    out_pw[OUTPUT_VARS::RHO]  = rho;
    out_pw[OUTPUT_VARS::EPS]  = eps;
    out_pw[OUTPUT_VARS::PRESS]= press;
    out_pw[OUTPUT_VARS::VEL1] = 0.0;
    out_pw[OUTPUT_VARS::VEL2] = 0.0;
    out_pw[OUTPUT_VARS::VEL3] = vphiU;
    output_base = OUTPUT_BASIS::SPHERICAL;
    return out_pw;
  }
   /**
   * @brief Interface to export an array of OUTPUT_VARS for a given point in a spherical
   * basis.  The logic for determining the EOS type is here and adds a bit of 
   * overhead as compared to a templated class where eos_t is known at compile time.  
   * Here we error on the side of convenience rather than speed since the cost is 
   * very small.
   * 
   * @tparam eos_t C++ Polytrope/Table type
   * @param x 
   * @param y 
   * @param z 
   * @return output_ary_t Interpolated solution at x,y,z
   */
  output_ary_t export_pointwise__spherical(double const & x, double const & y, double const & z);
};
}
