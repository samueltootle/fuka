#include "Solvers/exporter.hpp"
namespace Kadath::FUKA_Solvers {

struct CFMS_NS_ISO_Exporter : public Exporter<Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_ISO_NS_INFO>, Space_spheric_adapted> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_ISO_NS_INFO>;
  using space_t = Space_spheric_adapted;
  
  // Input ID grid functions that we interpolate on
  enum ISO_VARS : size_t {
    ISO_LAPSE,
    ISO_METRIC_A,
    ISO_METRIC_B,
    ISO_OMEGA,
    ISO_DOMEGA_DR,
    ISO_DOMEGA_DTHETA,
    ISO_H,
    ISO_CARTX,
    ISO_CARTY,
    NUM_ISO_VARS 
  };

  // Output physical grid functions
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
    RHO,
    EPS,
    PRESS,
    VELX,
    VELY,
    VELZ,
    NUM_OUTPUT_VARS
  };

  using interp_ary_t = std::array<double, NUM_ISO_VARS>; 
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>; 
  using grid_ary_t = std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  std::vector<ISO_VARS> xcts_spacetime_indicies{
    ISO_LAPSE,
    ISO_METRIC_A,
    ISO_METRIC_B,
    ISO_OMEGA,
  };
  std::vector<ISO_VARS> xcts_fluid_indicies{
    ISO_H,
    ISO_U
  };  

  // Types from Base class
  using Exporter<config_t, space_t>::base_space_t;
  using Exporter<config_t, space_t>::base_config_t;

  // Class members from base class to avoid using this->
  using Exporter<config_t, space_t>::space;
  using Exporter<config_t, space_t>::bconfig;
  using Exporter<config_t, space_t>::ndom;

  // CFMS_NS imported fields from file
  ptr_data_member(Scalar, lap_Aterm, shared);
  ptr_data_member(Scalar, lap_Bterm, shared);
  ptr_data_member(Scalar, Nu, shared);
  ptr_data_member(Scalar, lap_omega_term, shared);
  ptr_data_member(Scalar, logh, shared);


  // Constructed objects
  // ptr_data_member(Tensor, A, shared);
  
  // dphi = (-y, +x, 0)
  ptr_data_member(Scalar, X, shared);
  ptr_data_member(Scalar, Y, shared);

  ptr_data_member(Scalar, metric_A, shared);
  ptr_data_member(Scalar, metric_B, shared);
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Scalar, omega, shared);
  ptr_data_member(Scalar, domega_dr, shared);
  ptr_data_member(Scalar, domega_dt, shared);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  interp_ary_t quant_vals;
  output_ary_t out_pw;
  bool export_ready{false};
  int const ndim{3};
  
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

  public:
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  bool is_export_ready() const { return export_ready; }  
  const int & get_ndim() const { return ndim; }
  
  CFMS_NS_Exporter() : Exporter<config_t, space_t>(),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr), logh(nullptr), fluidvel(nullptr) {}
  
  CFMS_NS_Exporter(std::string config_filename) :
    Exporter<config_t, space_t>(config_filename),
        conformal_factor(nullptr), lapse(nullptr), shift(nullptr), logh(nullptr), fluidvel(nullptr) {

    load_solution_from_file();
    initialize_eos();
    extract_computed_grid_functions();
    populate_quants();
    
    // This is to avoid a "bug" where "something" in kadath is not
    // correctly initialized prior to copying to other threads resulting
    // in undefined behavior.  By running the interpolator once, this
    // bug seems to be avoided.
    this->export_pointwise(0.5, 0., 0.);
  }

  /**
   * @brief Construct a new object - a mutex is used to allow for thread safety in
   * copying.
   *  
   */
  CFMS_NS_Exporter(CFMS_NS_Exporter const & r);
  CFMS_NS_Exporter(CFMS_NS_Exporter&& b) noexcept = delete;
  CFMS_NS_Exporter& operator=(const CFMS_NS_Exporter& b);

  public:

  /**
   * @brief For a given coordinate, interpolate the ID solution
   * 
   * @param x 
   * @param y 
   * @param z 
   * @return interp_ary_t 
   */
  interp_ary_t interpolate_pointwise(double const & x, double const & y, double const & z);
  
  /**
   * @brief For a given coordinate, interpolate the ID solution for a subset of the ID variables
   * 
   * @param x 
   * @param y 
   * @param z 
   * @param slice indicies of array subset to interpolate
   * @return interp_ary_t 
   */
  interp_ary_t interpolate_pointwise_subset(double const & x, double const & y, double const & z,
    std::vector<ISO_VARS> slice);
  
  /**
   * @brief Export: Load only the spacetime variables
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
    quant_vals = interpolate_pointwise_subset(x, y, z, xcts_fluid_indicies);

    double const H = quant_vals[ISO_VARS::XCTS_H];
    double h = std::exp(H);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-12) {
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
    out_pw[OUTPUT_VARS::VELX]  = quant_vals[ISO_VARS::XCTS_UX];
    out_pw[OUTPUT_VARS::VELY]  = quant_vals[ISO_VARS::XCTS_UY];
    out_pw[OUTPUT_VARS::VELZ]  = quant_vals[ISO_VARS::XCTS_UZ];
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

    /*
    * Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis
    * Code generated from NrPyv2
    */
    auto ADM_Spherical_to_Cart =[&]() {
    
      const double xCart[3] = {x, y, z};
      // Perform the basis transform on ADM vectors/tensors from Spherical to Cartesian:

      // Set destination xx[3] based on desired xCart[3]
      double xx0, xx1, xx2;
      /*
      *  Original SymPy expressions:
      *  "[xx0 = sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2)]"
      *  "[xx1 = acos(xCart[2]/sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2))]"
      *  "[xx2 = atan2(xCart[1], xCart[0])]"
      */
      {
        const REAL tmp0 = sqrt(((xCart[0]) * (xCart[0])) + ((xCart[1]) * (xCart[1])) + ((xCart[2]) * (xCart[2])));
        xx0 = tmp0;
        xx1 = acos(xCart[2] / tmp0);
        xx2 = atan2(xCart[1], xCart[0]);
      }
      // Unpack initial_data for ADM vectors/tensors
      const double N = quant_vals[ISO_VARS::ISO_ALPHA];
      const double domega_dr = quant_vals[ISO_VARS::ISO_DOMEGA_DR];
      const double domega_dt = quant_vals[ISO_VARS::ISO_DOMEGA_DTHETA];

      const double betaSpherical0 = 0.0;                              // r
      const double betaSpherical1 = 0.0;                              // theta
      const double betaSpherical2 = -quant_vals[ISO_VARS::ISO_OMEGA]; // phi

      const double A = quant_vals[ISO_VARS::ISO_METRIC_A];
      const double B = quant_vals[ISO_VARS::ISO_METRIC_B];

      const double gammaSphericalDD01 = 0.0;
      const double gammaSphericalDD02 = 0.0;
      const double gammaSphericalDD12 = 0.0;
      const double gammaSphericalDD00 = A * A;
      const double gammaSphericalDD11 = A * A * xx0 * xx0;
      const double gammaSphericalDD22 = B * B * xx0 * xx0 * sin(xx2) * sin(xx2);

      const double KSphericalDD00 = 0.0;
      const double KSphericalDD01 = 0.0;
      const double KSphericalDD11 = 0.0;
      const double KSphericalDD22 = 0.0;
      const double KSphericalDD02 = -gammaSphericalDD22 / 2.0 / N * domega_dr;
      const double KSphericalDD12 = -gammaSphericalDD22 / 2.0 / N * domega_dt;
      const REAL tmp0 = cos(xx2);
      const REAL tmp1 = sin(xx1);
      const REAL tmp3 = cos(xx1);
      const REAL tmp6 = sin(xx2);
      const REAL tmp12 = ((xx0) * (xx0));
      const REAL tmp5 = betaSpherical1 * tmp3 * xx0;
      const REAL tmp7 = tmp6 * xx0;
      const REAL tmp9 = tmp0 * xx0;
      const REAL tmp10 = ((tmp0) * (tmp0));
      const REAL tmp11 = ((tmp6) * (tmp6));
      const REAL tmp13 = ((tmp1) * (tmp1) * (tmp1));
      const REAL tmp15 = ((tmp3) * (tmp3));
      const REAL tmp20 = ((tmp1) * (tmp1) * (tmp1) * (tmp1)) * ((xx0) * (xx0) * (xx0) * (xx0));
      const REAL tmp25 = ((tmp1) * (tmp1));
      const REAL tmp33 = tmp1 * tmp3;
      const REAL tmp17 = tmp1 * tmp12 * tmp15;
      const REAL tmp21 = gammaSphericalDD00 * tmp20;
      const REAL tmp23 = tmp13 * tmp3 * ((xx0) * (xx0) * (xx0));
      const REAL tmp26 = tmp12 * tmp25;
      const REAL tmp29 = -tmp15 * tmp7 - tmp25 * tmp7;
      const REAL tmp34 = gammaSphericalDD12 * tmp33;
      const REAL tmp41 = tmp15 * tmp9 + tmp25 * tmp9;
      const REAL tmp48 = tmp33 * tmp9;
      const REAL tmp55 = tmp1 * tmp12 * tmp3;
      const REAL tmp66 = tmp33 * tmp7;
      const REAL tmp71 = KSphericalDD00 * tmp20;
      const REAL tmp18 = (1.0 / ((tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17) *
                                (tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17)));
      const REAL tmp24 = 2 * gammaSphericalDD01 * tmp23;
      const REAL tmp36 = gammaSphericalDD02 * tmp26;
      const REAL tmp50 = -tmp10 * tmp25 * xx0 - tmp11 * tmp25 * xx0;
      const REAL tmp56 = tmp10 * tmp55 + tmp11 * tmp55;
      const REAL tmp73 = 2 * KSphericalDD01 * tmp23;
      const REAL tmp77 = KSphericalDD02 * tmp26;
      const REAL tmp19 = tmp10 * tmp18;
      const REAL tmp28 = gammaSphericalDD11 * tmp15 * tmp26;
      const REAL tmp31 = gammaSphericalDD22 * tmp18;
      const REAL tmp32 = tmp18 * tmp29;
      const REAL tmp39 = tmp0 * tmp18;
      const REAL tmp42 = tmp18 * tmp41;
      const REAL tmp57 = tmp18 * tmp56;
      const REAL tmp63 = tmp11 * tmp18;
      const REAL tmp68 = tmp18 * ((tmp50) * (tmp50));
      const REAL tmp69 = tmp18 * ((tmp56) * (tmp56));
      const REAL tmp74 = KSphericalDD11 * tmp15 * tmp26;
      const REAL tmp75 = KSphericalDD22 * tmp18;
      const REAL tmp40 = tmp39 * tmp6;
      const REAL tmp52 = gammaSphericalDD11 * tmp18 * tmp50;
      const REAL tmp54 = gammaSphericalDD01 * tmp26 * tmp50;
      const REAL tmp80 = KSphericalDD11 * tmp18 * tmp50;
      const REAL tmp81 = KSphericalDD01 * tmp26 * tmp50;
      const REAL tmp60 = gammaSphericalDD00 * tmp26 * tmp57;
      const REAL tmp83 = KSphericalDD00 * tmp26 * tmp57;
      out_pw[OUTPUT_VARS::BETAX] = betaSpherical0 * tmp0 * tmp1 - betaSpherical2 * tmp1 * tmp7 + tmp0 * tmp5;
      out_pw[OUTPUT_VARS::BETAY] = betaSpherical0 * tmp1 * tmp6 + betaSpherical2 * tmp1 * tmp9 + tmp5 * tmp6;
      out_pw[OUTPUT_VARS::BETAZ] = betaSpherical0 * tmp3 - betaSpherical1 * tmp1 * xx0;
      out_pw[OUTPUT_VARS::GXX] =
          2 * tmp0 * tmp32 * tmp36 + tmp19 * tmp21 + tmp19 * tmp24 + tmp19 * tmp28 + ((tmp29) * (tmp29)) * tmp31 + 2 * tmp32 * tmp34 * tmp9;
      out_pw[OUTPUT_VARS::GXY] = tmp0 * tmp36 * tmp42 + tmp21 * tmp40 + tmp24 * tmp40 + tmp28 * tmp40 + tmp29 * tmp31 * tmp41 + tmp32 * tmp34 * tmp7 +
                                tmp32 * tmp36 * tmp6 + tmp34 * tmp42 * tmp9;
      out_pw[OUTPUT_VARS::GXZ] = gammaSphericalDD01 * tmp48 * tmp57 + gammaSphericalDD02 * tmp32 * tmp56 + gammaSphericalDD12 * tmp32 * tmp50 +
                                tmp0 * tmp60 + tmp39 * tmp54 + tmp48 * tmp52;
      out_pw[OUTPUT_VARS::GYY] =
          tmp21 * tmp63 + tmp24 * tmp63 + tmp28 * tmp63 + tmp31 * ((tmp41) * (tmp41)) + 2 * tmp34 * tmp42 * tmp7 + 2 * tmp36 * tmp42 * tmp6;
      out_pw[OUTPUT_VARS::GYZ] = gammaSphericalDD01 * tmp57 * tmp66 + gammaSphericalDD02 * tmp42 * tmp56 + gammaSphericalDD12 * tmp42 * tmp50 +
                                tmp18 * tmp54 * tmp6 + tmp52 * tmp66 + tmp6 * tmp60;
      out_pw[OUTPUT_VARS::GZZ] = gammaSphericalDD00 * tmp69 + 2 * gammaSphericalDD01 * tmp50 * tmp57 + gammaSphericalDD11 * tmp68;
      out_pw[OUTPUT_VARS::KXX] =
          2 * KSphericalDD12 * tmp32 * tmp48 + 2 * tmp0 * tmp32 * tmp77 + tmp19 * tmp71 + tmp19 * tmp73 + tmp19 * tmp74 + ((tmp29) * (tmp29)) * tmp75;
      out_pw[OUTPUT_VARS::KXY] = KSphericalDD12 * tmp32 * tmp33 * tmp7 + KSphericalDD12 * tmp33 * tmp42 * tmp9 + tmp0 * tmp42 * tmp77 +
                                tmp29 * tmp41 * tmp75 + tmp32 * tmp6 * tmp77 + tmp40 * tmp71 + tmp40 * tmp73 + tmp40 * tmp74;
      out_pw[OUTPUT_VARS::KXZ] =
          KSphericalDD01 * tmp48 * tmp57 + KSphericalDD02 * tmp32 * tmp56 + KSphericalDD12 * tmp32 * tmp50 + tmp0 * tmp83 + tmp39 * tmp81 + tmp48 * tmp80;
      out_pw[OUTPUT_VARS::KYY] =
          2 * KSphericalDD12 * tmp42 * tmp66 + ((tmp41) * (tmp41)) * tmp75 + 2 * tmp42 * tmp6 * tmp77 + tmp63 * tmp71 + tmp63 * tmp73 + tmp63 * tmp74;
      out_pw[OUTPUT_VARS::KYZ] = KSphericalDD01 * tmp57 * tmp66 + KSphericalDD02 * tmp42 * tmp56 + KSphericalDD12 * tmp42 * tmp50 + tmp18 * tmp6 * tmp81 +
                                tmp6 * tmp83 + tmp66 * tmp80;
      out_pw[OUTPUT_VARS::KZZ] = KSphericalDD00 * tmp69 + 2 * KSphericalDD01 * tmp50 * tmp57 + KSphericalDD11 * tmp68;
    };

    double const N = quant_vals[ISO_VARS::ISO_ALPHA];
    double const omega = quant_vals[ISO_VARS::ISO_OMEGA];
    const double Omega = (*bconfig)(BCO_PARAMS::OMEGA);
    const double U_factor = (Omega - omega) / N;

    double const H = quant_vals[ISO_VARS::ISO_H];
    double h = std::exp(H);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-12) {
      rho = 0.;
      eps = 0.;
      press = 0.;
    }
    else {
      rho = EOS<eos_t, DENSITY>::get(h);
      eps = EOS<eos_t, EPSILON>::get(h);
      press = EOS<eos_t, PRESSURE>::get(h);
    }
    out_pw[OUTPUT_VARS::ALPHA] = N;
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
    out_pw[OUTPUT_VARS::VELX]  = - quant_vals[ISO_VARS::ISO_CARTY] * U_factor;
    out_pw[OUTPUT_VARS::VELY]  =   quant_vals[ISO_VARS::ISO_CARTX] * U_factor;
    out_pw[OUTPUT_VARS::VELZ]  = 0.;
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
};
}
