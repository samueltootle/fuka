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
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>; 
  using grid_ary_t = std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  std::vector<ISO_VARS> iso_spacetime_indicies{
    ISO_ALPHA,
    ISO_METRIC_A,
    ISO_METRIC_B,
    ISO_METRIC_OMEGA,
    ISO_DOMEGA_DR,
    ISO_DOMEGA_DTHETA,
  };
  std::vector<ISO_VARS> iso_fluid_indicies{
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
  ptr_data_member(Scalar, omega, shared);


  // Constructed objects
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Scalar, metric_A, shared);
  ptr_data_member(Scalar, metric_B, shared);  
  ptr_data_member(Scalar, metric_omega, shared);
  ptr_data_member(Scalar, domega_dr, shared);
  ptr_data_member(Scalar, domega_dt, shared);
  ptr_data_member(Scalar, fluidvel, shared);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  interp_ary_t quant_vals;
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

  public:
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  bool is_export_ready() const { return export_ready; }  
  const int & get_ndim() const { return ndim; }
  
  CFMS_NS_ISO_Exporter() : Exporter<config_t, space_t>(),
    lap_Aterm(nullptr), lap_Bterm(nullptr), Nu(nullptr), lap_omega_term(nullptr), logh(nullptr), 
      metric_A(nullptr), metric_B(nullptr), lapse(nullptr), metric_omega(nullptr), 
        domega_dr(nullptr), domega_dt(nullptr), fluidvel(nullptr) {}
  
  CFMS_NS_ISO_Exporter(std::string config_filename) :
    Exporter<config_t, space_t>(config_filename),
    lap_Aterm(nullptr), lap_Bterm(nullptr), Nu(nullptr), lap_omega_term(nullptr), logh(nullptr), 
      metric_A(nullptr), metric_B(nullptr), lapse(nullptr), metric_omega(nullptr), 
        domega_dr(nullptr), domega_dt(nullptr), fluidvel(nullptr) {

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
  
  // /**
  //  * @brief Export: Load only the spacetime variables
  //  * 
  //  * @param x 
  //  * @param y 
  //  * @param z
  //  * @return output_ary_t 
  //  */
  // output_ary_t export_pointwise_spacetime_vars(double const & x, double const & y, double const & z);

  // /**
  //  * @brief Export: Load only the fluid variables
  //  * 
  //  * @tparam eos_t EOS type
  //  * @param x 
  //  * @param y 
  //  * @param z
  //  * @return output_ary_t 
  //  */
  // template<class eos_t>
  // output_ary_t export_pointwise_fluid_vars_imp(double const & x, double const & y, double const & z) {
    
  //   // Reset to NAN
  //   quant_vals.fill(NAN);
  //   quant_vals = interpolate_pointwise_subset(x, y, z, iso_fluid_indicies);

  //   double const H = quant_vals[ISO_VARS::ISO_H];
  //   double h = std::exp(H);
  //   double rho, eps, press;

  //   // get quantities point-wise, since h is smoothest, and cut data at H=0
  //   if(std::fabs(H) <= 1e-12) {
  //     rho = 0.;
  //     eps = 0.;
  //     press = 0.;
  //   }
  //   else {
  //     rho = EOS<eos_t, DENSITY>::get(h);
  //     eps = EOS<eos_t, EPSILON>::get(h);
  //     press = EOS<eos_t, PRESSURE>::get(h);
  //   }
  //   out_pw[OUTPUT_VARS::RHO]   = rho;
  //   out_pw[OUTPUT_VARS::EPS]   = eps;
  //   out_pw[OUTPUT_VARS::PRESS] = press;
  //   // FIXME
  //   // out_pw[OUTPUT_VARS::VEL1]  = quant_vals[ISO_VARS::XCTS_UX];
  //   // out_pw[OUTPUT_VARS::VEL2]  = quant_vals[ISO_VARS::XCTS_UY];
  //   // out_pw[OUTPUT_VARS::VEL3]  = quant_vals[ISO_VARS::XCTS_UZ];
  //   return out_pw;
  // }

  // /**
  //  * @brief Export: Interface for only loading the fluid variables
  //  * 
  //  * @tparam eos_t EOS type
  //  * @param x 
  //  * @param y 
  //  * @param z
  //  * @return output_ary_t 
  //  */
  // output_ary_t export_pointwise_fluid_vars(double const & x, double const & y, double const & z);

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
    auto ADM_Spherical_to_Cart = [&] () {
      using REAL = double;
      const REAL xCart[3] = {x, y, z};
      // Perform the basis transform on ADM vectors/tensors from Spherical to Cartesian:

      // Set destination xx[3] based on desired xCart[3]
      REAL xx0, xx1, xx2;
      /*
      *  Original SymPy expressions:
      *  "[xx0 = sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2)]"
      *  "[xx1 = acos(xCart[2]/sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2))]"
      *  "[xx2 = atan2(xCart[1], xCart[0])]"
      */
      {
        const REAL tmp0 = sqrt(((xCart[0]) * (xCart[0])) + ((xCart[1]) * (xCart[1])) + ((xCart[2]) * (xCart[2])));
        xx0 = tmp0;
        if(std::fabs(xx0) < 1e-12) {
          xx0 = 1e-12;
          xx1 = 1e-12;
          xx2 = 1e-12;
        } else if(std::fabs(x) < 1e-12 && std::fabs(y) < 1e-12) {
          xx1 = 1e-12;
          xx2 = 1e-12;
        } else {
          xx1 = acos(xCart[2] / xx0); //theta (angle from Z to xy plane)
          xx2 = atan2(xCart[1], xCart[0]); // phi (angle from x to y axis)
        }
        // cout << "r: " << xx0 << ", t: " << xx1 << ", phi: " << xx2 << endl;
      }
      // Unpack initial_data for ADM vectors/tensors
      const REAL N = quant_vals[ISO_VARS::ISO_ALPHA];
      const REAL U_factor = quant_vals[ISO_VARS::ISO_U];
      const REAL domega_dr = quant_vals[ISO_VARS::ISO_DOMEGA_DR];
      const REAL domega_dt = quant_vals[ISO_VARS::ISO_DOMEGA_DTHETA];

      const REAL FluidVelU0 = 0.0;                                   // r
      const REAL FluidVelU1 = 0.0;                                   // theta
      const REAL FluidVelU2 = U_factor;                              // phi
      const REAL betaSphericalU0 = 0.0;                              // r
      const REAL betaSphericalU1 = 0.0;                              // theta
      const REAL betaSphericalU2 = -quant_vals[ISO_VARS::ISO_METRIC_OMEGA]; // phi

      const REAL A = quant_vals[ISO_VARS::ISO_METRIC_A];
      const REAL B = quant_vals[ISO_VARS::ISO_METRIC_B];

      const REAL gammaSphericalDD01 = 0.0;
      const REAL gammaSphericalDD02 = 0.0;
      const REAL gammaSphericalDD12 = 0.0;
      const REAL gammaSphericalDD00 = A * A;
      const REAL gammaSphericalDD11 = A * A * xx0 * xx0;
      const REAL gammaSphericalDD22 = B * B * xx0 * xx0 * sin(xx1) * sin(xx1);
      // cout << "g11: " << gammaSphericalDD00 << ", g22: " << gammaSphericalDD11 << ", g33: " << gammaSphericalDD22 << endl;

      const REAL KSphericalDD00 = 0.0;
      const REAL KSphericalDD01 = 0.0;
      const REAL KSphericalDD11 = 0.0;
      const REAL KSphericalDD22 = 0.0;
      const REAL KSphericalDD02 = -gammaSphericalDD22 / 2.0 / N * domega_dr;
      const REAL KSphericalDD12 = -gammaSphericalDD22 / 2.0 / N * domega_dt;
      const REAL tmp0 = cos(xx2);
      const REAL tmp1 = sin(xx1);
      const REAL tmp4 = cos(xx1);
      const REAL tmp6 = sin(xx2);
      const REAL tmp12 = ((xx0) * (xx0));
      const REAL tmp3 = tmp0 * xx0;
      const REAL tmp7 = tmp1 * xx0;
      const REAL tmp9 = tmp6 * xx0;
      const REAL tmp10 = ((tmp0) * (tmp0));
      const REAL tmp11 = ((tmp6) * (tmp6));
      const REAL tmp13 = ((tmp1) * (tmp1) * (tmp1));
      const REAL tmp15 = ((tmp4) * (tmp4));
      const REAL tmp20 = ((tmp1) * (tmp1) * (tmp1) * (tmp1)) * ((xx0) * (xx0) * (xx0) * (xx0));
      const REAL tmp25 = ((tmp1) * (tmp1));
      const REAL tmp17 = tmp1 * tmp12 * tmp15;
      const REAL tmp21 = gammaSphericalDD00 * tmp20;
      const REAL tmp23 = tmp13 * tmp4 * ((xx0) * (xx0) * (xx0));
      const REAL tmp26 = tmp12 * tmp25;
      const REAL tmp29 = -tmp15 * tmp9 - tmp25 * tmp9;
      const REAL tmp34 = tmp4 * tmp7;
      const REAL tmp41 = tmp15 * tmp3 + tmp25 * tmp3;
      const REAL tmp52 = tmp1 * tmp12 * tmp4;
      const REAL tmp68 = KSphericalDD00 * tmp20;
      const REAL tmp18 = (1.0 / ((tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17) *
                                (tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17)));
      const REAL tmp24 = 2 * gammaSphericalDD01 * tmp23;
      const REAL tmp35 = gammaSphericalDD12 * tmp34;
      const REAL tmp37 = gammaSphericalDD02 * tmp26;
      const REAL tmp47 = -tmp10 * tmp25 * xx0 - tmp11 * tmp25 * xx0;
      const REAL tmp53 = tmp10 * tmp52 + tmp11 * tmp52;
      const REAL tmp70 = 2 * KSphericalDD01 * tmp23;
      const REAL tmp73 = KSphericalDD12 * tmp34;
      const REAL tmp75 = KSphericalDD02 * tmp26;
      const REAL tmp19 = tmp10 * tmp18;
      const REAL tmp28 = gammaSphericalDD11 * tmp15 * tmp26;
      const REAL tmp31 = gammaSphericalDD22 * tmp18;
      const REAL tmp32 = tmp18 * tmp29;
      const REAL tmp39 = tmp0 * tmp18;
      const REAL tmp42 = tmp18 * tmp41;
      const REAL tmp54 = tmp18 * tmp53;
      const REAL tmp61 = tmp11 * tmp18;
      const REAL tmp64 = tmp18 * tmp6;
      const REAL tmp65 = tmp18 * ((tmp47) * (tmp47));
      const REAL tmp66 = tmp18 * ((tmp53) * (tmp53));
      const REAL tmp71 = KSphericalDD11 * tmp15 * tmp26;
      const REAL tmp72 = KSphericalDD22 * tmp18;
      const REAL tmp33 = tmp0 * tmp32;
      const REAL tmp40 = tmp39 * tmp6;
      const REAL tmp43 = tmp0 * tmp42;
      const REAL tmp44 = tmp32 * tmp6;
      const REAL tmp49 = gammaSphericalDD11 * tmp34 * tmp47;
      const REAL tmp51 = gammaSphericalDD01 * tmp26 * tmp47;
      const REAL tmp63 = tmp42 * tmp6;
      const REAL tmp77 = KSphericalDD11 * tmp34 * tmp47;
      const REAL tmp78 = KSphericalDD01 * tmp26 * tmp47;
      const REAL tmp56 = gammaSphericalDD01 * tmp34 * tmp54;
      const REAL tmp58 = gammaSphericalDD00 * tmp26 * tmp54;
      const REAL tmp79 = KSphericalDD01 * tmp34 * tmp54;
      const REAL tmp80 = KSphericalDD00 * tmp26 * tmp54;
      out_pw[OUTPUT_VARS::BETA1] = betaSphericalU0 * tmp0 * tmp1 + betaSphericalU1 * tmp3 * tmp4 - betaSphericalU2 * tmp6 * tmp7;
      out_pw[OUTPUT_VARS::BETA2] = betaSphericalU0 * tmp1 * tmp6 + betaSphericalU1 * tmp4 * tmp9 + betaSphericalU2 * tmp0 * tmp7;
      out_pw[OUTPUT_VARS::BETA3] = betaSphericalU0 * tmp4 - betaSphericalU1 * tmp7;
      out_pw[OUTPUT_VARS::G11] = tmp19 * tmp21 + tmp19 * tmp24 + tmp19 * tmp28 + ((tmp29) * (tmp29)) * tmp31 + 2 * tmp33 * tmp35 + 2 * tmp33 * tmp37;
      out_pw[OUTPUT_VARS::G12] =
          tmp21 * tmp40 + tmp24 * tmp40 + tmp28 * tmp40 + tmp29 * tmp31 * tmp41 + tmp35 * tmp43 + tmp35 * tmp44 + tmp37 * tmp43 + tmp37 * tmp44;
      out_pw[OUTPUT_VARS::G13] =
          gammaSphericalDD02 * tmp32 * tmp53 + gammaSphericalDD12 * tmp32 * tmp47 + tmp0 * tmp56 + tmp0 * tmp58 + tmp39 * tmp49 + tmp39 * tmp51;
      out_pw[OUTPUT_VARS::G22] = tmp21 * tmp61 + tmp24 * tmp61 + tmp28 * tmp61 + tmp31 * ((tmp41) * (tmp41)) + 2 * tmp35 * tmp63 + 2 * tmp37 * tmp63;
      out_pw[OUTPUT_VARS::G23] =
          gammaSphericalDD02 * tmp42 * tmp53 + gammaSphericalDD12 * tmp42 * tmp47 + tmp49 * tmp64 + tmp51 * tmp64 + tmp56 * tmp6 + tmp58 * tmp6;
      out_pw[OUTPUT_VARS::G33] = gammaSphericalDD00 * tmp66 + 2 * gammaSphericalDD01 * tmp47 * tmp54 + gammaSphericalDD11 * tmp65;
      out_pw[OUTPUT_VARS::K11] = tmp19 * tmp68 + tmp19 * tmp70 + tmp19 * tmp71 + ((tmp29) * (tmp29)) * tmp72 + 2 * tmp33 * tmp73 + 2 * tmp33 * tmp75;
      out_pw[OUTPUT_VARS::K12] =
          tmp29 * tmp41 * tmp72 + tmp40 * tmp68 + tmp40 * tmp70 + tmp40 * tmp71 + tmp43 * tmp73 + tmp43 * tmp75 + tmp44 * tmp73 + tmp44 * tmp75;
      out_pw[OUTPUT_VARS::K13] =
          KSphericalDD02 * tmp32 * tmp53 + KSphericalDD12 * tmp32 * tmp47 + tmp0 * tmp79 + tmp0 * tmp80 + tmp39 * tmp77 + tmp39 * tmp78;
      out_pw[OUTPUT_VARS::K22] = ((tmp41) * (tmp41)) * tmp72 + tmp61 * tmp68 + tmp61 * tmp70 + tmp61 * tmp71 + 2 * tmp63 * tmp73 + 2 * tmp63 * tmp75;
      out_pw[OUTPUT_VARS::K23] =
          KSphericalDD02 * tmp42 * tmp53 + KSphericalDD12 * tmp42 * tmp47 + tmp6 * tmp79 + tmp6 * tmp80 + tmp64 * tmp77 + tmp64 * tmp78;
      out_pw[OUTPUT_VARS::K33] = KSphericalDD00 * tmp66 + 2 * KSphericalDD01 * tmp47 * tmp54 + KSphericalDD11 * tmp65;
      out_pw[OUTPUT_VARS::VEL1] = FluidVelU0 * tmp0 * tmp1 + FluidVelU1 * tmp3 * tmp4 - FluidVelU2 * tmp6 * tmp7;
      out_pw[OUTPUT_VARS::VEL2] = FluidVelU0 * tmp1 * tmp6 + FluidVelU1 * tmp4 * tmp9 + FluidVelU2 * tmp0 * tmp7;
      out_pw[OUTPUT_VARS::VEL3] = FluidVelU0 * tmp4 - FluidVelU1 * tmp7;
    };
    ADM_Spherical_to_Cart();

    double const H = quant_vals[ISO_VARS::ISO_H];
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
    out_pw[OUTPUT_VARS::ALPHA] = quant_vals[ISO_VARS::ISO_ALPHA];
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
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
  output_ary_t export_pointwise__spherical_imp(double const & x, double const & y, double const & z) {        

    // Reset to NAN
    quant_vals.fill(NAN);
    quant_vals = interpolate_pointwise(x, y, z);

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

    const double domega_dr = quant_vals[ISO_VARS::ISO_DOMEGA_DR];
    const double domega_dt = quant_vals[ISO_VARS::ISO_DOMEGA_DTHETA];

    const double N = quant_vals[ISO_VARS::ISO_ALPHA];
    const double A = quant_vals[ISO_VARS::ISO_METRIC_A];
    const double B = quant_vals[ISO_VARS::ISO_METRIC_B];

    out_pw[OUTPUT_VARS::ALPHA] = N;

    out_pw[OUTPUT_VARS::BETA1] = 0.0;
    out_pw[OUTPUT_VARS::BETA2] = 0.0;
    out_pw[OUTPUT_VARS::BETA3] = -quant_vals[ISO_VARS::ISO_METRIC_OMEGA];

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

    double const H = quant_vals[ISO_VARS::ISO_H];
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
      vphiU = quant_vals[ISO_VARS::ISO_U];
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
