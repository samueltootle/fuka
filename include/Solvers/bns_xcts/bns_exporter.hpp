#include "Solvers/exporter.hpp"
#include "bin_ns.hpp"
namespace Kadath::FUKA_Solvers {

struct CFMS_BNS_Exporter : public Exporter<Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BIN_INFO>, Space_bin_ns> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BIN_INFO>;
  using space_t = Space_bin_ns;
  
  // Input ID grid functions that we interpolate on
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
    XCTS_H,
    XCTS_UX,
    XCTS_UY,
    XCTS_UZ,
    NUM_XCTS_VARS 
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

  using interp_ary_t = std::array<double, NUM_XCTS_VARS>; 
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>; 
  using grid_ary_t = std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  // Types from Base class
  using Exporter<config_t, space_t>::base_space_t;
  using Exporter<config_t, space_t>::base_config_t;

  // Class members from base class to avoid using this->
  using Exporter<config_t, space_t>::space;
  using Exporter<config_t, space_t>::bconfig;
  using Exporter<config_t, space_t>::ndom;

  // CFMS_BNS imported fields from file
  ptr_data_member(Scalar, conformal_factor, shared);
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Vector, shift, shared);
  ptr_data_member(Scalar, logh, shared);
  ptr_data_member(Scalar, velpotential, shared);

  // Constructed objects
  ptr_data_member(Tensor, A, shared);
  ptr_data_member(Vector, fluidvel, shared);

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
  
  CFMS_BNS_Exporter() : Exporter<config_t, space_t>(),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr), 
      logh(nullptr), velpotential(nullptr), fluidvel(nullptr) {}
  
  CFMS_BNS_Exporter(std::string config_filename) : Exporter<config_t, space_t>(config_filename),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr), 
      logh(nullptr), velpotential(nullptr), fluidvel(nullptr) {

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
  CFMS_BNS_Exporter(CFMS_BNS_Exporter const & r);
  CFMS_BNS_Exporter(CFMS_BNS_Exporter&& b) noexcept = delete;
  CFMS_BNS_Exporter& operator=(const CFMS_BNS_Exporter& b);

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
    
    quant_vals = interpolate_pointwise(x, y, z);
     
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

    double const H = quant_vals[XCTS_VARS::XCTS_H];
    double h = std::exp(H);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-12) {
      rho = 0.;
      eps = 0.;
      press = 0.;
    }
    else {
      rho = EOS<eos_t, eos_var_t::DENSITY>::get(h);
      eps = EOS<eos_t, eos_var_t::EPSILON>::get(h);
      press = EOS<eos_t, eos_var_t::PRESSURE>::get(h);
    }
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
    out_pw[OUTPUT_VARS::VELX]  = quant_vals[XCTS_VARS::XCTS_UX];
    out_pw[OUTPUT_VARS::VELY]  = quant_vals[XCTS_VARS::XCTS_UY];
    out_pw[OUTPUT_VARS::VELZ]  = quant_vals[XCTS_VARS::XCTS_UZ];
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
   * @tparam eos_t C++ Polytrope/Table type
   * @param npoints Number of coordinate points
   * @param xx 
   * @param yy 
   * @param zz 
   * @return grid_ary_t 
   */
  grid_ary_t export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz);
};
}
