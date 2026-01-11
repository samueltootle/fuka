#include "EOS/FUKA_EOS_Utilities.hh"
#include "Solvers/exporter.hpp"
#include "exporter_utilities.hpp"

namespace Kadath::FUKA_Solvers {

struct CFMS_NS_ISO_Exporter
    : public Exporter<Kadath::FUKA_Config::kadath_config_boost<
                          Kadath::FUKA_Config::BCO_ISO_NS_INFO>,
                      Space_polar_adapted> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<
      Kadath::FUKA_Config::BCO_ISO_NS_INFO>;
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

  // clang-format off
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
  // clang-format on

  using interp_ary_t = std::array<double, NUM_ISO_VARS>;
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>;
  using grid_ary_t =
      std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

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
   * @brief We initialize the FUKA EOS module here to avoid needing to handle
   * this by the an interface code on the importer side since the FUKA EOS
   * module is independent of the EOS module used in the evolution code.
   */
  void initialize_eos();

  /**
   * @brief Here we load the ID solution from file in order to initialize the ID
   * fields and Exporter variables
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
  template <class eos_t>
  void set_eos_ope(System_of_eqs& syst, Param& p) {
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
  }

  /**
   * @brief Compute needed grid functions (gf) within Kadath's System_of_eqs
   * before populating the gf pointers
   *
   */
  void extract_computed_grid_functions();

  /**
   * @brief Using the populated struct pointers, populate the "quants" array
   * with the associated gf references
   *
   */
  void populate_quants();

 public:
  std::vector<std::reference_wrapper<const Scalar>> const& get_quants() const {
    return quants;
  }

  bool is_export_ready() const { return export_ready; }

  const int& get_ndim() const { return ndim; }

  CFMS_NS_ISO_Exporter()
      : Exporter<config_t, space_t>(),
        lap_Aterm(nullptr),
        lap_Bterm(nullptr),
        Nu(nullptr),
        lap_omega_term(nullptr),
        logh(nullptr),
        metric_A(nullptr),
        metric_B(nullptr),
        lapse(nullptr),
        metric_omega(nullptr),
        domega_dr(nullptr),
        domega_dt(nullptr),
        fluidvel(nullptr) {}

  CFMS_NS_ISO_Exporter(std::string config_filename)
      : Exporter<config_t, space_t>(config_filename),
        lap_Aterm(nullptr),
        lap_Bterm(nullptr),
        Nu(nullptr),
        lap_omega_term(nullptr),
        logh(nullptr),
        metric_A(nullptr),
        metric_B(nullptr),
        lapse(nullptr),
        metric_omega(nullptr),
        domega_dr(nullptr),
        domega_dt(nullptr),
        fluidvel(nullptr) {
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
   * @brief Construct a new object - a mutex is used to allow for thread safety
   * in copying.
   *
   */
  CFMS_NS_ISO_Exporter(CFMS_NS_ISO_Exporter const& r);
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
  interp_ary_t interpolate_pointwise(double const& x,
                                     double const& y,
                                     double const& z);

  /**
   * @brief For a given coordinate, interpolate the ID solution for a subset of
   * the ID variables
   *
   * @param x
   * @param y
   * @param z
   * @param slice indicies of array subset to interpolate
   * @return interp_ary_t
   */
  interp_ary_t interpolate_pointwise_subset(double const& x,
                                            double const& y,
                                            double const& z,
                                            std::vector<ISO_VARS> slice);

  // /**
  //  * @brief Export: Load only the spacetime variables
  //  *
  //  * @param x
  //  * @param y
  //  * @param z
  //  * @return output_ary_t
  //  */
  // output_ary_t export_pointwise_spacetime_vars(double const & x, double const
  // & y, double const & z);

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
  // output_ary_t export_pointwise_fluid_vars_imp(double const & x, double const
  // & y, double const & z) {

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
  // output_ary_t export_pointwise_fluid_vars(double const & x, double const &
  // y, double const & z);

  /**
   * @brief Export an array of OUTPUT_VARS for a given point
   *
   * @tparam eos_t C++ Polytrope/Table type
   * @param x
   * @param y
   * @param z
   * @return output_ary_t Interpolated solution at x,y,z
   */
  template <class eos_t>
  struct export_pointwise_imp {
    friend CFMS_NS_ISO_Exporter;

    output_ary_t operator()(CFMS_NS_ISO_Exporter& base,
                            double const& x,
                            double const& y,
                            double const& z) {
      // Reset to NAN
      base.quant_vals.fill(NAN);
      base.quant_vals = base.interpolate_pointwise(x, y, z);

      export_utils::basis_transform_spherical_tofrom_cart basis_transform;
      basis_transform.set__cart(x, y, z);
      double const theta = basis_transform.get_theta();
      double const rsq = basis_transform.get_r_sq();

      double const N = base.quant_vals[ISO_VARS::ISO_ALPHA];
      double const U_factor = base.quant_vals[ISO_VARS::ISO_U];
      double const domega_dr = base.quant_vals[ISO_VARS::ISO_DOMEGA_DR];
      double const domega_dt = base.quant_vals[ISO_VARS::ISO_DOMEGA_DTHETA];

      double const A = base.quant_vals[ISO_VARS::ISO_METRIC_A];
      double const B = base.quant_vals[ISO_VARS::ISO_METRIC_B];

      double const Asq = A * A;
      double const Bsq = B * B;

      // clang-format off
      export_utils::basis_transform_spherical_tofrom_cart::vector_t FluidVelU_Sph = {
        0.0,
        0.0,
        U_factor
      };

      export_utils::basis_transform_spherical_tofrom_cart::vector_t betaU_Sph = {
        0.0,
        0.0,
        -base.quant_vals[ISO_VARS::ISO_METRIC_OMEGA]
      };

      export_utils::basis_transform_spherical_tofrom_cart::matrix_t gammaDD_Sph = {{
        {{ Asq, 0.0      , 0.0 }},
        {{ 0.0, Asq * rsq, 0.0 }},
        {{ 0.0, 0.0      , Bsq * rsq * sin(theta) * sin(theta) }}
      }};
      
      double const KDD_rphi  = -gammaDD_Sph[2][2] / 2.0 / N * domega_dr;
      double const KDD_thphi = -gammaDD_Sph[2][2] / 2.0 / N * domega_dt;
      export_utils::basis_transform_spherical_tofrom_cart::matrix_t KDD_Sph = {{
        {{ 0.0     , 0.0      , KDD_rphi }},
        {{ 0.0     , 0.0      , KDD_thphi}},
        {{ KDD_rphi, KDD_thphi, 0.0      }}
      }};
      // clang-format on

      auto const gammaDD_Cart =
          basis_transform.matrixDD__sph_to_cart(gammaDD_Sph);
      auto const KDD_Cart = basis_transform.matrixDD__sph_to_cart(KDD_Sph);
      auto const betaU_Cart = basis_transform.vectorU__sph_to_cart(betaU_Sph);
      auto const FluidVelU_Cart =
          basis_transform.vectorU__sph_to_cart(FluidVelU_Sph);
      base.out_pw[OUTPUT_VARS::BETA1] = betaU_Cart[0];
      base.out_pw[OUTPUT_VARS::BETA2] = betaU_Cart[1];
      base.out_pw[OUTPUT_VARS::BETA3] = betaU_Cart[2];
      base.out_pw[OUTPUT_VARS::G11] = gammaDD_Cart[0][0];
      base.out_pw[OUTPUT_VARS::G12] = gammaDD_Cart[0][1];
      base.out_pw[OUTPUT_VARS::G13] = gammaDD_Cart[0][2];
      base.out_pw[OUTPUT_VARS::G22] = gammaDD_Cart[1][1];
      base.out_pw[OUTPUT_VARS::G23] = gammaDD_Cart[1][2];
      base.out_pw[OUTPUT_VARS::G33] = gammaDD_Cart[2][2];
      base.out_pw[OUTPUT_VARS::K11] = KDD_Cart[0][0];
      base.out_pw[OUTPUT_VARS::K12] = KDD_Cart[0][1];
      base.out_pw[OUTPUT_VARS::K13] = KDD_Cart[0][2];
      base.out_pw[OUTPUT_VARS::K22] = KDD_Cart[1][1];
      base.out_pw[OUTPUT_VARS::K23] = KDD_Cart[1][2];
      base.out_pw[OUTPUT_VARS::K33] = KDD_Cart[2][2];
      base.out_pw[OUTPUT_VARS::VEL1] = FluidVelU_Cart[0];
      base.out_pw[OUTPUT_VARS::VEL2] = FluidVelU_Cart[1];
      base.out_pw[OUTPUT_VARS::VEL3] = FluidVelU_Cart[2];

      double const H = base.quant_vals[ISO_VARS::ISO_H];
      double h = std::exp(H);
      double rho, eps, press;

      // get quantities point-wise, since h is smoothest, and cut data at H=0
      if (std::fabs(H) <= 1e-14) {
        rho = 0.;
        eps = 0.;
        press = 0.;
        base.out_pw[OUTPUT_VARS::VEL1] = 0.;
        base.out_pw[OUTPUT_VARS::VEL2] = 0.;
        base.out_pw[OUTPUT_VARS::VEL3] = 0.;
      } else {
        rho = EOS<eos_t, DENSITY>::get(h);
        eps = EOS<eos_t, EPSILON>::get(h);
        press = EOS<eos_t, PRESSURE>::get(h);
      }
      base.out_pw[OUTPUT_VARS::ALPHA] = base.quant_vals[ISO_VARS::ISO_ALPHA];
      base.out_pw[OUTPUT_VARS::RHO] = rho;
      base.out_pw[OUTPUT_VARS::EPS] = eps;
      base.out_pw[OUTPUT_VARS::PRESS] = press;
      base.output_base = OUTPUT_BASIS::CARTESIAN;
      return base.out_pw;
    }
  };

  /**
   * @brief Interface to export an array of OUTPUT_VARS for a given point.  The
   * logic for determining the EOS type is here and adds a bit of overhead as
   * compared to a templated class where eos_t is known at compile time.  Here
   * we error on convenience rather than speed since the cost is very small.
   *
   * @tparam eos_t C++ Polytrope/Table type
   * @param x
   * @param y
   * @param z
   * @return output_ary_t Interpolated solution at x,y,z
   */
  output_ary_t export_pointwise(double const& x,
                                double const& y,
                                double const& z);

  /**
   * @brief Export an array of OUTPUT_VARS for an array of npoints
   *
   * @param npoints Number of coordinate points
   * @param xx
   * @param yy
   * @param zz
   * @return grid_ary_t
   */
  grid_ary_t export_coordinate_array(int const npoints,
                                     double const* xx,
                                     double const* yy,
                                     double const* zz);

  /**
   * @brief Export an array of OUTPUT_VARS for a given point in a Spherical
   * basis
   *
   * @tparam eos_t C++ Polytrope/Table type
   * @param x
   * @param y
   * @param z
   * @return output_ary_t Interpolated solution at x,y,z
   */
  template <class eos_t>
  struct export_pointwise__spherical_imp {
    friend CFMS_NS_ISO_Exporter;

    output_ary_t operator()(CFMS_NS_ISO_Exporter& base,
                            double const& x,
                            double const& y,
                            double const& z) {
      // Reset to NAN
      base.quant_vals.fill(NAN);
      base.quant_vals = base.interpolate_pointwise(x, y, z);

      const double r2 = x * x + y * y + z * z;
      const double r = std::sqrt(r2);

      auto get_theta = [&]() {
        double theta;
        if (r2 < 1e-12) {
          theta = acos(z / 1e-10);  // theta (angle from Z to xy plane)
        } else {
          theta = acos(z / r);  // theta (angle from Z to xy plane)
        }
        return theta;
      };
      const double theta = get_theta();
      const double sint = std::sin(theta);

      const double domega_dr = base.quant_vals[ISO_VARS::ISO_DOMEGA_DR];
      const double domega_dt = base.quant_vals[ISO_VARS::ISO_DOMEGA_DTHETA];

      const double N = base.quant_vals[ISO_VARS::ISO_ALPHA];
      const double A = base.quant_vals[ISO_VARS::ISO_METRIC_A];
      const double B = base.quant_vals[ISO_VARS::ISO_METRIC_B];

      base.out_pw[OUTPUT_VARS::ALPHA] = N;

      base.out_pw[OUTPUT_VARS::BETA1] = 0.0;
      base.out_pw[OUTPUT_VARS::BETA2] = 0.0;
      base.out_pw[OUTPUT_VARS::BETA3] =
          -base.quant_vals[ISO_VARS::ISO_METRIC_OMEGA];

      base.out_pw[OUTPUT_VARS::G11] = A * A;
      base.out_pw[OUTPUT_VARS::G12] = 0.0;
      base.out_pw[OUTPUT_VARS::G13] = 0.0;
      base.out_pw[OUTPUT_VARS::G22] = A * A * r2;
      base.out_pw[OUTPUT_VARS::G23] = 0.0;
      base.out_pw[OUTPUT_VARS::G33] = B * B * r2 * sint * sint;
      base.out_pw[OUTPUT_VARS::K11] = 0.0;
      base.out_pw[OUTPUT_VARS::K12] = 0.0;
      base.out_pw[OUTPUT_VARS::K13] =
          -base.out_pw[OUTPUT_VARS::G33] / 2.0 / N * domega_dr;
      base.out_pw[OUTPUT_VARS::K22] = 0.0;
      base.out_pw[OUTPUT_VARS::K23] =
          -base.out_pw[OUTPUT_VARS::G33] / 2.0 / N * domega_dt;
      base.out_pw[OUTPUT_VARS::K33] = 0.0;

      double const H = base.quant_vals[ISO_VARS::ISO_H];
      double h = std::exp(H);
      double rho, eps, press, vphiU;

      // get quantities point-wise, since h is smoothest, and cut data at H=0
      if (std::fabs(H) <= 1e-14) {
        rho = 0.;
        eps = 0.;
        press = 0.;
        vphiU = 0.;
      } else {
        rho = EOS<eos_t, DENSITY>::get(h);
        eps = EOS<eos_t, EPSILON>::get(h);
        press = EOS<eos_t, PRESSURE>::get(h);
        vphiU = base.quant_vals[ISO_VARS::ISO_U];
      }
      base.out_pw[OUTPUT_VARS::RHO] = rho;
      base.out_pw[OUTPUT_VARS::EPS] = eps;
      base.out_pw[OUTPUT_VARS::PRESS] = press;
      base.out_pw[OUTPUT_VARS::VEL1] = 0.0;
      base.out_pw[OUTPUT_VARS::VEL2] = 0.0;
      base.out_pw[OUTPUT_VARS::VEL3] = vphiU;
      base.output_base = OUTPUT_BASIS::SPHERICAL;
      return base.out_pw;
    }
  };

  /**
   * @brief Interface to export an array of OUTPUT_VARS for a given point in a
   * spherical basis.  The logic for determining the EOS type is here and adds a
   * bit of overhead as compared to a templated class where eos_t is known at
   * compile time. Here we error on the side of convenience rather than speed
   * since the cost is very small.
   *
   * @tparam eos_t C++ Polytrope/Table type
   * @param x
   * @param y
   * @param z
   * @return output_ary_t Interpolated solution at x,y,z
   */
  output_ary_t export_pointwise__spherical(double const& x,
                                           double const& y,
                                           double const& z);
};
}  // namespace Kadath::FUKA_Solvers
