#include "EOS/FUKA_EOS_Utilities.hh"
#include "Solvers/exporter.hpp"

namespace Kadath::FUKA_Solvers {

struct CFMS_NS_Exporter
    : public Exporter<Kadath::FUKA_Config::kadath_config_boost<
                          Kadath::FUKA_Config::BCO_NS_INFO>,
                      Space_spheric_adapted> {
    using config_t = Kadath::FUKA_Config::kadath_config_boost<
        Kadath::FUKA_Config::BCO_NS_INFO>;
    using space_t = Space_spheric_adapted;

    // Input ID grid functions that we interpolate on
    enum XCTS_VARS : size_t {
        XCTS_PSI,
        XCTS_ALPHA,
        XCTS_BETA1,
        XCTS_BETA2,
        XCTS_BETA3,
        XCTS_A11,
        XCTS_A12,
        XCTS_A13,
        XCTS_A22,
        XCTS_A23,
        XCTS_A33,
        XCTS_H,
        XCTS_UX,
        XCTS_UY,
        XCTS_UZ,
        NUM_XCTS_VARS
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
  std::vector<XCTS_VARS> xcts_spacetime_indicies{
    XCTS_PSI,
    XCTS_ALPHA,
    XCTS_BETA1,
    XCTS_BETA2,
    XCTS_BETA3,
    XCTS_A11,
    XCTS_A12,
    XCTS_A13,
    XCTS_A22,
    XCTS_A23,
    XCTS_A33
  };
  std::vector<XCTS_VARS> xcts_fluid_indicies{
    XCTS_H,
    XCTS_UX,
    XCTS_UY,
    XCTS_UZ
  };
    // clang-format on

    using interp_ary_t = std::array<double, NUM_XCTS_VARS>;
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
    ptr_data_member(Scalar, conformal_factor, shared);
    ptr_data_member(Scalar, lapse, shared);
    ptr_data_member(Vector, shift, shared);
    ptr_data_member(Scalar, logh, shared);
    ptr_data_member(Scalar, diff_omega, shared);

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
    std::vector<std::reference_wrapper<const Scalar>> const& get_quants()
        const {
        return quants;
    }

    bool is_export_ready() const { return export_ready; }

    const int& get_ndim() const { return ndim; }

    CFMS_NS_Exporter()
        : Exporter<config_t, space_t>(),
          conformal_factor(nullptr),
          lapse(nullptr),
          shift(nullptr),
          logh(nullptr),
          fluidvel(nullptr) {}

    CFMS_NS_Exporter(std::string config_filename)
        : Exporter<config_t, space_t>(config_filename),
          conformal_factor(nullptr),
          lapse(nullptr),
          shift(nullptr),
          logh(nullptr),
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
    CFMS_NS_Exporter(CFMS_NS_Exporter const& r);
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
                                              std::vector<XCTS_VARS> slice);

    /**
   * @brief Export: Load only the spacetime variables
   *
   * @param x
   * @param y
   * @param z
   * @return output_ary_t
   */
    output_ary_t export_pointwise_spacetime_vars(double const& x,
                                                 double const& y,
                                                 double const& z);

    /**
   * @brief Export: Load only the fluid variables
   *
   * @tparam eos_t EOS type
   * @param x
   * @param y
   * @param z
   * @return output_ary_t
   */
    template <class eos_t>
    struct export_pointwise_fluid_vars_imp {
        friend CFMS_NS_Exporter;

        output_ary_t operator()(CFMS_NS_Exporter& base,
                                double const& x,
                                double const& y,
                                double const& z) {
            // Reset to NAN
            base.quant_vals.fill(NAN);
            base.quant_vals =
                base.interpolate_pointwise_subset(x,
                                                  y,
                                                  z,
                                                  base.xcts_fluid_indicies);

            double const H = base.quant_vals[XCTS_VARS::XCTS_H];
            double h = std::exp(H);
            double rho, eps, press;

            // get quantities point-wise, since h is smoothest, and cut data at H=0
            if (std::fabs(H) <= 1e-12) {
                rho = 0.;
                eps = 0.;
                press = 0.;
            } else {
                rho = EOS<eos_t, DENSITY>::get(h);
                eps = EOS<eos_t, EPSILON>::get(h);
                press = EOS<eos_t, PRESSURE>::get(h);
            }
            base.out_pw[OUTPUT_VARS::RHO] = rho;
            base.out_pw[OUTPUT_VARS::EPS] = eps;
            base.out_pw[OUTPUT_VARS::PRESS] = press;
            base.out_pw[OUTPUT_VARS::VEL1] =
                base.quant_vals[XCTS_VARS::XCTS_UX];
            base.out_pw[OUTPUT_VARS::VEL2] =
                base.quant_vals[XCTS_VARS::XCTS_UY];
            base.out_pw[OUTPUT_VARS::VEL3] =
                base.quant_vals[XCTS_VARS::XCTS_UZ];
            return base.out_pw;
        }
    };

    /**
   * @brief Export: Interface for only loading the fluid variables
   *
   * @tparam eos_t EOS type
   * @param x
   * @param y
   * @param z
   * @return output_ary_t
   */
    output_ary_t export_pointwise_fluid_vars(double const& x,
                                             double const& y,
                                             double const& z);

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
        friend CFMS_NS_Exporter;

        output_ary_t operator()(CFMS_NS_Exporter& base,
                                double const& x,
                                double const& y,
                                double const& z) {
            base.quant_vals = base.interpolate_pointwise(x, y, z);

            // Fill output vector by storing non-conformal quantities
            auto const psi = base.quant_vals[XCTS_VARS::XCTS_PSI];
            auto const psi2 = psi * psi;
            auto const psi4 = psi2 * psi2;

            base.out_pw[OUTPUT_VARS::ALPHA] =
                base.quant_vals[XCTS_VARS::XCTS_ALPHA];

            base.out_pw[OUTPUT_VARS::BETA1] =
                base.quant_vals[XCTS_VARS::XCTS_BETA1];
            base.out_pw[OUTPUT_VARS::BETA2] =
                base.quant_vals[XCTS_VARS::XCTS_BETA2];
            base.out_pw[OUTPUT_VARS::BETA3] =
                base.quant_vals[XCTS_VARS::XCTS_BETA3];

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

            base.out_pw[OUTPUT_VARS::G11] = g[0][0];
            base.out_pw[OUTPUT_VARS::G12] = g[0][1];
            base.out_pw[OUTPUT_VARS::G13] = g[0][2];
            base.out_pw[OUTPUT_VARS::G22] = g[1][1];
            base.out_pw[OUTPUT_VARS::G23] = g[1][2];
            base.out_pw[OUTPUT_VARS::G33] = g[2][2];

            base.out_pw[OUTPUT_VARS::K11] =
                base.quant_vals[XCTS_VARS::XCTS_A11] * psi4;
            base.out_pw[OUTPUT_VARS::K12] =
                base.quant_vals[XCTS_VARS::XCTS_A12] * psi4;
            base.out_pw[OUTPUT_VARS::K13] =
                base.quant_vals[XCTS_VARS::XCTS_A13] * psi4;
            base.out_pw[OUTPUT_VARS::K22] =
                base.quant_vals[XCTS_VARS::XCTS_A22] * psi4;
            base.out_pw[OUTPUT_VARS::K23] =
                base.quant_vals[XCTS_VARS::XCTS_A23] * psi4;
            base.out_pw[OUTPUT_VARS::K33] =
                base.quant_vals[XCTS_VARS::XCTS_A33] * psi4;

            double const H = base.quant_vals[XCTS_VARS::XCTS_H];
            double h = std::exp(H);
            double rho, eps, press;

            // get quantities point-wise, since h is smoothest, and cut data at H=0
            if (std::fabs(H) <= 1e-12) {
                rho = 0.;
                eps = 0.;
                press = 0.;
            } else {
                rho = EOS<eos_t, eos_var_t::DENSITY>::get(h);
                eps = EOS<eos_t, eos_var_t::EPSILON>::get(h);
                press = EOS<eos_t, eos_var_t::PRESSURE>::get(h);
            }
            base.out_pw[OUTPUT_VARS::RHO] = rho;
            base.out_pw[OUTPUT_VARS::EPS] = eps;
            base.out_pw[OUTPUT_VARS::PRESS] = press;
            base.out_pw[OUTPUT_VARS::VEL1] =
                base.quant_vals[XCTS_VARS::XCTS_UX];
            base.out_pw[OUTPUT_VARS::VEL2] =
                base.quant_vals[XCTS_VARS::XCTS_UY];
            base.out_pw[OUTPUT_VARS::VEL3] =
                base.quant_vals[XCTS_VARS::XCTS_UZ];
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
};
}  // namespace Kadath::FUKA_Solvers
