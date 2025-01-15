#pragma once
#include "Solvers/exporter.hpp"
#include "adapted_bh.hpp"
#include <algorithm>
namespace Kadath::FUKA_Solvers {

struct CFMS_BH_Exporter
    : public Exporter<Kadath::FUKA_Config::kadath_config_boost<
                          Kadath::FUKA_Config::BCO_BH_INFO>,
                      Space_adapted_bh> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<
      Kadath::FUKA_Config::BCO_BH_INFO>;
  using space_t = Space_adapted_bh;

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
    NUM_XCTS_VARS
  };

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
    {"k33"  , OUTPUT_VARS::K33}
  };

  std::vector<CFMS_BH_Exporter::XCTS_VARS> xcts_solution_indicies{
    XCTS_PSI,
    XCTS_ALPHA,
    XCTS_BETA1,
    XCTS_BETA2,
    XCTS_BETA3
  };

  std::vector<CFMS_BH_Exporter::XCTS_VARS> xcts_all_indicies{
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
  // clang-format on

  using interp_ary_t = std::array<double, NUM_XCTS_VARS>;
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>;
  using grid_ary_t =
      std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  // Types
  using Exporter<config_t, space_t>::base_space_t;
  using Exporter<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Exporter<config_t, space_t>::space;
  using Exporter<config_t, space_t>::bconfig;
  using Exporter<config_t, space_t>::ndom;

  // CFMS_BH imported fields from file
  ptr_data_member(Scalar, conformal_factor, unique);
  ptr_data_member(Scalar, lapse, unique);
  ptr_data_member(Vector, shift, unique);

  // Constructed objects
  ptr_data_member(Tensor, A, unique);

 protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  interp_ary_t quant_vals;
  interp_ary_t quant_vals_origin;
  output_ary_t out_pw;
  bool export_ready{false};
  int const ndim{3};

  void load_solution_from_file() override;

  void extract_computed_grid_functions();

  void populate_quants();

  public:
  interp_ary_t const & get__quant_vals() const { return quant_vals; }
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  bool is_export_ready() const { return export_ready; }
  void set__export_ready(bool v) { export_ready = v; }
  const int & get_ndim() const { return ndim; }

  CFMS_BH_Exporter()
      : Exporter<config_t, space_t>(),
        conformal_factor(nullptr),
        lapse(nullptr),
        shift(nullptr) {}

  CFMS_BH_Exporter(std::string config_filename)
      : Exporter<config_t, space_t>(config_filename),
        conformal_factor(nullptr),
        lapse(nullptr),
        shift(nullptr) {
    load_solution_from_file();
    populate_quants();
    // Fill outer adapted domain with smooth junk
    export_utils::partial_fill_excision(*this, 1, 2);

    extract_computed_grid_functions();
    populate_quants();

    // Store origin value for use by excision filling of nucleus dom
    this->export_pointwise(0., 0., 0.);
    std::copy(quant_vals.begin(), quant_vals.end(), quant_vals_origin.begin());
    export_ready = true;
  }

  CFMS_BH_Exporter(CFMS_BH_Exporter const& r);
  CFMS_BH_Exporter(CFMS_BH_Exporter&& b) noexcept = delete;
  CFMS_BH_Exporter& operator=(const CFMS_BH_Exporter& b);

 public:
  interp_ary_t interpolate_pointwise_subset(double const & x, double const & y, double const & z,
    std::vector<XCTS_VARS> slice, double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3);
  interp_ary_t interpolate_pointwise__solution_gfs(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3);
  interp_ary_t interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3);

  output_ary_t export_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3);

  grid_ary_t export_coordinate_array(int const npoints,
                                     double const* xx,
                                     double const* yy,
                                     double const* zz,
                                     double const interpolation_offset = 0.,
                                     int const interp_order = 8,
                                     double const delta_r_rel = 0.3);
};
}  // namespace Kadath::FUKA_Solvers
