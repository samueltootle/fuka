#pragma once
#include "Solvers/exporter.hpp"
#include "bin_bh.hpp"

namespace Kadath::FUKA_Solvers {

struct CFMS_BBH_Exporter
    : public Exporter<Kadath::FUKA_Config::kadath_config_boost<
                          Kadath::FUKA_Config::BIN_INFO>,
                      Space_bin_bh> {
  using config_t =
      Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BIN_INFO>;
  using space_t = Space_bin_bh;

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

  // CFMS_BBH imported fields from file
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

  void load_solution_from_file() override;

  void extract_computed_grid_functions();

  void populate_quants();

 public:
  std::vector<std::reference_wrapper<const Scalar>> const& get_quants() const {
    return quants;
  }

  bool is_export_ready() const { return export_ready; }

  const int& get_ndim() const { return ndim; }

  CFMS_BBH_Exporter()
      : Exporter<config_t, space_t>(),
        conformal_factor(nullptr),
        lapse(nullptr),
        shift(nullptr) {}

  CFMS_BBH_Exporter(std::string config_filename)
      : Exporter<config_t, space_t>(config_filename),
        conformal_factor(nullptr),
        lapse(nullptr),
        shift(nullptr) {
    load_solution_from_file();
    extract_computed_grid_functions();
    populate_quants();

    // This is to avoid a "bug" where "something" in kadath is not
    // correctly initialized prior to copying to other threads resulting
    // in undefined behavior.  By running the interpolator once, this
    // bug seems to be avoided.
    this->export_pointwise(0.5, 0., 0.);
  }

  CFMS_BBH_Exporter(CFMS_BBH_Exporter const& r);
  CFMS_BBH_Exporter(CFMS_BBH_Exporter&& b) noexcept = delete;
  CFMS_BBH_Exporter& operator=(const CFMS_BBH_Exporter& b);

 public:
  interp_ary_t interpolate_pointwise(double const& x,
                                     double const& y,
                                     double const& z,
                                     double const interpolation_offset = 0.,
                                     int const interp_order = 8,
                                     double const delta_r_rel = 0.3);

  output_ary_t export_pointwise(double const& x,
                                double const& y,
                                double const& z,
                                double const interpolation_offset = 0.,
                                int const interp_order = 8,
                                double const delta_r_rel = 0.3);

  grid_ary_t export_coordinate_array(int const npoints,
                                     double const* xx,
                                     double const* yy,
                                     double const* zz,
                                     double const interpolation_offset = 0.,
                                     int const interp_order = 8,
                                     double const delta_r_rel = 0.3);
};
}  // namespace Kadath::FUKA_Solvers