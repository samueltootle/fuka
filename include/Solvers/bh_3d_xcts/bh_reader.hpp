#pragma once
#include "Solvers/reader.hpp"
#include "adapted_bh.hpp"
namespace Kadath::FUKA_Solvers {
using namespace Kadath::FUKA_Config;
template<class config_t, class space_t>
struct Readerv2 {
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;

  ptr_data_member(base_config_t, bconfig, unique);
  
  protected:
  int ndom{};

  public:
  Readerv2() {}
  Readerv2(std::string config_filename) : 
    bconfig(nullptr) {
      bconfig.reset(new base_config_t{config_filename});
      bconfig->open_config();
  }
};

// template<class config_t = , class Space_adapted_bh>
struct CFMS_BH_Readerv2 : public Readerv2<kadath_config_boost<BCO_BH_INFO>, Space_adapted_bh> {
  using config_t = kadath_config_boost<BCO_BH_INFO>;
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
  using pointwise_ary_t = std::vector<double>; 
  using grid_ary_t = std::array<pointwise_ary_t, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  using Readerv2<config_t, space_t>::base_config_t;
  using Readerv2<config_t, space_t>::bconfig;
  
  protected:
  using Readerv2<config_t, space_t>::ndom;
  std::vector<std::array<uint, 3>> ncoefs{};
  std::vector<uint> ncoefs1d{};
  std::vector<std::reference_wrapper<const Scalar>> quants;
  bool export_ready{false};
  int const ndim{3};

  std::vector<std::vector<double>> conformal_factor{};
  std::vector<std::vector<double>> lapse{};

  void load_solution_from_file();
  
  // super index for a given domain across all coefficients
  size_t IDX(size_t d, size_t r, size_t theta, size_t phi) {
    return ((r) + ncoefs[d][0] * ((theta) + ncoefs[d][1] * ((phi))))
  }

  public:
  bool is_export_ready() const { return export_ready; }
  CFMS_BH_Readerv2() : Readerv2<config_t, space_t>() {}
  CFMS_BH_Readerv2(std::string config_filename) :
    Readerv2<config_t, space_t>(config_filename) {
    load_solution_from_file();
  }
};
}