#pragma once
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"
#include "EOS/EOS.hh"
#include "name_tools.hpp"
#include "exporter_utilities.hpp"
#include "Solvers/fuka_syst/fuka_syst_setup.hpp"
#include "Solvers/fuka_syst/fuka_syst_tools.hpp"
#include "codes_utilities.hpp"
#include "bco_utilities.hpp"
#include "exporter_utilities.hpp"

#include <cstdlib>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

namespace Kadath::FUKA_Solvers {
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
  ptr_data_member(System_of_eqs, syst, unique);
  ptr_data_member(base_config_t, bconfig, unique);
  
  protected:
  int ndom{};

  public:
  Reader() : space(nullptr), syst(nullptr) {}
  Reader(std::string config_filename) : 
    space(nullptr), syst(nullptr), bconfig(nullptr) {
      bconfig.reset(new base_config_t{config_filename});
      bconfig->open_config();
  }
  
  virtual void load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space.reset(new space_t{ff1});
    syst.reset(new System_of_eqs{*space});
    fclose(ff1);
    ndom = space->get_nbr_domains();
  }
};

template<class config_t, class space_t>
struct CFMS_BH_Reader : public Reader<config_t, space_t> {
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
  constexpr static size_t nout = size_t(OUTPUT_VARS::NUM_OUTPUT_VARS);

  // Types
  using Reader<config_t, space_t>::base_space_t;
  using Reader<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Reader<config_t, space_t>::space;
  using Reader<config_t, space_t>::bconfig;
  using Reader<config_t, space_t>::syst;
  using Reader<config_t, space_t>::ndom;

  // CFMS_BH specific fields
  ptr_data_member(Scalar, conformal_factor, unique);
  ptr_data_member(Scalar, lapse, unique);
  ptr_data_member(Vector, shift, unique);
  ptr_data_member(Base_tensor, basis, unique);
  ptr_data_member(Metric_flat, fmet, unique);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  bool export_ready{false};

  void load_solution_from_file() override {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space.reset(new space_t{ff1});
    conformal_factor.reset( new Scalar(*space.get(), ff1)) ;
    lapse.reset( new Scalar(*space.get(), ff1)) ;
    shift.reset( new Vector(*space.get(), ff1)) ;
    fclose(ff1);
    basis.reset(new Base_tensor{shift->get_basis()});
    fmet.reset(new Metric_flat(*space, *basis));
    ndom = space->get_nbr_domains();
  }

  void extract_xcts_grid_functions() {
    if(quants.capacity() != XCTS_VARS::NUM_XCTS_VARS) {
      for (int i = 0; i < XCTS_VARS::NUM_XCTS_VARS; ++i)
        quants.push_back(std::cref(*conformal_factor));
    }
    syst.reset(new System_of_eqs (*space));    
    fmet->set_system(*syst, "f") ;

    // Fields - must be initialized before common setup
    syst->add_cst("N"  , *lapse) ;
    syst->add_cst("bet", *shift) ;

    // syst->add_def("Ts^i = N * bet^i");
    syst->add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    Tensor A(syst->give_val_def("A"));

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
      XCTS_VARS::XCTS_AZZ}, A);
    export_ready = true;
  }

  public:
  bool is_export_ready() const { return export_ready; }

  CFMS_BH_Reader(std::string config_filename) :
    Reader<config_t, space_t>(config_filename),
      basis(nullptr), fmet(nullptr),
        conformal_factor(nullptr), lapse(nullptr), shift(nullptr) {
  
    load_solution_from_file();
    extract_xcts_grid_functions();
  }

  CFMS_BH_Reader(CFMS_BH_Reader const & r) : fmet(nullptr) {
    space.reset(new space_t((*r.get_space())));
    basis.reset(new Base_tensor(*space, r.get_basis()->get_basis(0)));
    fmet.reset(new Metric_flat(*space, *basis));
    bconfig.reset(new config_t(*r.bconfig));
    
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    bconfig->set_filename("test");
    export_ready = false;

    extract_xcts_grid_functions();
    
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
  }

  

  public:
  std::array<std::vector<double>, nout> export_array(int const npoints, double const * xx, double const * yy, double const * zz,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {

    //sdf

  }
};
}