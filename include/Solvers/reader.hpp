#pragma once
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"
#include "EOS/EOS.hh"
#include "name_tools.hpp"
#include "exporter_utilities.hpp"
#include <cstdlib>
#include <string>
#include <filesystem>
#include "codes_utilities.hpp"
#include "bco_utilities.hpp"
namespace fs = std::filesystem;
using namespace export_utils;

namespace Kadath::FUKA_Solvers {

template<class config_t, class space_t>
struct Reader {
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;

  ptr_data_member(space_t, space, unique);
  ptr_data_member(System_of_eqs, syst, unique);
  ptr_data_member(base_config_t, bconfig, unique);

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
  }
};

template<class config_t, class space_t>
struct CFMS_BH_Reader : public Reader<config_t, space_t> {
  // Types
  using Reader<config_t, space_t>::base_space_t;
  using Reader<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Reader<config_t, space_t>::space;
  using Reader<config_t, space_t>::bconfig;
  using Reader<config_t, space_t>::syst;

  // CFMS_BH specific fields
  ptr_data_member(Scalar, conformal_factor, unique);
  ptr_data_member(Scalar, lapse, unique);
  ptr_data_member(Vector, shift, unique);
  ptr_data_member(Base_tensor, basis, unique);

  CFMS_BH_Reader(std::string config_filename) :
    Reader<config_t, space_t>(config_filename),
      basis(nullptr),
        conformal_factor(nullptr), lapse(nullptr), shift(nullptr) {
  
    load_solution_from_file();
  }

  CFMS_BH_Reader(CFMS_BH_Reader const & r) {
    space.reset(new space_t((*r.get_space())));
    basis.reset(new Base_tensor(*space, r.get_basis()->get_basis(0)));
    syst.reset(new System_of_eqs{*space});
    bconfig.reset(new config_t(*r.bconfig));
    
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    bconfig->set_filename("test");
    Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
  }

  void load_solution_from_file() override {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space.reset(new space_t{ff1});
    syst.reset(new System_of_eqs{*space});
    conformal_factor.reset( new Scalar(*space.get(), ff1)) ;
    lapse.reset( new Scalar(*space.get(), ff1)) ;
    shift.reset( new Vector(*space.get(), ff1)) ;
    fclose(ff1);
    basis.reset(new Base_tensor{shift->get_basis()});
  }

};
}