#pragma once
#include "adapted_bh.hpp"
#include "bco_utilities.hpp"
#include "codes_utilities.hpp"
#include "exporter_utilities.hpp"
#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "name_tools.hpp"

#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "Configurator/configurator_boost.hpp"

#include "EOS/EOS.hh"

#include "coord_fields.hpp"

#include <cstdlib>
#include <string>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif
#include <mutex>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kadath::FUKA_Solvers {
static std::mutex copy_mutex;

/**
 * @brief The following "Exporter" is really a bandage until the FUKA_Solvers
 * can be rewritten.  Essentially many of the utilities her are duplicate to
 * the solvers, however, pointers are used instead of reference thereby allowing
 * the possibility for dynamic allocation for, e.g. multi-threaded importing of
 * the initial data.  It would have no impact on the solving of the system of
 * equations
 *
 * @tparam config_t
 * @tparam space_t
 */
template <class config_t, class space_t>
struct Exporter {
  using base_config_t = std::decay_t<config_t>;
  using base_space_t = std::decay_t<space_t>;

  ptr_data_member(space_t, space, unique);
  ptr_data_member(base_config_t, bconfig, unique);

 protected:
  int ndom{};

 public:
  Exporter() : space(nullptr) {}

  Exporter(std::string config_filename) : space(nullptr), bconfig(nullptr) {
    bconfig.reset(new base_config_t{config_filename});
    bconfig->open_config();
  }

  virtual void load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen(spacein.c_str(), "r");
    space.reset(new space_t{ff1});
    fclose(ff1);
    ndom = space->get_nbr_domains();
  }
};
}  // namespace Kadath::FUKA_Solvers
