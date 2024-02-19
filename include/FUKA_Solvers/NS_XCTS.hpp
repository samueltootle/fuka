#include "kadath.hpp"
#include "kadath_adapted.hpp"
#include "codes_utilities.hpp"

#include "Configurator/config_bco.hpp"
#include "coord_fields.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"

#include "Solvers/sequences/parameter_sequence.hpp"
#include "Solvers/sequences/ns_sequence.hpp"

/**
 * \addtogroup NS_XCTS
 * \ingroup FUKA
 * @{*/
namespace Kadath::FUKA_Solvers {
// using ::Kadath::FUKA_Config;
// using ::Kadath::FUKA_Config_Utils;

struct NS_XCTS_DIFFROT {
  using base_space_t = Space_spheric_adapted;
  using base_config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_NS_INFO>;

  using cfgen_t  = CoordFields<base_space_t>;
  using cfary_t  = std::array<std::optional<Vector>, NUM_VECTORS>;
  
  ptr_data_member(base_space_t, space, unique);
  ptr_data_member(base_config_t, bconfig, unique);

  protected:
  static constexpr int ndom{3};
  int verbosity;
  ::Kadath::FUKA_Config::STAGES solver_stage;

  // Support containers
  ptr_data_member(Base_tensor, basis, unique);
  ptr_data_member(Metric_flat, fmet, unique);
  ptr_data_member(cfgen_t, cfields, unique);
  ptr_data_member(cfary_t, coord_vectors, unique);
  ptr_data_member(ns_sequence, seq, unique);

  // Variable fields - i.e. Solution
  ptr_data_member(Scalar, conformal_factor, unique);
  ptr_data_member(Scalar, lapse, unique);
  ptr_data_member(Vector, shift, unique);
  ptr_data_member(Scalar, logh, unique);
  ptr_data_member(Scalar, diff_omega, unique);

  public:
  void save_to_file() const;

  NS_XCTS_DIFFROT() = default;
  NS_XCTS_DIFFROT(std::string filename);
  NS_XCTS_DIFFROT(base_config_t& config_);
};
/** @}*/
};