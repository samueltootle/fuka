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

struct NS_XCTS_BASE {
  using base_space_t = Space_spheric_adapted;
  using base_config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BCO_NS_INFO>;
  
  ptr_data_member(base_space_t, space, unique);
  ptr_data_member(base_config_t, bconfig, unique);

  using cfgen_t  = CoordFields<base_space_t>;
  using cfary_t  = std::array<std::optional<Vector>, NUM_VECTORS>;

  protected:
  internal_variable(int, rank)
  internal_variable(int, verbosity)
  internal_variable(int, ndom)
  internal_variable(std::string, outputdir)
  ::Kadath::FUKA_Config::STAGES solver_stage{::Kadath::FUKA_Config::STAGES::NUM_STAGES};

  // EOS Parameters - Perhaps this should be a container?
  internal_variable(double, h_cut);
  internal_variable(std::string, eos_file);
  internal_variable(std::string, eos_type);

  // Support containers
  ptr_data_member(Base_tensor, basis, unique);
  ptr_data_member(Metric_flat, fmet, unique);
  ptr_data_member(System_of_eqs, syst, unique);
  ptr_data_member(cfgen_t, cfields, unique);
  ptr_data_member(cfary_t, coord_vectors, unique);

  // Sequence containers
  ptr_data_member(ns_sequence, seq, unique);
  ptr_data_member(Parameter_sequence<BCO_PARAMS>, resolution, unique);

  // Variable fields - i.e. Solution
  ptr_data_member(Scalar, conformal_factor, unique);
  ptr_data_member(Scalar, lapse, unique);
  ptr_data_member(Vector, shift, unique);
  ptr_data_member(Scalar, logh, unique);

  NS_XCTS_BASE();
  NS_XCTS_BASE(base_config_t& config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, 
      int const rank_ = 0);

  protected:
  void initialize_support_containers();
};

template<class eos_t>
struct NS_XCTS_NOROT : NS_XCTS_BASE {

  private:
  void syst_init();
  void print_diagnostics(const int ite, const double conv) const;
  void update_config_quantities();
  void load_solution_from_file();  

  public:
  void save_to_file() const;
  void solve();
  std::string converged_filename(const std::string stage) const;

  NS_XCTS_NOROT() = default;
  NS_XCTS_NOROT(std::string filename);
  NS_XCTS_NOROT(base_config_t& config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, 
      int const rank_ = 0);
};


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
#include "NS_XCTS_BASE.cpp"
#include "NS_XCTS_NOROT.cpp"