#include "utilities.hpp"
namespace Kadath::FUKA_Solvers {

  inline NS_XCTS_BASE::NS_XCTS_BASE() :
      rank(0.), verbosity(0), ndom(-1), bconfig(nullptr), basis(nullptr), fmet(nullptr), 
        cfields(nullptr), coord_vectors(nullptr), seq(nullptr), 
          resolution(nullptr), syst(nullptr), conformal_factor(nullptr), 
            lapse(nullptr), shift(nullptr), logh(nullptr), outputdir("") {}

  inline NS_XCTS_BASE::NS_XCTS_BASE(NS_XCTS_BASE::base_config_t& config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, int const rank_) :
      rank(rank_), verbosity(0), ndom(-1), bconfig(new base_config_t(config_)), basis(nullptr), fmet(nullptr), 
        cfields(nullptr), coord_vectors(nullptr), syst(nullptr), conformal_factor(nullptr), 
          lapse(nullptr), shift(nullptr), logh(nullptr), seq(new ns_sequence(seq_)), 
            resolution(new Parameter_sequence<BCO_PARAMS>(res_)), outputdir(outputdir_) {
    
    std::array<bool, NUM_STAGES>& stage_enabled = bconfig->return_stages();
    auto [ last_stage_, last_stage_idx_ ] = get_last_enabled(MSTAGE, stage_enabled);
    last_stage_idx = last_stage_idx_;
  }

  inline void NS_XCTS_BASE::initialize_support_containers() {

    basis.reset(new Base_tensor(shift->get_basis()));
    fmet.reset(new Metric_flat(*space, *basis));

    cfields.reset(new cfgen_t(*space));
    coord_vectors = std::make_unique<cfary_t>(default_co_vector_ary(*space));
  }

  inline void NS_XCTS_BASE::save_to_file() const {
    Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift, *logh);
  }

  inline void NS_XCTS_BASE::reset_all_ptrs() {
    // Fields
    conformal_factor.reset(nullptr) ;
    lapse.reset(nullptr);
    shift.reset(nullptr);
    logh.reset(nullptr);

    // Containers
    basis.reset(nullptr);
    fmet.reset(nullptr);
    syst.reset(nullptr);
    cfields.reset(nullptr);
    coord_vectors.reset(nullptr);

    // Space
    space.reset(nullptr);
  }
}