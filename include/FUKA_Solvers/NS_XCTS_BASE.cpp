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
            resolution(new Parameter_sequence<BCO_PARAMS>(res_)), outputdir(outputdir_) {}

  inline void NS_XCTS_BASE::initialize_support_containers() {

    basis.reset(new Base_tensor(shift->get_basis()));
    fmet.reset(new Metric_flat(*space, *basis));

    cfields.reset(new cfgen_t(*space));
    coord_vectors = std::make_unique<cfary_t>(default_co_vector_ary(*space));
  }
}