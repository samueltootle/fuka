#include "FUKA_Solvers/NS_XCTS.hpp"
namespace Kadath::FUKA_Solvers {
  void NS_XCTS_DIFFROT::save_to_file() const {
    Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift, *logh, *diff_omega);
  }
  
  NS_XCTS_DIFFROT::NS_XCTS_DIFFROT(std::string filename) {
    bconfig.reset(new base_config_t(filename));
  }

  NS_XCTS_DIFFROT::NS_XCTS_DIFFROT(NS_XCTS_DIFFROT::base_config_t& config_) {
    bconfig.reset(new base_config_t(config_));
  }
}