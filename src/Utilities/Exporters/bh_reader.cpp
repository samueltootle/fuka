#include "Solvers/bh_3d_xcts/bh_reader.hpp"

namespace Kadath::FUKA_Solvers {
double summation_1d (int base, double xx, const Array<double>& tab);

void CFMS_BH_Readerv2::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space_t _space{ff1};
    Scalar _conf(_space, ff1) ;
    Scalar _lapse( _space, ff1) ;
    Vector _shift( _space, ff1) ;
    
    fclose(ff1);
    _conf.coef();
    _lapse.coef();
    _shift.coef();

    ndom = _space.get_nbr_domains();
    ncoefs.resize(ndom);
    ncoefs1d.resize(ndom);
  
    for(auto& V : id_vars) {
      V.resize(ndom);
    }

    for(auto d = 0; d < ndom; ++d) {
      auto nbr_coefs = _space.get_domain(d)->get_nbr_coefs();
      ncoefs[d][0] = nbr_coefs(0); 
      ncoefs[d][1] = nbr_coefs(1); 
      ncoefs[d][2] = nbr_coefs(2);
      ncoefs1d[d] = nbr_coefs(0) * nbr_coefs(1) * nbr_coefs(2);
      
      for(auto& V : id_vars) {
        V[d].resize(ncoefs1d[d]);
      }
      
      Array<double> const _conf_d(_conf(d).get_coef());
      Array<double> const _lapse_d(_lapse(d).get_coef());
      Array<double> const _shift1_d(_shift(1)(d).get_coef());
      Array<double> const _shift2_d(_shift(2)(d).get_coef());
      Array<double> const _shift3_d(_shift(3)(d).get_coef());
      Index I(_conf_d.get_dimensions());

      bool loop=true;
      while(loop) {
        const auto i = I(0);
        const auto j = I(1);
        const auto k = I(2);
        auto idx = IDX(d, i, j, k);
        
        
        id_vars[XCTS_PSI][d][idx] = _conf_d(I);
        id_vars[XCTS_ALPHA][d][idx] = _lapse_d(I);
        id_vars[XCTS_BETAX][d][idx] = _shift1_d(I);
        id_vars[XCTS_BETAY][d][idx] = _shift2_d(I);
        id_vars[XCTS_BETAZ][d][idx] = _shift3_d(I);
        loop = I.inc();
      }
    }
    
    export_ready = true;
    // basis.reset(new Base_tensor{shift->get_basis()});
    // fmet.reset(new Metric_flat(*space, *basis));
  }
}