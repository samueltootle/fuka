#include "Solvers/bh_3d_xcts/bh_reader.hpp"

namespace Kadath::FUKA_Solvers {

void CFMS_BH_Readerv2::initialize_containers() {
    ncoefs.resize(ndom);
    ncoefs1d.resize(ndom);
  
    for(auto& V : id_vars) {
      V.resize(ndom);
    }
}

void CFMS_BH_Readerv2::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    space_t _space{ff1};
    Scalar _conf(_space, ff1) ;
    Scalar _lapse( _space, ff1) ;
    Vector _shift( _space, ff1) ;
    fclose(ff1);
    std::unique_ptr<Tensor> A{nullptr};
    {
      Base_tensor basis(Base_tensor{_shift.get_basis()});
      Metric_flat fmet(_space, basis);
      System_of_eqs syst(_space);
      fmet.set_system(syst, "f");

      syst.add_cst("N"  , _lapse) ;
      syst.add_cst("bet", _shift) ;

      // syst->add_def("Ts^i = N * bet^i");
      syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
      A.reset(new Tensor(syst.give_val_def("A")));
    }
    _conf.coef();
    _lapse.coef();
    _shift.coef();
    A->coef();

    ndom = _space.get_nbr_domains();
    initialize_containers();

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
        {
          uint c = 0;
          const auto A_IDXS = {
            XCTS_VARS::XCTS_AXX, 
            XCTS_VARS::XCTS_AXY, 
            XCTS_VARS::XCTS_AXZ, 
            XCTS_VARS::XCTS_AYY, 
            XCTS_VARS::XCTS_AYZ, 
            XCTS_VARS::XCTS_AZZ};
          for(auto i : A_IDXS) {
              auto tidx = export_utils::R2TensorSymmetricIndices[c];
              Kadath::Array<int> ind (A->indices(tidx));
              // quants[i] = std::cref(field(ind));
              Array<double> const A_coef = (*A)(ind)(d).get_coef();
              id_vars[i][d][idx] = A_coef(I);
            }
        }
        loop = I.inc();
      }
    }
    export_ready = true;
  }
}