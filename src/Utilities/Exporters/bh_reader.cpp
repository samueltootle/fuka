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
    ncoefs.resize(ndom);
    conformal_factor.resize(ndom);

    for(auto d = 0; d < ndom; ++d) {
        auto nbr_coefs = _space.get_domain(d)->get_nbr_coefs();
        ncoefs[d][0] = nbr_coefs(0); 
        ncoefs[d][1] = nbr_coefs(1); 
        ncoefs[d][2] = nbr_coefs(2);
        ncoefs1d[d] = nbr_coefs(0) * nbr_coefs(1) * nbr_coefs(2);
        conformal_factor[d].resize(ncoefs1d[d]);
    }
    is_ready = true;
    // basis.reset(new Base_tensor{shift->get_basis()});
    // fmet.reset(new Metric_flat(*space, *basis));
  }
}