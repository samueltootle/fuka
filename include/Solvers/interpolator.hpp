#include <memory>
namespace Kadath::FUKA_Solvers {

double summation_1d (int base, double xx, const Array<double>& tab);

template<class space_t>
struct find_dom {
  int dom{-1};
  find_dom(std::unique_ptr<space_t>& space, Point& p, int dom_min, int dom_max) {
    for(auto d = dom_max-1; d >= dom_min; ++d) {
      if(space->get_domain(d)->is_in(p)) {
        dom = d;
        break;
      }
    }
    if(dom == -1) {
      std::stringstream msg;
      msg << "Point " << p << " not found in the numerical space. ";
      msg << space.get() << ", " << bco_utils::get_radius(space->get_domain(2), INNER_BC) << endl;
      // cout << *(space->get_domain(d)) << endl;
      // throw std::runtime_error(msg.str().c_str());
      cout << msg.str() << endl;
    }
  }
  int operator()() { return dom; }
};

template<class FUKA_READER, class space_ptr_t>
struct Interpolator {
    protected:
    std::unique_ptr<FUKA_READER> solution{nullptr};
    space_ptr_t& space;

    public:
    Interpolator(FUKA_READER const & _sol, space_ptr_t& _space) : space(_space) {
        solution.reset(new FUKA_READER{_sol});
    }

    // DEBUG ONLY
    void print_field_coefs(size_t FIELD_IDX) {
        auto& f = solution->get_field_vector(FIELD_IDX);
        
        uint cnt = 0;
        for(const auto & e : f[2]) {
            std::cout << e << ", ";
            if(cnt > solution->get_ncoefs_i(2, 0)) {
                cout << '\n';
                cnt = 0;
            }
            else cnt++;
        }
        cout << '\n';
    }
    std::vector<double> interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3); 

    double summation (const Point& num, const std::vector<double>& cf) const;
};

template<class FUKA_READER, class space_ptr_t>
std::vector<double> Interpolator<FUKA_READER, space_ptr_t>::interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {
    
    std::vector<double> quant_vals(FUKA_READER::XCTS_VARS::NUM_XCTS_VARS);

    // // Initial guess of the excision radius - needed for filling
    // double rbh = bco_utils::get_radius(space->get_domain(2), INNER_BC);

    // double r2yz = y * y + z * z;

    // double r = std::sqrt(x * x + r2yz);

    // // lambda function for filling excised region
    // auto interp_f = [&](auto& ah_r, auto& extrap_r, auto bh_ori) {
    //   // Avoid division by "0"
    //   extrap_r = (extrap_r <= 1e-12) ? 1e-12 : extrap_r;
    //   double x_shifted = (x == 0.) ? 1e-14 : x;
    //   double theta = std::acos(z / extrap_r);
      
    //   // atan2 is needed here
    //   double phi = std::atan2(y, (x_shifted - bh_ori)); 

    //   // Where the filling takes places
    //   export_utils::spherical_turduck(
    //     quants, quant_vals, interp_order, delta_r_rel, interpolation_offset, 
    //     rbh, extrap_r, theta, phi, 2, bh_ori
    //   );
    // };

    // if (r <= (1. + interpolation_offset) * rbh) {
    //   interp_f(rbh, r, 0.);
    //   // For testing only
    //   // quant_vals[XCTS_VARS::XCTS_ALPHA] = fd();
    // } else { 
      Point abs_coords(solution.get_ndim());
      abs_coords.set(1) = x;
      abs_coords.set(2) = y;
      abs_coords.set(3) = z;
      find_dom fd(space, abs_coords, 0, solution.get_ndom());
      Point num(space->get_domain(fd())->absol_to_num(abs_coords)) ;
    //   cout << "Dom: " << fd() << endl;
      // For testing only
      // quant_vals[XCTS_VARS::XCTS_ALPHA] = fd();
    //   for (int k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
    //     quant_vals[k] = quants[k].get().val_point(abs_coords);
    //   }
    // }
    return quant_vals;
  }

template<class FUKA_READER, class space_ptr_t>
double Interpolator<FUKA_READER, space_ptr_t>::summation (const Point& num, const std::vector<double>& cf) const {

    std::vector<double>* courant = new std::vector<double>(cf) ;
    auto const & ndim = solution.get_ndim();
    // Dim_array nbr_coefs (cf.get_dimensions()) ;
    
    // Loop on dimensions (except the last):
    for (int d=0 ; d<ndim-1 ; d++) {
        int dim_output = ndim-1-d ;
    //     Dim_array nbr_output (dim_output) ;
    //     for (int k=0 ; k<dim_output ; k++)
    //         nbr_output.set(k) = nbr_coefs(k+d+1) ;
    //     Array<double> output (nbr_output) ;
        
    //     Index inout (nbr_output) ;
    //     Array<double> tab_1d (cf.get_size(d)) ;
    //     Index incourant (courant->get_dimensions()) ;
        
    //     // Loop on the points :
    //     bool loop = true ;
    //     while (loop) {
    //         int base = (*bases_1d[d])(inout) ;
            
    //         for (int k=0 ; k<dim_output ; k++)
    //             incourant.set(k+1) = inout(k) ;
    //         for (int k=0 ; k<cf.get_size(d) ; k++) {
    //             incourant.set(0) = k ;
    //             tab_1d.set(k) = (*courant)(incourant) ;
    //         }
    //         output.set(inout) = summation_1d (base, num(d+1), tab_1d) ;
    //         loop = inout.inc() ;
    //         }
    //     delete courant ;
    //     courant = new Array<double> (output) ;
    }
    
    // // Last base :
    // assert (courant->get_ndim()==1) ;
    // assert (bases_1d[ndim-1]->get_ndim()==1) ;
    // assert (bases_1d[ndim-1]->get_size(0)==1) ;
    // double res = summation_1d ((*bases_1d[ndim-1])(0), num(ndim), *courant) ;
    // delete courant ;
    double res = 0;
    return res ;
}
}