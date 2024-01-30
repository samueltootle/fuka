#include "Solvers/bh_3d_xcts/bh_exporter.hpp"
namespace Kadath::FUKA_Solvers {
  CFMS_BH_Exporter::CFMS_BH_Exporter(CFMS_BH_Exporter const & r) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;

    space.reset(new space_t((*r.get_space())));
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    A.reset(new Tensor(*space, *r.A.get()));

    bconfig.reset(new config_t(*r.bconfig));
    
    export_ready = false;

    populate_quants();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
    // std::cout << "copy\n";
  }
  
  CFMS_BH_Exporter& CFMS_BH_Exporter::operator=(const CFMS_BH_Exporter& b) {
    if (this == &b) return *this;

    CFMS_BH_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
  }
  
  void CFMS_BH_Exporter::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;
    
    space.reset(new space_t{ff1});
    conformal_factor.reset( new Scalar(*space.get(), ff1)) ;
    lapse.reset( new Scalar(*space.get(), ff1)) ;
    shift.reset( new Vector(*space.get(), ff1)) ;
    
    fclose(ff1);

    ndom = space->get_nbr_domains();
  }

  void CFMS_BH_Exporter::extract_computed_grid_functions() {
    System_of_eqs syst(*space);
    Base_tensor basis(shift->get_basis());
    Metric_flat fmet(*space, basis);
    fmet.set_system(syst, "f") ;

    // Fields - must be initialized before common setup
    syst.add_cst("N"  , *lapse) ;
    syst.add_cst("bet", *shift) ;

    syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    A.reset(new Tensor(syst.give_val_def("A")));
    A->coef();
  }

  void CFMS_BH_Exporter::populate_quants() {
    if(quants.capacity() != XCTS_VARS::NUM_XCTS_VARS) {
      for (size_t i = 0; i < XCTS_VARS::NUM_XCTS_VARS; ++i)
        quants.push_back(std::cref(*conformal_factor));
    }
    quants[XCTS_VARS::XCTS_PSI] = std::cref(*conformal_factor);
    quants[XCTS_VARS::XCTS_ALPHA] = std::cref(*lapse);
    quants[XCTS_VARS::XCTS_BETA1] = std::cref((*shift)(1));
    quants[XCTS_VARS::XCTS_BETA2] = std::cref((*shift)(2));
    quants[XCTS_VARS::XCTS_BETA3] = std::cref((*shift)(3));

    export_utils::add_tensor_refs(quants, {
      XCTS_VARS::XCTS_A11, 
      XCTS_VARS::XCTS_A12, 
      XCTS_VARS::XCTS_A13, 
      XCTS_VARS::XCTS_A22, 
      XCTS_VARS::XCTS_A23, 
      XCTS_VARS::XCTS_A33}, *A);
    export_ready = true;
  }

  CFMS_BH_Exporter::interp_ary_t CFMS_BH_Exporter::interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {

    // Initial guess of the excision radius - needed for filling
    // FIXME make excision generic
    double rbh = bco_utils::get_radius(space->get_domain(2), INNER_BC);

    double r2yz = y * y + z * z;

    double r = std::sqrt(x * x + r2yz);

    // lambda function for filling excised region
    auto interp_f = [&](auto& ah_r, auto& extrap_r, auto bh_ori) {
      // Avoid division by "0"
      extrap_r = (extrap_r <= 1e-12) ? 1e-12 : extrap_r;
      double x_shifted = (x == 0.) ? 1e-14 : x;
      double theta = std::acos(z / extrap_r);
      
      // atan2 is needed here
      double phi = std::atan2(y, (x_shifted - bh_ori)); 

      // Where the filling takes places
      export_utils::spherical_turduck(
        quants, quant_vals, interp_order, delta_r_rel, interpolation_offset, 
        ah_r, extrap_r, theta, phi, 2, bh_ori
      );
    };

    if (r <= (1. + interpolation_offset) * rbh) {
      interp_f(rbh, r, 0.);
    } else { 
      Point abs_coords(ndim);
      abs_coords.set(1) = x;
      abs_coords.set(2) = y;
      abs_coords.set(3) = z;

      for (size_t k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
      }
    }
    return quant_vals;
  }

  CFMS_BH_Exporter::output_ary_t CFMS_BH_Exporter::export_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {
      
    quant_vals = interpolate_pointwise(x, y, z, interpolation_offset, interp_order, delta_r_rel);
    
    // Fill output vector by storing non-conformal quantities
    auto const psi = quant_vals[XCTS_VARS::XCTS_PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    out_pw[OUTPUT_VARS::ALPHA] = quant_vals[XCTS_VARS::XCTS_ALPHA];

    out_pw[OUTPUT_VARS::BETA1] = quant_vals[XCTS_VARS::XCTS_BETA1];
    out_pw[OUTPUT_VARS::BETA2] = quant_vals[XCTS_VARS::XCTS_BETA2];
    out_pw[OUTPUT_VARS::BETA3] = quant_vals[XCTS_VARS::XCTS_BETA3];

    double g[3][3];
    g[0][0] = psi4;
    g[0][1] = 0.0;
    g[0][2] = 0.0;
    g[1][1] = psi4;
    g[1][2] = 0.0;
    g[2][2] = psi4;
    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];

    out_pw[OUTPUT_VARS::G11] = g[0][0];
    out_pw[OUTPUT_VARS::G12] = g[0][1];
    out_pw[OUTPUT_VARS::G13] = g[0][2];
    out_pw[OUTPUT_VARS::G22] = g[1][1];
    out_pw[OUTPUT_VARS::G23] = g[1][2];
    out_pw[OUTPUT_VARS::G33] = g[2][2];

    out_pw[OUTPUT_VARS::K11] = quant_vals[XCTS_VARS::XCTS_A11] * psi4;
    out_pw[OUTPUT_VARS::K12] = quant_vals[XCTS_VARS::XCTS_A12] * psi4;
    out_pw[OUTPUT_VARS::K13] = quant_vals[XCTS_VARS::XCTS_A13] * psi4;
    out_pw[OUTPUT_VARS::K22] = quant_vals[XCTS_VARS::XCTS_A22] * psi4;
    out_pw[OUTPUT_VARS::K23] = quant_vals[XCTS_VARS::XCTS_A23] * psi4;
    out_pw[OUTPUT_VARS::K33] = quant_vals[XCTS_VARS::XCTS_A33] * psi4;
    return out_pw;
  }

  CFMS_BH_Exporter::grid_ary_t CFMS_BH_Exporter::export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {
    
    grid_ary_t out;
    for(auto& v : out) {
      v.resize(npoints);
    }
    
    for (int i = 0; i < npoints; ++i) {
      export_pointwise(xx[i], yy[i], zz[i], interpolation_offset, interp_order, delta_r_rel);
      
      out[OUTPUT_VARS::ALPHA][i] = out_pw[OUTPUT_VARS::ALPHA];

      out[OUTPUT_VARS::BETA1][i] = out_pw[OUTPUT_VARS::BETA1];
      out[OUTPUT_VARS::BETA2][i] = out_pw[OUTPUT_VARS::BETA2];
      out[OUTPUT_VARS::BETA3][i] = out_pw[OUTPUT_VARS::BETA3];

      out[OUTPUT_VARS::G11][i] = out_pw[OUTPUT_VARS::G11];
      out[OUTPUT_VARS::G12][i] = out_pw[OUTPUT_VARS::G12];
      out[OUTPUT_VARS::G13][i] = out_pw[OUTPUT_VARS::G13];
      out[OUTPUT_VARS::G22][i] = out_pw[OUTPUT_VARS::G22];
      out[OUTPUT_VARS::G23][i] = out_pw[OUTPUT_VARS::G23];
      out[OUTPUT_VARS::G33][i] = out_pw[OUTPUT_VARS::G33];

      out[OUTPUT_VARS::K11][i] = out_pw[OUTPUT_VARS::K11];
      out[OUTPUT_VARS::K12][i] = out_pw[OUTPUT_VARS::K12];
      out[OUTPUT_VARS::K13][i] = out_pw[OUTPUT_VARS::K13];
      out[OUTPUT_VARS::K22][i] = out_pw[OUTPUT_VARS::K22];
      out[OUTPUT_VARS::K23][i] = out_pw[OUTPUT_VARS::K23];
      out[OUTPUT_VARS::K33][i] = out_pw[OUTPUT_VARS::K33];
    }
    return out;
  }
}