#include "Solvers/bhns_xcts/bhns_exporter.hpp"
namespace Kadath::FUKA_Solvers {
  CFMS_BHNS_Exporter::CFMS_BHNS_Exporter(CFMS_BHNS_Exporter const & r) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;
    h_cut= r.h_cut;
    eos_file = r.eos_file;
    eos_type = r.eos_type;

    space.reset(new space_t((*r.get_space())));
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    logh.reset(new Scalar(*space, *r.logh.get()));
    velpotential.reset(new Scalar(*space, *r.velpotential.get()));
    fluidvel.reset(new Vector(*space, *r.fluidvel.get()));
    A.reset(new Tensor(*space, *r.A.get()));

    bconfig.reset(new config_t(*r.bconfig));
    
    export_ready = false;

    populate_quants();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
    // std::cout << "copy\n";
  }
  
  CFMS_BHNS_Exporter& CFMS_BHNS_Exporter::operator=(const CFMS_BHNS_Exporter& b) {
    if (this == &b) return *this;

    CFMS_BHNS_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
  }
  
  void CFMS_BHNS_Exporter::initialize_eos() {
    // Initialize EOS
    h_cut    = (*bconfig).template eos<double>(Kadath::FUKA_Config::EOS_PARAMS::HCUT, Kadath::FUKA_Config::NODES::BCO1);
    eos_file = (*bconfig).template eos<std::string>(Kadath::FUKA_Config::EOS_PARAMS::EOSFILE, Kadath::FUKA_Config::NODES::BCO1);
    eos_type = (*bconfig).template eos<std::string>(Kadath::FUKA_Config::EOS_PARAMS::EOSTYPE, Kadath::FUKA_Config::NODES::BCO1);

    if(eos_type == "Cold_Table") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = ((*bconfig).template eos<int>(Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS, Kadath::FUKA_Config::NODES::BCO1) == 0) ? \
            2000 : (*bconfig).template eos<int>(Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS, Kadath::FUKA_Config::NODES::BCO1);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    } else if(eos_type == "Cold_PWPoly") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    }  else {
      throw std::invalid_argument("\nInvalid EOS Type\n)");
    }
    // end adding EOS OPEs
  }
  
  void CFMS_BHNS_Exporter::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;

    space.reset(new space_t{ff1});
    conformal_factor.reset( new Scalar(*space.get(), ff1)) ;
    lapse.reset( new Scalar(*space.get(), ff1)) ;
    shift.reset( new Vector(*space.get(), ff1)) ;
    logh.reset( new Scalar(*space.get(), ff1)) ;
    velpotential.reset( new Scalar(*space.get(), ff1)) ;
    
    fclose(ff1);
    
    ndom = space->get_nbr_domains();
  }

  void CFMS_BHNS_Exporter::extract_computed_grid_functions() {
    // fields depending on the coords
    CoordFields<space_t> cfields(*space);
    vec_ary_t coord_vectors {default_binary_vector_ary(*space)};

    // get origin of the system and initialize coordinate fields
    double xc1 = Kadath::bco_utils::get_center(*space,space->NS);
    double xc2 = Kadath::bco_utils::get_center(*space,space->BH);
    double xo  = Kadath::bco_utils::get_center(*space,ndom-1);
    update_fields(cfields, coord_vectors, {}, xo, xc1, xc2);

    // Initialize Flat Metric
    Base_tensor basis(shift->get_basis());
    Metric_flat fmet(*space, basis);

    // Start - Setup System of equations
    System_of_eqs syst(*space);
    fmet.set_system(syst, "f") ;

    Param p;
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      set_eos_ope<eos_t>(syst, p);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      set_eos_ope<eos_t>(syst, p);
    } // end adding EOS OPEs

    // Fields - must be initialized before common setup
    syst.add_cst("P"  , *conformal_factor);
    syst.add_cst("N"  , *lapse) ;
    syst.add_cst("bet", *shift) ;
    syst.add_cst("H"  , *logh);
    syst.add_cst("phi", *velpotential);

    syst.add_cst("omes1" , (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA, Kadath::FUKA_Config::NODES::BCO1));
    syst.add_cst("omes2" , (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA, Kadath::FUKA_Config::NODES::BCO2));
    
    syst.add_cst("mg", *coord_vectors[Kadath::coord_vector::GLOBAL_ROT]);
    syst.add_cst("mm", *coord_vectors[Kadath::coord_vector::BCO1_ROT]) ;
    syst.add_cst("mp", *coord_vectors[Kadath::coord_vector::BCO2_ROT]) ;

    syst.add_def("h = exp(H)");
    for(int d = space->NS; d <= space->ADAPTEDNS; ++d){
      syst.add_def(d, "s^i  = omes1 * mm^i");
    }
    for(int d = 0; d < ndom; ++d) {
      if(d <= space->ADAPTEDNS)
        syst.add_def(d, "eta_i = D_i phi + P^4 * s_i");
      else
        syst.add_def(d, "eta_i = D_i phi");
    }

    syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    A.reset(new Tensor(syst.give_val_def("A")));
    A->coef();
  
    // definitions for the fluid 3-velocity
    syst.add_def("Wsquare = eta^i * eta_i / h^2 / P^4 + 1.");
    syst.add_def("W = sqrt(Wsquare)");

    syst.add_def("U^i = eta^i / P^4 / h / W");
    fluidvel.reset(new Vector(syst.give_val_def("U")));
    fluidvel->coef();
  }


  void CFMS_BHNS_Exporter::populate_quants() {
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
    
    // Fluid related quantities
    quants[XCTS_VARS::XCTS_H] = std::cref(*logh);
    quants[XCTS_VARS::XCTS_UX] = std::cref((*fluidvel)(1));
    quants[XCTS_VARS::XCTS_UY] = std::cref((*fluidvel)(2));
    quants[XCTS_VARS::XCTS_UZ] = std::cref((*fluidvel)(3));
    export_ready = true;
  }

  CFMS_BHNS_Exporter::interp_ary_t CFMS_BHNS_Exporter::interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {

    double const xBH = Kadath::bco_utils::get_center(*space,space->BH);
    double const rBH = Kadath::bco_utils::get_radius(space->get_domain(space->ADAPTEDBH + 1),INNER_BC);
    
    double const & xcom_shift = (*bconfig)(Kadath::FUKA_Config::BIN_PARAMS::COM);
    double const & ycom_shift = (*bconfig)(Kadath::FUKA_Config::BIN_PARAMS::COMY);
    
    double const x_shifted = x - xcom_shift;
    double const y_shifted = y - ycom_shift;

    // Relative x-coordinate to BH center
    double const rel_xBH = x_shifted - xBH;

    double r2yz = y_shifted * y_shifted + z * z;

    // Relative radius to BH center
    double rel_rBH = std::sqrt(rel_xBH * rel_xBH + r2yz);

    // lambda function for filling excised region
    auto interp_f = [&](auto& ah_r, auto& extrap_r, auto bh_ori, auto BH_INNER_ADAPTED_IDX) {
      // Avoid division by "0"
      if(extrap_r == 0.) extrap_r = 1e-14;
      
      double xs = x_shifted - bh_ori;
      if(xs == 0.) xs = 1e-14;
      double theta = std::acos(z / extrap_r);
      
      // atan2 is needed here
      double phi = std::atan2(y, xs); 

      // Where the filling takes places
      export_utils::spherical_turduck(
        quants, quant_vals, interp_order, delta_r_rel, interpolation_offset, 
        ah_r, extrap_r, theta, phi, BH_INNER_ADAPTED_IDX, bh_ori
      );
    };
    if (rel_rBH <= (1. + interpolation_offset) * rBH) {
      interp_f(rBH, rel_rBH, xBH, space->ADAPTEDBH+1);
    } else {
      Point abs_coords(ndim);
      abs_coords.set(1) = x - xcom_shift;
      abs_coords.set(2) = y - ycom_shift;
      abs_coords.set(3) = z;
      
      for (size_t k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
          quant_vals[k] = quants[k].get().val_point(abs_coords);
      }
    }
    
    return quant_vals;
  }
  CFMS_BHNS_Exporter::output_ary_t CFMS_BHNS_Exporter::export_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise_imp<eos_t>(x, y, z, interpolation_offset, interp_order, delta_r_rel);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise_imp<eos_t>(x, y, z, interpolation_offset, interp_order, delta_r_rel);
    } // end adding EOS OPEs
    throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  }

  CFMS_BHNS_Exporter::grid_ary_t CFMS_BHNS_Exporter::export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz,
    double const interpolation_offset, int const interp_order, double const delta_r_rel) {
    
    grid_ary_t out;
    for(auto& v : out) {
      v.resize(npoints);
    }
    
    for (size_t i = 0; i < npoints; ++i) {
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

      out[OUTPUT_VARS::RHO][i] = out_pw[OUTPUT_VARS::RHO];
      out[OUTPUT_VARS::EPS][i] = out_pw[OUTPUT_VARS::EPS];
      out[OUTPUT_VARS::PRESS][i] = out_pw[OUTPUT_VARS::PRESS];
      out[OUTPUT_VARS::VEL1][i]  = out_pw[OUTPUT_VARS::VEL1];
      out[OUTPUT_VARS::VEL2][i]  = out_pw[OUTPUT_VARS::VEL2];
      out[OUTPUT_VARS::VEL3][i]  = out_pw[OUTPUT_VARS::VEL3];
    }
    return out;
  }
}