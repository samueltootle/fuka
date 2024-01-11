#include "Solvers/bns_xcts/bns_exporter.hpp"
namespace Kadath::FUKA_Solvers {
  CFMS_BNS_Exporter::CFMS_BNS_Exporter(CFMS_BNS_Exporter const & r) {
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
  
  CFMS_BNS_Exporter& CFMS_BNS_Exporter::operator=(const CFMS_BNS_Exporter& b) {
    if (this == &b) return *this;

    CFMS_BNS_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
  }
  
  void CFMS_BNS_Exporter::initialize_eos() {
    using namespace Kadath::FUKA_Config;
    // Initialize EOS
    h_cut    = (*bconfig).template eos<double>(EOS_PARAMS::HCUT, NODES::BCO1);
    eos_file = (*bconfig).template eos<std::string>(EOS_PARAMS::EOSFILE, NODES::BCO1);
    eos_type = (*bconfig).template eos<std::string>(EOS_PARAMS::EOSTYPE, NODES::BCO1);

    if(eos_type == "Cold_Table") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = ((*bconfig).template eos<int>(EOS_PARAMS::INTERP_PTS, NODES::BCO1) == 0) ? \
                              2000 : (*bconfig).template eos<int>(EOS_PARAMS::INTERP_PTS, NODES::BCO1);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    } else if(eos_type == "Cold_PWPoly") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    }  else {
      throw std::invalid_argument("\nInvalid EOS Type\n)");
    }// end adding EOS OPEs
  }
  
  void CFMS_BNS_Exporter::load_solution_from_file() {
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

  void CFMS_BNS_Exporter::extract_computed_grid_functions() {
    // fields depending on the coords
    CoordFields<space_t> cfields(*space);
    vec_ary_t coord_vectors {default_binary_vector_ary(*space)};

    // get origin of the system and initialize coordinate fields
    double xc1 = Kadath::bco_utils::get_center(*space,space->NS1);
    double xc2 = Kadath::bco_utils::get_center(*space,space->NS2);
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
    syst.add_cst("N"  , *lapse) ;
    syst.add_cst("bet", *shift) ;
    syst.add_cst("P"  , *conformal_factor);
    syst.add_cst("H"  , *logh);
    syst.add_cst("phi", *velpotential);

    syst.add_cst("omes1" , (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA, Kadath::FUKA_Config::NODES::BCO1));
    syst.add_cst("omes2" , (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA, Kadath::FUKA_Config::NODES::BCO2));
    
    syst.add_cst("mg", *coord_vectors[Kadath::coord_vector::GLOBAL_ROT]);
    syst.add_cst("mm", *coord_vectors[Kadath::coord_vector::BCO1_ROT]) ;
    syst.add_cst("mp", *coord_vectors[Kadath::coord_vector::BCO2_ROT]) ;

    syst.add_def("h = exp(H)");
    for(int d = space->NS1; d <= space->ADAPTED1; ++d){
      syst.add_def(d, "s^i  = omes1 * mm^i");
    }
    for(int d = space->NS2; d <= space->ADAPTED2; ++d){
      syst.add_def(d, "s^i  = omes2 * mp^i");
    }
    for(int d = 0; d < ndom; ++d) {
      if((d <= space->ADAPTED1) || (d <= space->ADAPTED2 && d >= space->NS2))
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

  void CFMS_BNS_Exporter::populate_quants() {
    if(quants.capacity() != XCTS_VARS::NUM_XCTS_VARS) {
      for (size_t i = 0; i < XCTS_VARS::NUM_XCTS_VARS; ++i)
        quants.push_back(std::cref(*conformal_factor));
    }
    quants[XCTS_VARS::XCTS_PSI] = std::cref(*conformal_factor);
    quants[XCTS_VARS::XCTS_ALPHA] = std::cref(*lapse);
    quants[XCTS_VARS::XCTS_BETAX] = std::cref((*shift)(1));
    quants[XCTS_VARS::XCTS_BETAY] = std::cref((*shift)(2));
    quants[XCTS_VARS::XCTS_BETAZ] = std::cref((*shift)(3));

    export_utils::add_tensor_refs(quants, {
      XCTS_VARS::XCTS_AXX, 
      XCTS_VARS::XCTS_AXY, 
      XCTS_VARS::XCTS_AXZ, 
      XCTS_VARS::XCTS_AYY, 
      XCTS_VARS::XCTS_AYZ, 
      XCTS_VARS::XCTS_AZZ}, *A);
    
    // Fluid related quantities
    quants[XCTS_VARS::XCTS_H] = std::cref(*logh);
    quants[XCTS_VARS::XCTS_UX] = std::cref((*fluidvel)(1));
    quants[XCTS_VARS::XCTS_UY] = std::cref((*fluidvel)(2));
    quants[XCTS_VARS::XCTS_UZ] = std::cref((*fluidvel)(3));
    export_ready = true;
  }

  CFMS_BNS_Exporter::interp_ary_t CFMS_BNS_Exporter::interpolate_pointwise(double const & x, double const & y, double const & z) {
    
    double const & xcom_shift = (*bconfig)(Kadath::FUKA_Config::BIN_PARAMS::COM);
    double const & ycom_shift = (*bconfig)(Kadath::FUKA_Config::BIN_PARAMS::COMY);
    Point abs_coords(ndim);
    abs_coords.set(1) = x - xcom_shift;
    abs_coords.set(2) = y - ycom_shift;
    abs_coords.set(3) = z;
    
    for (size_t k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals;
  }

  CFMS_BNS_Exporter::interp_ary_t CFMS_BNS_Exporter::interpolate_pointwise_subset(double const & x, double const & y, double const & z,
    std::vector<CFMS_BNS_Exporter::XCTS_VARS> slice) {
    
    double const & xcom_shift = (*bconfig)(Kadath::FUKA_Config::BIN_PARAMS::COM);
    double const & ycom_shift = (*bconfig)(Kadath::FUKA_Config::BIN_PARAMS::COMY);
    Point abs_coords(ndim);
    abs_coords.set(1) = x - xcom_shift;
    abs_coords.set(2) = y - ycom_shift;
    abs_coords.set(3) = z;
    
    for (const auto k : slice) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals;
  }

  CFMS_BNS_Exporter::output_ary_t CFMS_BNS_Exporter::export_pointwise(double const & x, double const & y, double const & z) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise_imp<eos_t>(x, y, z);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise_imp<eos_t>(x, y, z);
    } // end adding EOS OPEs
    throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  }

  CFMS_BNS_Exporter::output_ary_t CFMS_BNS_Exporter::export_pointwise_fluid_vars(double const & x, double const & y, double const & z) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise_fluid_vars_imp<eos_t>(x, y, z);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise_fluid_vars_imp<eos_t>(x, y, z);
    } // end adding EOS OPEs
    throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  }

  CFMS_BNS_Exporter::output_ary_t CFMS_BNS_Exporter::export_pointwise_spacetime_vars(double const & x, double const & y, double const & z) {
    
    // Reset to NAN
    for(auto& e : quant_vals) {
      e = NAN;
    }
    quant_vals = interpolate_pointwise_subset(x, y, z, xcts_spacetime_indicies);

    // Fill output vector by storing non-conformal quantities
    auto const psi = quant_vals[XCTS_VARS::XCTS_PSI];
    auto const psi2 = psi * psi;
    auto const psi4 = psi2 * psi2;

    out_pw[OUTPUT_VARS::ALPHA] = quant_vals[XCTS_VARS::XCTS_ALPHA];

    out_pw[OUTPUT_VARS::BETAX] = quant_vals[XCTS_VARS::XCTS_BETAX];
    out_pw[OUTPUT_VARS::BETAY] = quant_vals[XCTS_VARS::XCTS_BETAY];
    out_pw[OUTPUT_VARS::BETAZ] = quant_vals[XCTS_VARS::XCTS_BETAZ];

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

    out_pw[OUTPUT_VARS::GXX] = g[0][0];
    out_pw[OUTPUT_VARS::GXY] = g[0][1];
    out_pw[OUTPUT_VARS::GXZ] = g[0][2];
    out_pw[OUTPUT_VARS::GYY] = g[1][1];
    out_pw[OUTPUT_VARS::GYZ] = g[1][2];
    out_pw[OUTPUT_VARS::GZZ] = g[2][2];

    out_pw[OUTPUT_VARS::KXX] = quant_vals[XCTS_VARS::XCTS_AXX] * psi4;
    out_pw[OUTPUT_VARS::KXY] = quant_vals[XCTS_VARS::XCTS_AXY] * psi4;
    out_pw[OUTPUT_VARS::KXZ] = quant_vals[XCTS_VARS::XCTS_AXZ] * psi4;
    out_pw[OUTPUT_VARS::KYY] = quant_vals[XCTS_VARS::XCTS_AYY] * psi4;
    out_pw[OUTPUT_VARS::KYZ] = quant_vals[XCTS_VARS::XCTS_AYZ] * psi4;
    out_pw[OUTPUT_VARS::KZZ] = quant_vals[XCTS_VARS::XCTS_AZZ] * psi4;
    
    return out_pw;
  }

  CFMS_BNS_Exporter::grid_ary_t CFMS_BNS_Exporter::export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz) {
    
    grid_ary_t out;
    for(auto& v : out) {
      v.resize(npoints);
    }
    
    for (size_t i = 0; i < npoints; ++i) {
      export_pointwise(xx[i], yy[i], zz[i]);

      out[OUTPUT_VARS::ALPHA][i] = out_pw[OUTPUT_VARS::ALPHA];

      out[OUTPUT_VARS::BETAX][i] = out_pw[OUTPUT_VARS::BETAX];
      out[OUTPUT_VARS::BETAY][i] = out_pw[OUTPUT_VARS::BETAY];
      out[OUTPUT_VARS::BETAZ][i] = out_pw[OUTPUT_VARS::BETAZ];

      out[OUTPUT_VARS::GXX][i] = out_pw[OUTPUT_VARS::GXX];
      out[OUTPUT_VARS::GXY][i] = out_pw[OUTPUT_VARS::GXY];
      out[OUTPUT_VARS::GXZ][i] = out_pw[OUTPUT_VARS::GXZ];
      out[OUTPUT_VARS::GYY][i] = out_pw[OUTPUT_VARS::GYY];
      out[OUTPUT_VARS::GYZ][i] = out_pw[OUTPUT_VARS::GYZ];
      out[OUTPUT_VARS::GZZ][i] = out_pw[OUTPUT_VARS::GZZ];

      out[OUTPUT_VARS::KXX][i] = out_pw[OUTPUT_VARS::KXX];
      out[OUTPUT_VARS::KXY][i] = out_pw[OUTPUT_VARS::KXY];
      out[OUTPUT_VARS::KXZ][i] = out_pw[OUTPUT_VARS::KXZ];
      out[OUTPUT_VARS::KYY][i] = out_pw[OUTPUT_VARS::KYY];
      out[OUTPUT_VARS::KYZ][i] = out_pw[OUTPUT_VARS::KYZ];
      out[OUTPUT_VARS::KZZ][i] = out_pw[OUTPUT_VARS::KZZ];

      out[OUTPUT_VARS::RHO][i] = out_pw[OUTPUT_VARS::RHO];
      out[OUTPUT_VARS::EPS][i] = out_pw[OUTPUT_VARS::EPS];
      out[OUTPUT_VARS::PRESS][i] = out_pw[OUTPUT_VARS::PRESS];
      out[OUTPUT_VARS::VELX][i]  = out_pw[OUTPUT_VARS::VELX];
      out[OUTPUT_VARS::VELY][i]  = out_pw[OUTPUT_VARS::VELY];
      out[OUTPUT_VARS::VELZ][i]  = out_pw[OUTPUT_VARS::VELZ];
    }
    return out;
  }
}