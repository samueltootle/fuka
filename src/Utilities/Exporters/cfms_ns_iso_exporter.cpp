#include "Solvers/ns_isotropic/ns_isotropic_exporter.hpp"
namespace Kadath::FUKA_Solvers {
  CFMS_NS_ISO_Exporter::CFMS_NS_ISO_Exporter(CFMS_NS_ISO_Exporter const & r) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;
    h_cut= r.h_cut;
    eos_file = r.eos_file;
    eos_type = r.eos_type;

    space.reset(new space_t((*r.get_space())));
    // Copy ID fields
    lap_Aterm.reset(new Scalar(*space, *r.lap_Aterm.get()));
    Nu.reset(new Scalar(*space, *r.Nu.get()));
    logh.reset(new Scalar(*space, *r.logh.get()));
    lap_Bterm.reset(new Scalar(*space, *r.lap_Bterm.get()));

    // For uniform rotation solutions
    if(r.lap_omega_term)
      lap_omega_term.reset(new Scalar(*space, *r.lap_omega_term.get()));
    else
      lap_omega_term = nullptr;
    // For differential rotation solutions
    if(r.omega)
      omega.reset(new Scalar(*space, *r.omega.get()));
    else
      omega = nullptr;

    // Copy Computed Terms    
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    metric_A.reset(new Scalar(*space, *r.metric_A.get()));
    metric_B.reset(new Scalar(*space, *r.metric_B.get()));
    metric_omega.reset(new Scalar(*space, *r.metric_omega.get()));
    domega_dr.reset(new Scalar(*space, *r.domega_dr.get()));
    domega_dt.reset(new Scalar(*space, *r.domega_dt.get()));
    fluidvel.reset(new Scalar(*space, *r.fluidvel.get()));

    bconfig.reset(new config_t(*r.bconfig));
    
    export_ready = false;

    gxx.reset(new Scalar(*space, r.gxx.get()));
    gxy.reset(new Scalar(*space, r.gxy.get()));
    gxz.reset(new Scalar(*space, r.gxz.get()));
    gyy.reset(new Scalar(*space, r.gyy.get()));
    gyz.reset(new Scalar(*space, r.gyz.get()));
    gzz.reset(new Scalar(*space, r.gzz.get()));
    Kxx.reset(new Scalar(*space, r.Kxx.get()));
    Kxy.reset(new Scalar(*space, r.Kxy.get()));
    Kxz.reset(new Scalar(*space, r.Kxz.get()));
    Kyy.reset(new Scalar(*space, r.Kyy.get()));
    Kyz.reset(new Scalar(*space, r.Kyz.get()));
    Kzz.reset(new Scalar(*space, r.Kzz.get()));
    velUx.reset(new Scalar(*space, r.velUx.get()));
    velUy.reset(new Scalar(*space, r.velUy.get()));
    velUz.reset(new Scalar(*space, r.velUz.get()));

    populate_quants();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
    // std::cout << "copy\n";
  }

  CFMS_NS_ISO_Exporter& CFMS_NS_ISO_Exporter::operator=(const CFMS_NS_ISO_Exporter& b) {
    if (this == &b) return *this;

    CFMS_NS_ISO_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
  }

  void CFMS_NS_ISO_Exporter::initialize_eos() {
    // Initialize EOS
    h_cut    = (*bconfig).template eos<double>(Kadath::FUKA_Config::EOS_PARAMS::HCUT);
    eos_file = (*bconfig).template eos<std::string>(Kadath::FUKA_Config::EOS_PARAMS::EOSFILE);
    eos_type = (*bconfig).template eos<std::string>(Kadath::FUKA_Config::EOS_PARAMS::EOSTYPE);

    if(eos_type == "Cold_Table") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = ((*bconfig).template eos<int>(Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS) == 0) ? \
                              2000 : (*bconfig).template eos<int>(Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
    } else if(eos_type == "Cold_PWPoly") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    } else {
      throw std::invalid_argument("\nInvalid EOS Type\n)");
    }// end adding EOS OPEs
  }

  void CFMS_NS_ISO_Exporter::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;

    space.reset(new space_t{ff1});
    lap_Aterm.reset( new Scalar(*space.get(), ff1)) ;
    Nu.reset( new Scalar(*space.get(), ff1)) ;
    logh.reset( new Scalar(*space.get(), ff1)) ;
    lap_Bterm.reset( new Scalar(*space.get(), ff1)) ;
    if(bconfig->set_field(Kadath::FUKA_Config::BCO_FIELDS::LAP_WTERM)) {
      lap_omega_term.reset(new Scalar(*space.get(), ff1));
     
      if(bconfig->field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)) {
        omega.reset(new Scalar(*space.get(), ff1));
        #ifdef DEBUG
        std::cout << "**** Reading Differentially rotating solution ****\n";
        #endif
      } else{
        #ifdef DEBUG
        std::cout << "**** Reading Uniformly rotating solution ****\n";
        #endif
      }
    } else {
      #ifdef DEBUG
      std::cout << "**** Reading non-rotating, spherical solution ****\n";
      #endif
    }
    
    fclose(ff1);
    ndom = space->get_nbr_domains();
  }

  void CFMS_NS_ISO_Exporter::extract_computed_grid_functions() {
    System_of_eqs syst(*space);

    Param p;
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      set_eos_ope<eos_t>(syst, p);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      set_eos_ope<eos_t>(syst, p);
    } // end adding EOS OPEs

    // Fields - must be initialized before common setup
    syst.add_cst("nu", *Nu);
    syst.add_cst("lapAterm", *lap_Aterm);
    syst.add_cst("lapBterm", *lap_Bterm);

    if(bconfig->field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)) {
      #ifdef DEBUG
      std::cout << "**** Importing differential rotation profile ****\n";
      #endif
      syst.add_cst("Omega", *omega);
    } else {
      #ifdef DEBUG
      std::cout << "**** Importing uniform rotation profile ****\n";
      #endif
      syst.add_cst("Omega", (*bconfig)(BCO_PARAMS::OMEGA));
    }
    syst.add_def("N = exp(nu)");
    syst.add_def("A = exp(lapAterm - nu)");
    syst.add_def("B = (divrsint(lapBterm) + 1) / N");

    if(bconfig->set_field(Kadath::FUKA_Config::BCO_FIELDS::LAP_WTERM)) {
      syst.add_cst("wrsint"  , *lap_omega_term);
      syst.add_def("w = divrsint(wrsint)");
      syst.add_def("drw = dr(w)");
      syst.add_def("dtw = dt(w)");

      syst.add_def("UphiU = 1 / N * (Omega - w)");

      fluidvel.reset(new Scalar(syst.give_val_def("UphiU")));
      fluidvel->coef();

      metric_omega.reset(new Scalar(syst.give_val_def("w")));
      metric_omega->coef();
      domega_dr.reset(new Scalar(syst.give_val_def("drw")));
      domega_dr->coef();
      domega_dt.reset(new Scalar(syst.give_val_def("dtw")));
      domega_dt->coef();
    } else {
      #ifdef DEBUG
      std::cout << "**** Importing non-rotation profile ****\n";
      #endif
      fluidvel.reset(new Scalar(*space));
      fluidvel->annule_hard();
      fluidvel->std_base();

      metric_omega.reset(new Scalar(*space));
      metric_omega->annule_hard();      
      metric_omega->std_base();
      domega_dr.reset(new Scalar(*space));
      domega_dr->annule_hard();
      domega_dr->std_base();
      domega_dt.reset(new Scalar(*space));
      domega_dt->annule_hard();
      domega_dt->std_base();
    }
    metric_A.reset(new Scalar(syst.give_val_def("A")));
    metric_A->coef();
    metric_B.reset(new Scalar(syst.give_val_def("B")));
    metric_B->coef();
    lapse.reset(new Scalar(syst.give_val_def("N")));
    lapse->coef();
  }

  void CFMS_NS_ISO_Exporter::populate_quants() {
    if(quants.capacity() != ISO_VARS::NUM_ISO_VARS) {
      for (size_t i = 0; i < ISO_VARS::NUM_ISO_VARS; ++i)
        quants.push_back(std::cref(*Nu));
    }
    quants[ISO_VARS::ISO_METRIC_A] = std::cref(*metric_A);
    quants[ISO_VARS::ISO_METRIC_B] = std::cref(*metric_B);
    quants[ISO_VARS::ISO_ALPHA] = std::cref(*lapse);
    quants[ISO_VARS::ISO_METRIC_OMEGA] = std::cref(*metric_omega);
    quants[ISO_VARS::ISO_DOMEGA_DR] = std::cref(*domega_dr);
    quants[ISO_VARS::ISO_DOMEGA_DTHETA] = std::cref(*domega_dt);

    
    // Fluid related quantities
    quants[ISO_VARS::ISO_H] = std::cref(*logh);
    quants[ISO_VARS::ISO_U] = std::cref(*fluidvel);

    if(output_quants.capacity() != CART_VARS::NUM_CART_VARS) {
      for (size_t i = 0; i < CART_VARS::NUM_CART_VARS; ++i)
        output_quants.push_back(std::cref(*lapse));
    }

    output_quants[CART_VARS::CART_ALPHA] = std::cref(*lapse);
    output_quants[CART_VARS::CART_H] = std::cref(*logh);
    output_quants[CART_VARS::CART_G11] = std::cref(*gxx);
    output_quants[CART_VARS::CART_G12] = std::cref(*gxy);
    output_quants[CART_VARS::CART_G13] = std::cref(*gxz);
    output_quants[CART_VARS::CART_G22] = std::cref(*gyy);
    output_quants[CART_VARS::CART_G23] = std::cref(*gyz);
    output_quants[CART_VARS::CART_G33] = std::cref(*gzz);
    output_quants[CART_VARS::CART_K11] = std::cref(*Kxx);
    output_quants[CART_VARS::CART_K12] = std::cref(*Kxy);
    output_quants[CART_VARS::CART_K13] = std::cref(*Kxz);
    output_quants[CART_VARS::CART_K22] = std::cref(*Kyy);
    output_quants[CART_VARS::CART_K23] = std::cref(*Kyz);
    output_quants[CART_VARS::CART_K33] = std::cref(*Kzz);
    output_quants[CART_VARS::CART_VEL1] = std::cref(*velUx);
    output_quants[CART_VARS::CART_VEL2] = std::cref(*velUy);
    output_quants[CART_VARS::CART_VEL3] = std::cref(*velUz);    

    export_ready = true;
  }

  CFMS_NS_ISO_Exporter::interp_cart_ary_t CFMS_NS_ISO_Exporter::interpolate_pointwise(double const & x, double const & y, double const & z) {
    
    double r2_xy = x * x + y * y;
    double r_xy = std::sqrt(r2_xy);
    Point abs_coords(ndim);
    abs_coords.set(1) = r_xy;
    abs_coords.set(2) = z;
    
    for (size_t k = 0; k < CART_VARS::NUM_CART_VARS; ++k) {
        quant_vals[k] = output_quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals;
  }

  CFMS_NS_ISO_Exporter::interp_ary_t CFMS_NS_ISO_Exporter::interpolate_pointwise__spherical(double const & x, double const & y, double const & z) {
    
    double r2_xy = x * x + y * y;
    double r_xy = std::sqrt(r2_xy);
    Point abs_coords(ndim);
    abs_coords.set(1) = r_xy;
    abs_coords.set(2) = z;

    CFMS_NS_ISO_Exporter::interp_ary_t quant_vals__sph;
    
    for (size_t k = 0; k < ISO_VARS::NUM_ISO_VARS; ++k) {
        quant_vals__sph[k] = quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals__sph;
  }

  CFMS_NS_ISO_Exporter::interp_cart_ary_t CFMS_NS_ISO_Exporter::interpolate_pointwise_subset(double const & x, double const & y, double const & z,
    std::vector<CFMS_NS_ISO_Exporter::CART_VARS> slice) {
    
    double r2_xy = x * x + y * y;
    double r_xy = std::sqrt(r2_xy);
    Point abs_coords(ndim);
    abs_coords.set(1) = r_xy;
    abs_coords.set(2) = z;
    
    for (const auto k : slice) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals;
  }

  CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise(double const & x, double const & y, double const & z) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise_imp<eos_t>(x, y, z);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise_imp<eos_t>(x, y, z);
    } // end adding EOS OPEs
    std::string msg("\nExport: Invalid EOS Type: " + eos_type + "\n");
    throw std::invalid_argument(msg.c_str());
  }

  CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise_fluid_vars(double const & x, double const & y, double const & z) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise_fluid_vars_imp<eos_t>(x, y, z);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise_fluid_vars_imp<eos_t>(x, y, z);
    } // end adding EOS OPEs
    throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  }

  CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise_spacetime_vars(double const & x, double const & y, double const & z) {
    // Reset to NAN
    
    quant_vals.fill(NAN);
    quant_vals = interpolate_pointwise_subset(x, y, z, spacetime_indicies__cart);

    out_pw[OUTPUT_VARS::ALPHA] = quant_vals[CART_VARS::CART_ALPHA];

    out_pw[OUTPUT_VARS::BETA1] = 0.0;
    out_pw[OUTPUT_VARS::BETA2] = 0.0;
    out_pw[OUTPUT_VARS::BETA3] = 0.0;

    out_pw[OUTPUT_VARS::G11] = quant_vals[CART_VARS::CART_G11];
    out_pw[OUTPUT_VARS::G12] = quant_vals[CART_VARS::CART_G12];
    out_pw[OUTPUT_VARS::G13] = quant_vals[CART_VARS::CART_G13];
    out_pw[OUTPUT_VARS::G22] = quant_vals[CART_VARS::CART_G22];
    out_pw[OUTPUT_VARS::G23] = quant_vals[CART_VARS::CART_G23];
    out_pw[OUTPUT_VARS::G33] = quant_vals[CART_VARS::CART_G33];
    out_pw[OUTPUT_VARS::K11] = quant_vals[CART_VARS::CART_K11];
    out_pw[OUTPUT_VARS::K12] = quant_vals[CART_VARS::CART_K12];
    out_pw[OUTPUT_VARS::K13] = quant_vals[CART_VARS::CART_K13];
    out_pw[OUTPUT_VARS::K22] = quant_vals[CART_VARS::CART_K22];
    out_pw[OUTPUT_VARS::K23] = quant_vals[CART_VARS::CART_K23];
    out_pw[OUTPUT_VARS::K33] = quant_vals[CART_VARS::CART_K33];

    output_base = OUTPUT_BASIS::CARTESIAN;
    return out_pw;    
  }

  // CFMS_NS_ISO_Exporter::grid_ary_t CFMS_NS_ISO_Exporter::export_coordinate_array(
  //   int const npoints, double const * xx, double const * yy, double const * zz) {
    
  //   grid_ary_t out;
  //   for(auto& v : out) {
  //     v.resize(npoints);
  //   }
    
  //   for (size_t i = 0; i < npoints; ++i) {
  //     export_pointwise(xx[i], yy[i], zz[i]);

  //     out[OUTPUT_VARS::ALPHA][i] = out_pw[OUTPUT_VARS::ALPHA];

  //     out[OUTPUT_VARS::BETA1][i] = out_pw[OUTPUT_VARS::BETA1];
  //     out[OUTPUT_VARS::BETA2][i] = out_pw[OUTPUT_VARS::BETA2];
  //     out[OUTPUT_VARS::BETA3][i] = out_pw[OUTPUT_VARS::BETA3];

  //     out[OUTPUT_VARS::G11][i] = out_pw[OUTPUT_VARS::G11];
  //     out[OUTPUT_VARS::G12][i] = out_pw[OUTPUT_VARS::G12];
  //     out[OUTPUT_VARS::G13][i] = out_pw[OUTPUT_VARS::G13];
  //     out[OUTPUT_VARS::G22][i] = out_pw[OUTPUT_VARS::G22];
  //     out[OUTPUT_VARS::G23][i] = out_pw[OUTPUT_VARS::G23];
  //     out[OUTPUT_VARS::G33][i] = out_pw[OUTPUT_VARS::G33];

  //     out[OUTPUT_VARS::K11][i] = out_pw[OUTPUT_VARS::K11];
  //     out[OUTPUT_VARS::K12][i] = out_pw[OUTPUT_VARS::K12];
  //     out[OUTPUT_VARS::K13][i] = out_pw[OUTPUT_VARS::K13];
  //     out[OUTPUT_VARS::K22][i] = out_pw[OUTPUT_VARS::K22];
  //     out[OUTPUT_VARS::K23][i] = out_pw[OUTPUT_VARS::K23];
  //     out[OUTPUT_VARS::K33][i] = out_pw[OUTPUT_VARS::K33];

  //     out[OUTPUT_VARS::RHO][i] = out_pw[OUTPUT_VARS::RHO];
  //     out[OUTPUT_VARS::EPS][i] = out_pw[OUTPUT_VARS::EPS];
  //     out[OUTPUT_VARS::PRESS][i] = out_pw[OUTPUT_VARS::PRESS];
  //     out[OUTPUT_VARS::VEL1][i]  = out_pw[OUTPUT_VARS::VEL1];
  //     out[OUTPUT_VARS::VEL2][i]  = out_pw[OUTPUT_VARS::VEL2];
  //     out[OUTPUT_VARS::VEL3][i]  = out_pw[OUTPUT_VARS::VEL3];
  //   }
  //   return out;
  // }

  CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise__spherical(double const & x, double const & y, double const & z) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise__spherical_imp<eos_t>(x, y, z);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise__spherical_imp<eos_t>(x, y, z);
    } // end adding EOS OPEs
    std::string msg("\nExport: Invalid EOS Type: " + eos_type + "\n");
    throw std::invalid_argument(msg.c_str());
  }

  // FIXME replace codegen code with something readable - preferably
  // generic
  void CFMS_NS_ISO_Exporter::basis_transform__Spherical_to_Cart() {
    for(int d = 0; d < ndom; ++d) {
      auto dom  = space->get_domain(d);
      auto XX = dom->get_cart(1);
      auto ZZ = dom->get_cart(2);
      
      auto npts = space->get_domain(d)->get_nbr_points();
      Index pos(npts);
      
      auto ADM_Spherical_to_Cart = [&] (auto const x, auto const y, auto const z) {
        using REAL = double;
        const REAL xCart[3] = {x, y, z};
        // Perform the basis transform on ADM vectors/tensors from Spherical to Cartesian:

        // Set destination xx[3] based on desired xCart[3]
        REAL xx0, xx1, xx2;
        /*
        *  Original SymPy expressions:
        *  "[xx0 = sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2)]"
        *  "[xx1 = acos(xCart[2]/sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2))]"
        *  "[xx2 = atan2(xCart[1], xCart[0])]"
        */
        {
          const REAL tmp0 = sqrt(((xCart[0]) * (xCart[0])) + ((xCart[1]) * (xCart[1])) + ((xCart[2]) * (xCart[2])));
          xx0 = tmp0;
          if(std::fabs(xx0) < 1e-12) {
            xx0 = 1e-12;
            xx1 = acos(xCart[2] / xx0); //theta (angle from Z to xy plane)
            const double xfixed = std::copysign(1e-12, xCart[0]);
            xx2 = atan2(xCart[1], xfixed); // phi (angle from x to y axis)
          } else if(std::fabs(x) < 1e-12 && std::fabs(y) < 1e-12) {
            xx1 = 1e-10; //theta (angle from Z to xy plane)
            
            const double xfixed = std::copysign(1e-12, xCart[0]);
            xx2 = atan2(xCart[1], xfixed); // phi (angle from x to y axis)
          } else {
            xx1 = acos(xCart[2] / xx0); //theta (angle from Z to xy plane)
            xx2 = atan2(xCart[1], xCart[0]); // phi (angle from x to y axis)
          }
          // cout << "r: " << xx0 << ", t: " << xx1 << ", phi: " << xx2 << endl;
        }
        // Unpack initial_data for ADM vectors/tensors
        const REAL N = (*lapse)(d)(pos);
        const REAL U_factor = (*fluidvel)(d)(pos);
        const REAL domega_dr_ = (*domega_dr)(d)(pos);
        const REAL domega_dt_ = (*domega_dt)(d)(pos);

        const REAL FluidVelU0 = 0.0;                                   // r
        const REAL FluidVelU1 = 0.0;                                   // theta
        const REAL FluidVelU2 = U_factor;                              // phi
        const REAL betaSphericalU0 = 0.0;                              // r
        const REAL betaSphericalU1 = 0.0;                              // theta
        const REAL betaSphericalU2 = -(*metric_omega)(d)(pos); // phi

        const REAL A = (*metric_A)(d)(pos);
        const REAL B = (*metric_B)(d)(pos);

        const REAL gammaSphericalDD01 = 0.0;
        const REAL gammaSphericalDD02 = 0.0;
        const REAL gammaSphericalDD12 = 0.0;
        const REAL gammaSphericalDD00 = A * A;
        const REAL gammaSphericalDD11 = A * A * xx0 * xx0;
        const REAL gammaSphericalDD22 = B * B * xx0 * xx0 * sin(xx1) * sin(xx1);
        // cout << "g11: " << gammaSphericalDD00 << ", g22: " << gammaSphericalDD11 << ", g33: " << gammaSphericalDD22 << endl;

        const REAL KSphericalDD00 = 0.0;
        const REAL KSphericalDD01 = 0.0;
        const REAL KSphericalDD11 = 0.0;
        const REAL KSphericalDD22 = 0.0;
        const REAL KSphericalDD02 = -gammaSphericalDD22 / 2.0 / N * domega_dr_;
        const REAL KSphericalDD12 = -gammaSphericalDD22 / 2.0 / N * domega_dt_;
        const REAL tmp0 = cos(xx2);
        const REAL tmp1 = sin(xx1);
        const REAL tmp4 = cos(xx1);
        const REAL tmp6 = sin(xx2);
        const REAL tmp12 = ((xx0) * (xx0));
        const REAL tmp3 = tmp0 * xx0;
        const REAL tmp7 = tmp1 * xx0;
        const REAL tmp9 = tmp6 * xx0;
        const REAL tmp10 = ((tmp0) * (tmp0));
        const REAL tmp11 = ((tmp6) * (tmp6));
        const REAL tmp13 = ((tmp1) * (tmp1) * (tmp1));
        const REAL tmp15 = ((tmp4) * (tmp4));
        const REAL tmp20 = ((tmp1) * (tmp1) * (tmp1) * (tmp1)) * ((xx0) * (xx0) * (xx0) * (xx0));
        const REAL tmp25 = ((tmp1) * (tmp1));
        const REAL tmp17 = tmp1 * tmp12 * tmp15;
        const REAL tmp21 = gammaSphericalDD00 * tmp20;
        const REAL tmp23 = tmp13 * tmp4 * ((xx0) * (xx0) * (xx0));
        const REAL tmp26 = tmp12 * tmp25;
        const REAL tmp29 = -tmp15 * tmp9 - tmp25 * tmp9;
        const REAL tmp34 = tmp4 * tmp7;
        const REAL tmp41 = tmp15 * tmp3 + tmp25 * tmp3;
        const REAL tmp52 = tmp1 * tmp12 * tmp4;
        const REAL tmp68 = KSphericalDD00 * tmp20;
        const REAL tmp18 = (1.0 / ((tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17) *
                                  (tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17)));
        const REAL tmp24 = 2 * gammaSphericalDD01 * tmp23;
        const REAL tmp35 = gammaSphericalDD12 * tmp34;
        const REAL tmp37 = gammaSphericalDD02 * tmp26;
        const REAL tmp47 = -tmp10 * tmp25 * xx0 - tmp11 * tmp25 * xx0;
        const REAL tmp53 = tmp10 * tmp52 + tmp11 * tmp52;
        const REAL tmp70 = 2 * KSphericalDD01 * tmp23;
        const REAL tmp73 = KSphericalDD12 * tmp34;
        const REAL tmp75 = KSphericalDD02 * tmp26;
        const REAL tmp19 = tmp10 * tmp18;
        const REAL tmp28 = gammaSphericalDD11 * tmp15 * tmp26;
        const REAL tmp31 = gammaSphericalDD22 * tmp18;
        const REAL tmp32 = tmp18 * tmp29;
        const REAL tmp39 = tmp0 * tmp18;
        const REAL tmp42 = tmp18 * tmp41;
        const REAL tmp54 = tmp18 * tmp53;
        const REAL tmp61 = tmp11 * tmp18;
        const REAL tmp64 = tmp18 * tmp6;
        const REAL tmp65 = tmp18 * ((tmp47) * (tmp47));
        const REAL tmp66 = tmp18 * ((tmp53) * (tmp53));
        const REAL tmp71 = KSphericalDD11 * tmp15 * tmp26;
        const REAL tmp72 = KSphericalDD22 * tmp18;
        const REAL tmp33 = tmp0 * tmp32;
        const REAL tmp40 = tmp39 * tmp6;
        const REAL tmp43 = tmp0 * tmp42;
        const REAL tmp44 = tmp32 * tmp6;
        const REAL tmp49 = gammaSphericalDD11 * tmp34 * tmp47;
        const REAL tmp51 = gammaSphericalDD01 * tmp26 * tmp47;
        const REAL tmp63 = tmp42 * tmp6;
        const REAL tmp77 = KSphericalDD11 * tmp34 * tmp47;
        const REAL tmp78 = KSphericalDD01 * tmp26 * tmp47;
        const REAL tmp56 = gammaSphericalDD01 * tmp34 * tmp54;
        const REAL tmp58 = gammaSphericalDD00 * tmp26 * tmp54;
        const REAL tmp79 = KSphericalDD01 * tmp34 * tmp54;
        const REAL tmp80 = KSphericalDD00 * tmp26 * tmp54;
        // output_quants[OUTPUT_VARS::BETA1]->set_domain(d).set(pos) = betaSphericalU0 * tmp0 * tmp1 + betaSphericalU1 * tmp3 * tmp4 - betaSphericalU2 * tmp6 * tmp7;
        // output_quants[OUTPUT_VARS::BETA2]->set_domain(d).set(pos) = betaSphericalU0 * tmp1 * tmp6 + betaSphericalU1 * tmp4 * tmp9 + betaSphericalU2 * tmp0 * tmp7;
        // output_quants[OUTPUT_VARS::BETA3]->set_domain(d).set(pos) = betaSphericalU0 * tmp4 - betaSphericalU1 * tmp7;
        (*gxx).set_domain(d).set(pos) = tmp19 * tmp21 + tmp19 * tmp24 + tmp19 * tmp28 + ((tmp29) * (tmp29)) * tmp31 + 2 * tmp33 * tmp35 + 2 * tmp33 * tmp37;
        (*gxy).set_domain(d).set(pos) =
            tmp21 * tmp40 + tmp24 * tmp40 + tmp28 * tmp40 + tmp29 * tmp31 * tmp41 + tmp35 * tmp43 + tmp35 * tmp44 + tmp37 * tmp43 + tmp37 * tmp44;
        (*gxz).set_domain(d).set(pos) =
            gammaSphericalDD02 * tmp32 * tmp53 + gammaSphericalDD12 * tmp32 * tmp47 + tmp0 * tmp56 + tmp0 * tmp58 + tmp39 * tmp49 + tmp39 * tmp51;
        (*gyy).set_domain(d).set(pos) = tmp21 * tmp61 + tmp24 * tmp61 + tmp28 * tmp61 + tmp31 * ((tmp41) * (tmp41)) + 2 * tmp35 * tmp63 + 2 * tmp37 * tmp63;
        (*gyz).set_domain(d).set(pos) =
            gammaSphericalDD02 * tmp42 * tmp53 + gammaSphericalDD12 * tmp42 * tmp47 + tmp49 * tmp64 + tmp51 * tmp64 + tmp56 * tmp6 + tmp58 * tmp6;
        (*gzz).set_domain(d).set(pos) = gammaSphericalDD00 * tmp66 + 2 * gammaSphericalDD01 * tmp47 * tmp54 + gammaSphericalDD11 * tmp65;
        (*Kxx).set_domain(d).set(pos) = tmp19 * tmp68 + tmp19 * tmp70 + tmp19 * tmp71 + ((tmp29) * (tmp29)) * tmp72 + 2 * tmp33 * tmp73 + 2 * tmp33 * tmp75;
        (*Kxy).set_domain(d).set(pos) =
            tmp29 * tmp41 * tmp72 + tmp40 * tmp68 + tmp40 * tmp70 + tmp40 * tmp71 + tmp43 * tmp73 + tmp43 * tmp75 + tmp44 * tmp73 + tmp44 * tmp75;
        (*Kxz).set_domain(d).set(pos) =
            KSphericalDD02 * tmp32 * tmp53 + KSphericalDD12 * tmp32 * tmp47 + tmp0 * tmp79 + tmp0 * tmp80 + tmp39 * tmp77 + tmp39 * tmp78;
        (*Kyy).set_domain(d).set(pos) = ((tmp41) * (tmp41)) * tmp72 + tmp61 * tmp68 + tmp61 * tmp70 + tmp61 * tmp71 + 2 * tmp63 * tmp73 + 2 * tmp63 * tmp75;
        (*Kyz).set_domain(d).set(pos) =
            KSphericalDD02 * tmp42 * tmp53 + KSphericalDD12 * tmp42 * tmp47 + tmp6 * tmp79 + tmp6 * tmp80 + tmp64 * tmp77 + tmp64 * tmp78;
        (*Kzz).set_domain(d).set(pos) = KSphericalDD00 * tmp66 + 2 * KSphericalDD01 * tmp47 * tmp54 + KSphericalDD11 * tmp65;
        (*velUx).set_domain(d).set(pos) = FluidVelU0 * tmp0 * tmp1 + FluidVelU1 * tmp3 * tmp4 - FluidVelU2 * tmp6 * tmp7;
        (*velUy).set_domain(d).set(pos) = FluidVelU0 * tmp1 * tmp6 + FluidVelU1 * tmp4 * tmp9 + FluidVelU2 * tmp0 * tmp7;
        (*velUz).set_domain(d).set(pos) = FluidVelU0 * tmp4 - FluidVelU1 * tmp7;
      };
      do {
        auto const x = XX(pos);
        auto const z = ZZ(pos);
        auto const y = 0.;
        ADM_Spherical_to_Cart(x, y, z);

      }while(pos.inc());      
    }
    // Ensure velocity is zero outside the star
    for(int d = space->ADAPTED_INNER; d < ndom; ++d) {
      velUx->set_domain(d).annule_hard();
      velUy->set_domain(d).annule_hard();
      velUz->set_domain(d).annule_hard();
    }
    {
      // enforce outer boundary conditions (flat space)
      int d = ndom - 1;
      auto npts = space->get_domain(d)->get_nbr_points();
      Index pos(npts);
      pos.set(0) = npts(0) - 1;
      for(int k = 0; k < npts(1); ++k) {
        pos.set(1) = k;
        gxx->set_domain(d).set(pos) = 1.0;
        gxy->set_domain(d).set(pos) = 0.0;
        gxz->set_domain(d).set(pos) = 0.0;
        gyy->set_domain(d).set(pos) = 1.0;
        gyz->set_domain(d).set(pos) = 0.0;
        gzz->set_domain(d).set(pos) = 1.0;
        Kxx->set_domain(d).set(pos) = 0.0;
        Kxy->set_domain(d).set(pos) = 0.0;
        Kxz->set_domain(d).set(pos) = 0.0;
        Kyy->set_domain(d).set(pos) = 0.0;
        Kyz->set_domain(d).set(pos) = 0.0;
        Kzz->set_domain(d).set(pos) = 0.0;
      }
    }
    // Compute spectral coefficients for all fields
    gxx->coef();
    gxy->coef();
    gxz->coef();
    gyy->coef();
    gyz->coef();
    gzz->coef();
    Kxx->coef();
    Kxy->coef();
    Kxz->coef();
    Kyy->coef();
    Kyz->coef();
    Kzz->coef();
    velUx->coef();
    velUy->coef();
    velUz->coef();
  }

  void CFMS_NS_ISO_Exporter::initialize_Cartesian_fields() {
    // Initialize pointers
    gxx.reset(new Scalar(*space));
    gxy.reset(new Scalar(*space));
    gxz.reset(new Scalar(*space));
    gyy.reset(new Scalar(*space));
    gyz.reset(new Scalar(*space));
    gzz.reset(new Scalar(*space));
    Kxx.reset(new Scalar(*space));
    Kxy.reset(new Scalar(*space));
    Kxz.reset(new Scalar(*space));
    Kyy.reset(new Scalar(*space));
    Kyz.reset(new Scalar(*space));
    Kzz.reset(new Scalar(*space));
    velUx.reset(new Scalar(*space));
    velUy.reset(new Scalar(*space));
    velUz.reset(new Scalar(*space));
    
    // Set fields to zero
    gxx->annule_hard();
    gxy->annule_hard();
    gxz->annule_hard();
    gyy->annule_hard();
    gyz->annule_hard();
    gzz->annule_hard();
    Kxx->annule_hard();
    Kxy->annule_hard();
    Kxz->annule_hard();
    Kyy->annule_hard();
    Kyz->annule_hard();
    Kzz->annule_hard();
    velUx->annule_hard();
    velUy->annule_hard();
    velUz->annule_hard();

    // Set to standard spectral basis
    gxx->std_base();
    gxy->std_base();
    gxz->std_base();
    gyy->std_base();
    gyz->std_base();
    gzz->std_base();
    Kxx->std_base();
    Kxy->std_base();
    Kxz->std_base();
    Kyy->std_base();
    Kyz->std_base();
    Kzz->std_base();
    velUx->std_base();
    velUy->std_base();
    velUz->std_base();
  }
}
