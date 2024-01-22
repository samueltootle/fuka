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
    lap_Bterm.reset(new Scalar(*space, *r.lap_Bterm.get()));
    Nu.reset(new Scalar(*space, *r.Nu.get()));
    lap_omega_term.reset(new Scalar(*space, *r.lap_omega_term.get()));
    logh.reset(new Scalar(*space, *r.logh.get()));
    omega.reset(new Scalar(*space, *r.omega.get()));

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
        std::cout << "**** Reading Differentially rotating solution ****\n"
      } else{
        std::cout << "**** Reading Uniformly rotating solution ****\n"
      }
    } else {
      std::cout << "**** Reading non-rotating, spherical solution ****\n"
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
    syst.add_cst("wrsint"  , *lap_omega_term);
    if(bconfig->field(Kadath::FUKA_Config::BCO_FIELDS::DIFF_OMEGA)) {
      std::cout << "**** Importing differential rotation profile ****\n"
      syst.add_cst("Omega", *omega);
    } else {
      std::cout << "**** Importing uniform rotation profile ****\n"
      syst.add_cst("Omega", (*bconfig)(BCO_PARAMS::OMEGA));
    }
    syst.add_def("N = exp(nu)");
    syst.add_def("A = exp(lapAterm - nu)");
    syst.add_def("B = (divrsint(lapBterm) + 1) / N");
    syst.add_def("w = divrsint(wrsint)");
    syst.add_def("drw = dr(w)");
    syst.add_def("dtw = dt(w)");

    if(bconfig->set_field(Kadath::FUKA_Config::BCO_FIELDS::LAP_WTERM)) {
      syst.add_def("UphiU = 1 / N * (Omega - w)");

      fluidvel.reset(new Scalar(syst.give_val_def("UphiU")));
      fluidvel->coef();
    } else {
      fluidvel.reset(new Scalar(*space));
      fluidvel->annule_hard();
    }
    metric_A.reset(new Scalar(syst.give_val_def("A")));
    metric_A->coef();
    metric_B.reset(new Scalar(syst.give_val_def("B")));
    metric_B->coef();
    lapse.reset(new Scalar(syst.give_val_def("N")));
    lapse->coef();
    metric_omega.reset(new Scalar(syst.give_val_def("w")));
    metric_omega->coef();
    domega_dr.reset(new Scalar(syst.give_val_def("drw")));
    domega_dr->coef();
    domega_dt.reset(new Scalar(syst.give_val_def("dtw")));
    domega_dt->coef();
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
    export_ready = true;
  }

  CFMS_NS_ISO_Exporter::interp_ary_t CFMS_NS_ISO_Exporter::interpolate_pointwise(double const & x, double const & y, double const & z) {
    
    double r2_xy = x * x + y * y;
    double r_xy = std::sqrt(r2_xy);
    Point abs_coords(ndim);
    abs_coords.set(1) = r_xy;
    abs_coords.set(2) = z;
    
    for (size_t k = 0; k < ISO_VARS::NUM_ISO_VARS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals;
  }

  CFMS_NS_ISO_Exporter::interp_ary_t CFMS_NS_ISO_Exporter::interpolate_pointwise_subset(double const & x, double const & y, double const & z,
    std::vector<CFMS_NS_ISO_Exporter::ISO_VARS> slice) {
    
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
    throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  }

  // CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise_fluid_vars(double const & x, double const & y, double const & z) {
  //   if(eos_type == "Cold_Table") {
  //     using eos_t = Kadath::Margherita::Cold_Table;
  //     return export_pointwise_fluid_vars_imp<eos_t>(x, y, z);
  //   }else if(eos_type == "Cold_PWPoly") {
  //     using eos_t = Kadath::Margherita::Cold_PWPoly;
  //     return export_pointwise_fluid_vars_imp<eos_t>(x, y, z);
  //   } // end adding EOS OPEs
  //   throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  // }

  // CFMS_NS_ISO_Exporter::output_ary_t CFMS_NS_ISO_Exporter::export_pointwise_spacetime_vars(double const & x, double const & y, double const & z) {
    
  //   // Reset to NAN
  //   for(auto& e : quant_vals) {
  //     e = NAN;
  //   }
  //   quant_vals = interpolate_pointwise_subset(x, y, z, iso_spacetime_indicies);

  //   /*
  //   * Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis
  //   * Code generated from NrPyv2
  //   */
  //   auto ADM_Spherical_to_Cart =[&]() {
  //     using REAL = double; 
  //     const double xCart[3] = {x, y, z};
  //     // Perform the basis transform on ADM vectors/tensors from Spherical to Cartesian:

  //     // Set destination xx[3] based on desired xCart[3]
  //     double xx0, xx1, xx2;
  //     /*
  //     *  Original SymPy expressions:
  //     *  "[xx0 = sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2)]"
  //     *  "[xx1 = acos(xCart[2]/sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2))]"
  //     *  "[xx2 = atan2(xCart[1], xCart[0])]"
  //     */
  //     {
  //       const double tmp0 = sqrt(((xCart[0]) * (xCart[0])) + ((xCart[1]) * (xCart[1])) + ((xCart[2]) * (xCart[2])));
  //       xx0 = tmp0;
  //       if(std::fabs(xx0) < 1e-12) xx0 = 1e-12;
  //       const REAL X = (std::fabs(xCart[0]) < 1e-12) ? 1e-12 : xCart[0];
  //       xx1 = acos(xCart[2] / xx0); //theta (angle from Z to xy plane)
  //       xx2 = atan2(xCart[1], X); // phi (angle from x to y axis)
  //     }
  //     // Unpack initial_data for ADM vectors/tensors
  //     const double N = quant_vals[ISO_VARS::ISO_ALPHA];
  //     const double domega_dr = quant_vals[ISO_VARS::ISO_DOMEGA_DR];
  //     const double domega_dt = quant_vals[ISO_VARS::ISO_DOMEGA_DTHETA];

  //     const double betaSphericalU0 = 0.0;                              // r
  //     const double betaSphericalU1 = 0.0;                              // theta
  //     const double betaSphericalU2 = -quant_vals[ISO_VARS::ISO_METRIC_OMEGA]; // phi

  //     const double A = quant_vals[ISO_VARS::ISO_METRIC_A];
  //     const double B = quant_vals[ISO_VARS::ISO_METRIC_B];

  //     const double gammaSphericalDD01 = 0.0;
  //     const double gammaSphericalDD02 = 0.0;
  //     const double gammaSphericalDD12 = 0.0;
  //     const double gammaSphericalDD00 = A * A;
  //     const double gammaSphericalDD11 = A * A * xx0 * xx0;
  //     const double gammaSphericalDD22 = B * B * xx0 * xx0 * sin(xx1) * sin(xx1);

  //     const double KSphericalDD00 = 0.0;
  //     const double KSphericalDD01 = 0.0;
  //     const double KSphericalDD11 = 0.0;
  //     const double KSphericalDD22 = 0.0;
  //     const double KSphericalDD02 = -gammaSphericalDD22 / 2.0 / N * domega_dr;
  //     const double KSphericalDD12 = -gammaSphericalDD22 / 2.0 / N * domega_dt;
  //     const double tmp0 = cos(xx2);
  //     const double tmp1 = sin(xx1);
  //     const double tmp4 = cos(xx1);
  //     const double tmp6 = sin(xx2);
  //     const double tmp12 = ((xx0) * (xx0));
  //     const double tmp3 = tmp0 * xx0;
  //     const double tmp7 = tmp1 * xx0;
  //     const double tmp9 = tmp6 * xx0;
  //     const double tmp10 = ((tmp0) * (tmp0));
  //     const double tmp11 = ((tmp6) * (tmp6));
  //     const double tmp13 = ((tmp1) * (tmp1) * (tmp1));
  //     const double tmp15 = ((tmp4) * (tmp4));
  //     const double tmp20 = ((tmp1) * (tmp1) * (tmp1) * (tmp1)) * ((xx0) * (xx0) * (xx0) * (xx0));
  //     const double tmp25 = ((tmp1) * (tmp1));
  //     const double tmp17 = tmp1 * tmp12 * tmp15;
  //     const double tmp21 = gammaSphericalDD00 * tmp20;
  //     const double tmp23 = tmp13 * tmp4 * ((xx0) * (xx0) * (xx0));
  //     const double tmp26 = tmp12 * tmp25;
  //     const double tmp29 = -tmp15 * tmp9 - tmp25 * tmp9;
  //     const double tmp34 = tmp4 * tmp7;
  //     const double tmp41 = tmp15 * tmp3 + tmp25 * tmp3;
  //     const double tmp52 = tmp1 * tmp12 * tmp4;
  //     const double tmp68 = KSphericalDD00 * tmp20;
  //     const double tmp18 = (1.0 / ((tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17) *
  //                               (tmp10 * tmp12 * tmp13 + tmp10 * tmp17 + tmp11 * tmp12 * tmp13 + tmp11 * tmp17)));
  //     const double tmp24 = 2 * gammaSphericalDD01 * tmp23;
  //     const double tmp35 = gammaSphericalDD12 * tmp34;
  //     const double tmp37 = gammaSphericalDD02 * tmp26;
  //     const double tmp47 = -tmp10 * tmp25 * xx0 - tmp11 * tmp25 * xx0;
  //     const double tmp53 = tmp10 * tmp52 + tmp11 * tmp52;
  //     const double tmp70 = 2 * KSphericalDD01 * tmp23;
  //     const double tmp73 = KSphericalDD12 * tmp34;
  //     const double tmp75 = KSphericalDD02 * tmp26;
  //     const double tmp19 = tmp10 * tmp18;
  //     const double tmp28 = gammaSphericalDD11 * tmp15 * tmp26;
  //     const double tmp31 = gammaSphericalDD22 * tmp18;
  //     const double tmp32 = tmp18 * tmp29;
  //     const double tmp39 = tmp0 * tmp18;
  //     const double tmp42 = tmp18 * tmp41;
  //     const double tmp54 = tmp18 * tmp53;
  //     const double tmp61 = tmp11 * tmp18;
  //     const double tmp64 = tmp18 * tmp6;
  //     const double tmp65 = tmp18 * ((tmp47) * (tmp47));
  //     const double tmp66 = tmp18 * ((tmp53) * (tmp53));
  //     const double tmp71 = KSphericalDD11 * tmp15 * tmp26;
  //     const double tmp72 = KSphericalDD22 * tmp18;
  //     const double tmp33 = tmp0 * tmp32;
  //     const double tmp40 = tmp39 * tmp6;
  //     const double tmp43 = tmp0 * tmp42;
  //     const double tmp44 = tmp32 * tmp6;
  //     const double tmp49 = gammaSphericalDD11 * tmp34 * tmp47;
  //     const double tmp51 = gammaSphericalDD01 * tmp26 * tmp47;
  //     const double tmp63 = tmp42 * tmp6;
  //     const double tmp77 = KSphericalDD11 * tmp34 * tmp47;
  //     const double tmp78 = KSphericalDD01 * tmp26 * tmp47;
  //     const double tmp56 = gammaSphericalDD01 * tmp34 * tmp54;
  //     const double tmp58 = gammaSphericalDD00 * tmp26 * tmp54;
  //     const double tmp79 = KSphericalDD01 * tmp34 * tmp54;
  //     const double tmp80 = KSphericalDD00 * tmp26 * tmp54;
  //     out_pw[OUTPUT_VARS::BETA1] = betaSphericalU0 * tmp0 * tmp1 + betaSphericalU1 * tmp3 * tmp4 - betaSphericalU2 * tmp6 * tmp7;
  //     out_pw[OUTPUT_VARS::BETA2] = betaSphericalU0 * tmp1 * tmp6 + betaSphericalU1 * tmp4 * tmp9 + betaSphericalU2 * tmp0 * tmp7;
  //     out_pw[OUTPUT_VARS::BETA3] = betaSphericalU0 * tmp4 - betaSphericalU1 * tmp7;
  //     out_pw[OUTPUT_VARS::G11] = tmp19 * tmp21 + tmp19 * tmp24 + tmp19 * tmp28 + ((tmp29) * (tmp29)) * tmp31 + 2 * tmp33 * tmp35 + 2 * tmp33 * tmp37;
  //     out_pw[OUTPUT_VARS::G12] =
  //         tmp21 * tmp40 + tmp24 * tmp40 + tmp28 * tmp40 + tmp29 * tmp31 * tmp41 + tmp35 * tmp43 + tmp35 * tmp44 + tmp37 * tmp43 + tmp37 * tmp44;
  //     out_pw[OUTPUT_VARS::G13] =
  //         gammaSphericalDD02 * tmp32 * tmp53 + gammaSphericalDD12 * tmp32 * tmp47 + tmp0 * tmp56 + tmp0 * tmp58 + tmp39 * tmp49 + tmp39 * tmp51;
  //     out_pw[OUTPUT_VARS::G22] = tmp21 * tmp61 + tmp24 * tmp61 + tmp28 * tmp61 + tmp31 * ((tmp41) * (tmp41)) + 2 * tmp35 * tmp63 + 2 * tmp37 * tmp63;
  //     out_pw[OUTPUT_VARS::G23] =
  //         gammaSphericalDD02 * tmp42 * tmp53 + gammaSphericalDD12 * tmp42 * tmp47 + tmp49 * tmp64 + tmp51 * tmp64 + tmp56 * tmp6 + tmp58 * tmp6;
  //     out_pw[OUTPUT_VARS::G33] = gammaSphericalDD00 * tmp66 + 2 * gammaSphericalDD01 * tmp47 * tmp54 + gammaSphericalDD11 * tmp65;
  //     out_pw[OUTPUT_VARS::K11] = tmp19 * tmp68 + tmp19 * tmp70 + tmp19 * tmp71 + ((tmp29) * (tmp29)) * tmp72 + 2 * tmp33 * tmp73 + 2 * tmp33 * tmp75;
  //     out_pw[OUTPUT_VARS::K12] =
  //         tmp29 * tmp41 * tmp72 + tmp40 * tmp68 + tmp40 * tmp70 + tmp40 * tmp71 + tmp43 * tmp73 + tmp43 * tmp75 + tmp44 * tmp73 + tmp44 * tmp75;
  //     out_pw[OUTPUT_VARS::K13] =
  //         KSphericalDD02 * tmp32 * tmp53 + KSphericalDD12 * tmp32 * tmp47 + tmp0 * tmp79 + tmp0 * tmp80 + tmp39 * tmp77 + tmp39 * tmp78;
  //     out_pw[OUTPUT_VARS::K22] = ((tmp41) * (tmp41)) * tmp72 + tmp61 * tmp68 + tmp61 * tmp70 + tmp61 * tmp71 + 2 * tmp63 * tmp73 + 2 * tmp63 * tmp75;
  //     out_pw[OUTPUT_VARS::K23] =
  //         KSphericalDD02 * tmp42 * tmp53 + KSphericalDD12 * tmp42 * tmp47 + tmp6 * tmp79 + tmp6 * tmp80 + tmp64 * tmp77 + tmp64 * tmp78;
  //     out_pw[OUTPUT_VARS::K33] = KSphericalDD00 * tmp66 + 2 * KSphericalDD01 * tmp47 * tmp54 + KSphericalDD11 * tmp65;
  //   };
  //   ADM_Spherical_to_Cart();
    
  //   out_pw[OUTPUT_VARS::ALPHA] = quant_vals[ISO_VARS::ISO_ALPHA];
  //   return out_pw;
  // }

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
    throw std::invalid_argument("\nExport: Invalid EOS Type\n)");
  }
}
