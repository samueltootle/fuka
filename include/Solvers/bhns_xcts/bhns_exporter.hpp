#include "Solvers/exporter.hpp"
#include "bhns.hpp"
namespace Kadath::FUKA_Solvers {

template<class eos_t>
struct CFMS_BHNS_Exporter : public Exporter<Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BIN_INFO>, Space_bhns> {
  using config_t = Kadath::FUKA_Config::kadath_config_boost<Kadath::FUKA_Config::BIN_INFO>;
  using space_t = Space_bhns;

  enum XCTS_VARS : size_t {  
    XCTS_PSI,
    XCTS_ALPHA,
    XCTS_BETAX,
    XCTS_BETAY,
    XCTS_BETAZ,
    XCTS_AXX,
    XCTS_AXY,
    XCTS_AXZ,
    XCTS_AYY,
    XCTS_AYZ,
    XCTS_AZZ,
    XCTS_H,
    XCTS_UX,
    XCTS_UY,
    XCTS_UZ,
    NUM_XCTS_VARS 
  };

  enum OUTPUT_VARS : size_t {
    ALPHA,
    BETAX,
    BETAY,
    BETAZ,
    GXX,
    GXY,
    GXZ,
    GYY,
    GYZ,
    GZZ,
    KXX,
    KXY,
    KXZ,
    KYY,
    KYZ,
    KZZ,
    RHO,
    EPS,
    PRESS,
    VELX,
    VELY,
    VELZ,
    NUM_OUTPUT_VARS
  };

  using interp_ary_t = std::array<double, NUM_XCTS_VARS>; 
  using output_ary_t = std::array<double, NUM_OUTPUT_VARS>; 
  using grid_ary_t = std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  // Types
  using Exporter<config_t, space_t>::base_space_t;
  using Exporter<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Exporter<config_t, space_t>::space;
  using Exporter<config_t, space_t>::bconfig;
  using Exporter<config_t, space_t>::ndom;

  // CFMS_BHNS imported fields from file
  ptr_data_member(Scalar, conformal_factor, shared);
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Vector, shift, shared);
  ptr_data_member(Scalar, logh, shared);
  ptr_data_member(Scalar, velpotential, shared);

  // Constructed objects
  ptr_data_member(Tensor, A, shared);
  ptr_data_member(Vector, fluidvel, shared);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  interp_ary_t quant_vals;
  output_ary_t out_pw;
  bool export_ready{false};
  int const ndim{3};

  void load_solution_from_file() override {
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
  
  void extract_computed_grid_functions() {
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
    syst.add_ope("eps"  , &EOS<eos_t, eos_var_t::EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, eos_var_t::PRESSURE>::action, &p);
    syst.add_ope("rho"  , &EOS<eos_t, eos_var_t::DENSITY>::action, &p);
    // end adding EOS OPEs

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

  void populate_quants() {
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

  public:
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  bool is_export_ready() const { return export_ready; }  
  const int & get_ndim() const { return ndim; }
  
  CFMS_BHNS_Exporter() : Exporter<config_t, space_t>(),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr), 
      logh(nullptr), velpotential(nullptr), fluidvel(nullptr) {}
  
  CFMS_BHNS_Exporter(std::string config_filename) : Exporter<config_t, space_t>(config_filename),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr), 
      logh(nullptr), velpotential(nullptr), fluidvel(nullptr) {

    load_solution_from_file();
    extract_computed_grid_functions();
    populate_quants();
    
    // This is to avoid a "bug" where "something" in kadath is not
    // correctly initialized prior to copying to other threads resulting
    // in undefined behavior.  By running the interpolator once, this
    // bug seems to be avoided.
    this->export_pointwise(0.5, 0., 0.);
  }

  CFMS_BHNS_Exporter(CFMS_BHNS_Exporter const & r) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;

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

  CFMS_BHNS_Exporter(CFMS_BHNS_Exporter&& b) noexcept = delete;
  CFMS_BHNS_Exporter& operator=(const CFMS_BHNS_Exporter& b)
  {
    if (this == &b) return *this;

    CFMS_BHNS_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
  }

  public:

  interp_ary_t interpolate_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {

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

  
  output_ary_t export_pointwise(double const & x, double const & y, double const & z,
    double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
    
    quant_vals = interpolate_pointwise(x, y, z, interpolation_offset, interp_order, delta_r_rel);
     
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

    double const H = quant_vals[XCTS_VARS::XCTS_H];
    double h = std::exp(H);
    double rho, eps, press;

    // get quantities point-wise, since h is smoothest, and cut data at H=0
    if(std::fabs(H) <= 1e-12) {
      rho = 0.;
      eps = 0.;
      press = 0.;
    }
    else {
      rho = EOS<eos_t, eos_var_t::DENSITY>::get(h);
      eps = EOS<eos_t, eos_var_t::EPSILON>::get(h);
      press = EOS<eos_t, eos_var_t::PRESSURE>::get(h);
    }
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
    out_pw[OUTPUT_VARS::VELX]  = quant_vals[XCTS_VARS::XCTS_UX];
    out_pw[OUTPUT_VARS::VELY]  = quant_vals[XCTS_VARS::XCTS_UY];
    out_pw[OUTPUT_VARS::VELZ]  = quant_vals[XCTS_VARS::XCTS_UZ];
    return out_pw;
  }

  grid_ary_t export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz,
      double const interpolation_offset = 0., int const interp_order = 8, double const delta_r_rel = 0.3) {
    
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
};
}
