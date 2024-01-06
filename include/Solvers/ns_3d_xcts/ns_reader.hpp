#include "Solvers/reader.hpp"
namespace Kadath::FUKA_Solvers {
template<class eos_t, class config_t, class space_t>
struct CFMS_NS_Reader : public Reader<config_t, space_t> {
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

  using pointwise_ary_t = std::vector<double>; 
  using grid_ary_t = std::array<pointwise_ary_t, OUTPUT_VARS::NUM_OUTPUT_VARS>;

  // Types
  using Reader<config_t, space_t>::base_space_t;
  using Reader<config_t, space_t>::base_config_t;

  // Reuse base class members
  using Reader<config_t, space_t>::space;
  using Reader<config_t, space_t>::bconfig;
  using Reader<config_t, space_t>::syst;
  using Reader<config_t, space_t>::ndom;

  // CFMS_BH imported fields from file
  ptr_data_member(Scalar, conformal_factor, shared);
  ptr_data_member(Scalar, lapse, shared);
  ptr_data_member(Vector, shift, shared);
  ptr_data_member(Scalar, logh, shared);
  ptr_data_member(Vector, fluidvel, shared);

  // Constructed objects
  ptr_data_member(Base_tensor, basis, shared);
  ptr_data_member(Metric_flat, fmet, shared);
  ptr_data_member(Tensor, A, shared);

  protected:
  std::vector<std::reference_wrapper<const Scalar>> quants;
  std::vector<double> quant_vals;
  pointwise_ary_t out_pw;
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
    fclose(ff1);
    
    basis.reset(new Base_tensor{shift->get_basis()});
    fmet.reset(new Metric_flat(*space, *basis));
    ndom = space->get_nbr_domains();
  }
  
  void extract_computed_grid_functions() {
    syst.reset(new System_of_eqs (*space));
    fmet->set_system(*syst, "f") ;
    
    // fields depending on the coords
    CoordFields<Space_spheric_adapted> cf_generator(*space);
    vec_ary_t coord_vectors {default_co_vector_ary(*space)};

    // get origin of the system and initialize coordinate fields
    double xo = Kadath::bco_utils::get_center(*space,0);
    update_fields_co(cf_generator, coord_vectors, {}, xo);

    Param p;
    syst->add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst->add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst->add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
    // end adding EOS OPEs

    // Fields - must be initialized before common setup
    syst->add_cst("N"  , *lapse) ;
    syst->add_cst("bet", *shift) ;
    syst->add_cst("ome" , (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA));
    syst->add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
    syst->add_def("omega^i = bet^i + ome * mg^i");

    syst->add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    A.reset(new Tensor(syst->give_val_def("A")));
    A->coef();
  
    // definitions for the fluid 3-velocity
    syst->add_def("U^i = omega^i / N");
    fluidvel.reset(new Vector(syst->give_val_def("U")));
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
  bool is_export_ready() const { return export_ready; }
  std::vector<std::reference_wrapper<const Scalar>> const & get_quants() const { return quants; }
  const int & get_ndim() const { return ndim; }
  
  CFMS_NS_Reader() : Reader<config_t, space_t>(),
    basis(nullptr), fmet(nullptr),
    conformal_factor(nullptr), lapse(nullptr), shift(nullptr), logh(nullptr), fluidvel(nullptr) {}
  CFMS_NS_Reader(std::string config_filename) :
    Reader<config_t, space_t>(config_filename),
      basis(nullptr), fmet(nullptr),
        conformal_factor(nullptr), lapse(nullptr), shift(nullptr), logh(nullptr), fluidvel(nullptr) {
  
    quant_vals.resize(XCTS_VARS::NUM_XCTS_VARS);
    out_pw.resize(OUTPUT_VARS::NUM_OUTPUT_VARS);

    load_solution_from_file();
    extract_computed_grid_functions();
    populate_quants();
    
    // This is to avoid a "bug" where "something" in kadath is not
    // correctly initialized prior to copying to other threads resulting
    // in undefined behavior.  By running the interpolator once, this
    // bug seems to be avoided.
    this->export_pointwise(0.5, 0., 0.);

    //for(auto& e : out_pw)
    //  cout << e << ", ";
    //cout << endl;
  }

  CFMS_NS_Reader(CFMS_NS_Reader const & r) : fmet(nullptr) {
    std::lock_guard<std::mutex> lock(copy_mutex);
    ndom = r.ndom;
    
    quant_vals.resize(XCTS_VARS::NUM_XCTS_VARS);
    out_pw.resize(OUTPUT_VARS::NUM_OUTPUT_VARS);

    space.reset(new space_t((*r.get_space())));
    conformal_factor.reset(new Scalar(*space.get(), *r.conformal_factor.get()));
    lapse.reset(new Scalar(*space, *r.lapse.get()));
    shift.reset(new Vector(*space, *r.shift.get()));
    logh.reset(new Scalar(*space, *r.logh.get()));
    fluidvel.reset(new Vector(*space, *r.fluidvel.get()));
    A.reset(new Tensor(*space, *r.A.get()));

    basis.reset(new Base_tensor(*space, r.get_basis()->get_basis(0)));
    fmet.reset(new Metric_flat(*space, *basis));
    bconfig.reset(new config_t(*r.bconfig));
    
    export_ready = false;

    populate_quants();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
    // std::cout << "copy\n";
  }

  CFMS_NS_Reader(CFMS_NS_Reader&& b) noexcept = delete;
  CFMS_NS_Reader& operator=(const CFMS_NS_Reader& b)
  {
    if (this == &b) return *this;

    CFMS_NS_Reader tmp(b);
    *this = std::move(tmp);
    return *this;
  }

  public:

  std::vector<double> interpolate_pointwise(double const & x, double const & y, double const & z) {
    
    Point abs_coords(ndim);
    abs_coords.set(1) = x;
    abs_coords.set(2) = y;
    abs_coords.set(3) = z;
    
    // For testing only
    for (size_t k = 0; k < XCTS_VARS::NUM_XCTS_VARS; ++k) {
        quant_vals[k] = quants[k].get().val_point(abs_coords);
    }
    
    return quant_vals;
  }

  
  pointwise_ary_t export_pointwise(
    double const & x, double const & y, double const & z) {
    
    if(quant_vals.size() != XCTS_VARS::NUM_XCTS_VARS)
      quant_vals.resize(XCTS_VARS::NUM_XCTS_VARS);
    
    if(out_pw.size() != OUTPUT_VARS::NUM_OUTPUT_VARS)
      out_pw.resize(OUTPUT_VARS::NUM_OUTPUT_VARS);
    
    quant_vals = interpolate_pointwise(x, y, z);
     
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
      rho = EOS<eos_t, DENSITY>::get(h);
      eps = EOS<eos_t, EPSILON>::get(h);
      press = EOS<eos_t, PRESSURE>::get(h);
    }
    out_pw[OUTPUT_VARS::RHO]   = rho;
    out_pw[OUTPUT_VARS::EPS]   = eps;
    out_pw[OUTPUT_VARS::PRESS] = press;
    out_pw[OUTPUT_VARS::VELX]  = quant_vals[XCTS_VARS::XCTS_UX];
    out_pw[OUTPUT_VARS::VELY]  = quant_vals[XCTS_VARS::XCTS_UY];
    out_pw[OUTPUT_VARS::VELZ]  = quant_vals[XCTS_VARS::XCTS_UZ];
    return out_pw;
  }

  std::array<std::vector<double>, OUTPUT_VARS::NUM_OUTPUT_VARS> export_coordinate_array(
    int const npoints, double const * xx, double const * yy, double const * zz) {
    
    std::array<std::vector<double>,OUTPUT_VARS::NUM_OUTPUT_VARS> out;
    for(auto& v : out)
      v.resize(npoints);
    
    for (size_t i = 0; i < npoints; ++i) {
      out_pw = export_pointwise(xx[i], yy[i], zz[i]);

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
    }
    return out;
  }
};
}
