#include "Solvers/ns_3d_xcts/ns_exporter.hpp"
namespace Kadath::FUKA_Solvers {
  CFMS_NS_Exporter::CFMS_NS_Exporter(CFMS_NS_Exporter const & r) {
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
    fluidvel.reset(new Vector(*space, *r.fluidvel.get()));
    A.reset(new Tensor(*space, *r.A.get()));

    bconfig.reset(new config_t(*r.bconfig));
    
    export_ready = false;

    populate_quants();
    // For testing only
    // Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift);
    // std::cout << "copy\n";
  }

  CFMS_NS_Exporter& CFMS_NS_Exporter::operator=(const CFMS_NS_Exporter& b) {
    if (this == &b) return *this;

    CFMS_NS_Exporter tmp(b);
    *this = std::move(tmp);
    return *this;
  }

  void CFMS_NS_Exporter::initialize_eos() {
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
    }

    if(eos_type == "Cold_PWPoly") {
      using namespace Kadath::Margherita;
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
    } // end adding EOS OPEs
  }

  void CFMS_NS_Exporter::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;

    space.reset(new space_t{ff1});
    conformal_factor.reset( new Scalar(*space.get(), ff1)) ;
    lapse.reset( new Scalar(*space.get(), ff1)) ;
    shift.reset( new Vector(*space.get(), ff1)) ;
    logh.reset( new Scalar(*space.get(), ff1)) ;
    
    fclose(ff1);
        
    ndom = space->get_nbr_domains();
  }

  void CFMS_NS_Exporter::extract_computed_grid_functions() {
    System_of_eqs syst(*space);
    Base_tensor basis(shift->get_basis());
    Metric_flat fmet(*space, basis);
    fmet.set_system(syst, "f") ;
    
    // fields depending on the coords
    CoordFields<Space_spheric_adapted> cf_generator(*space);
    vec_ary_t coord_vectors {default_co_vector_ary(*space)};

    // get origin of the system and initialize coordinate fields
    double xo = Kadath::bco_utils::get_center(*space,0);
    update_fields_co(cf_generator, coord_vectors, {}, xo);

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
    syst.add_cst("ome" , (*bconfig)(Kadath::FUKA_Config::BCO_PARAMS::OMEGA));
    syst.add_cst("mg"  , *coord_vectors[GLOBAL_ROT]);
    syst.add_def("omega^i = bet^i + ome * mg^i");

    syst.add_def("A_ij = (D_i bet_j + D_j bet_i - 2. / 3.* D^k bet_k * f_ij) /2. / N");
    A.reset(new Tensor(syst.give_val_def("A")));
    A->coef();
  
    // definitions for the fluid 3-velocity
    syst.add_def("U^i = omega^i / N");
    fluidvel.reset(new Vector(syst.give_val_def("U")));
    fluidvel->coef();
  }

  void CFMS_NS_Exporter::populate_quants() {
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

  CFMS_NS_Exporter::interp_ary_t CFMS_NS_Exporter::interpolate_pointwise(double const & x, double const & y, double const & z) {
    
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


  CFMS_NS_Exporter::output_ary_t CFMS_NS_Exporter::export_pointwise(double const & x, double const & y, double const & z) {
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_pointwise_imp<eos_t>(x, y, z);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_pointwise_imp<eos_t>(x, y, z);
    } // end adding EOS OPEs
    throw std::invalid_argument("\nInvalid EOS Type\n)");
  }

  CFMS_NS_Exporter::grid_ary_t CFMS_NS_Exporter::export_coordinate_array(int const npoints, 
    double const * xx, double const * yy, double const * zz) {
    
    if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;
      return export_coordinate_array_imp<eos_t>(npoints, xx, yy, zz);
    }else if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;
      return export_coordinate_array_imp<eos_t>(npoints, xx, yy, zz);
    } // end adding EOS OPEs
    throw std::invalid_argument("\nInvalid EOS Type\n)");
  }
}