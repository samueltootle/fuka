#include "utilities.hpp"
namespace Kadath::FUKA_Solvers {

  inline NS_XCTS_BASE::NS_XCTS_BASE() :
      rank(0.), verbosity(0), ndom(-1), bconfig(nullptr), basis(nullptr), fmet(nullptr), 
        cfields(nullptr), coord_vectors(nullptr), seq(nullptr), 
          resolution(nullptr), syst(nullptr), conformal_factor(nullptr), 
            lapse(nullptr), shift(nullptr), logh(nullptr), outputdir("") {}

  inline NS_XCTS_BASE::NS_XCTS_BASE(NS_XCTS_BASE::base_config_t* config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, int const rank_) :
      rank(rank_), verbosity(0), ndom(-1), bconfig(config_), basis(nullptr), fmet(nullptr), 
        cfields(nullptr), coord_vectors(nullptr), syst(nullptr), conformal_factor(nullptr), 
          lapse(nullptr), shift(nullptr), logh(nullptr), seq(new ns_sequence(seq_)), 
            resolution(new Parameter_sequence<BCO_PARAMS>(res_)), outputdir(outputdir_) {
    
    std::array<bool, NUM_STAGES>& stage_enabled = bconfig->return_stages();
    auto [ last_stage_, last_stage_idx_ ] = get_last_enabled(MSTAGE, stage_enabled);
    last_stage_idx = last_stage_idx_;
  }

  inline void NS_XCTS_BASE::initialize_support_containers() {

    basis.reset(new Base_tensor(shift->get_basis()));
    fmet.reset(new Metric_flat(*space, *basis));

    cfields.reset(new cfgen_t(*space));
    coord_vectors = std::make_unique<cfary_t>(default_co_vector_ary(*space));
  }

  inline void NS_XCTS_BASE::save_to_file() const {
    Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift, *logh);
  }

  inline void NS_XCTS_BASE::reset_all_ptrs() {
    // Fields
    conformal_factor.reset(nullptr) ;
    lapse.reset(nullptr);
    shift.reset(nullptr);
    logh.reset(nullptr);

    // Containers
    basis.reset(nullptr);
    fmet.reset(nullptr);
    syst.reset(nullptr);
    cfields.reset(nullptr);
    coord_vectors.reset(nullptr);

    // Space
    space.reset(nullptr);
  }

  inline void NS_XCTS_BASE::regrid() {
    std::string outputfile{"ns_regrid"};
    
    if(rank == 0) {
    std::cout << "Resolution of old space: "
      << space->get_domain(0)->get_nbr_points()(0) << " (r), "
      << space->get_domain(0)->get_nbr_points()(1) << " (theta), "
      << space->get_domain(0)->get_nbr_points()(2) << " (phi)" << std::endl;
  
    int ndim = 3;
    // get the adapted domain and cast it to its correct type to be able to call its member functions
    const Domain_shell_outer_adapted* old_outer_adapted =
        dynamic_cast<const Domain_shell_outer_adapted*>(space->get_domain(1));
    
    // setup a scalar field representing the old radius
    Scalar old_space_radius(*space);
    old_space_radius = 0.;

    // get the radius from each domain
    for(int i = 0; i < space->get_nbr_domains(); ++i)
      old_space_radius.set_domain(i) = space->get_domain(i)->get_radius();
    // get the adapted radius of the adapted domain
    old_space_radius.set_domain(1) = old_outer_adapted->get_outer_radius();

    // define a standard decomposition, compatible with the parity of this field
    old_space_radius.std_base();
    //end setup old radius field

    // get the minimal and maximal radius from the adapted domain
    auto [r_min, r_max] = Kadath::bco_utils::get_rmin_rmax(*space, 1);

    std::cout << "Rmin/max: " << r_min << " " << r_max << std::endl;

    // set new resolutions in each spatial dimension
    Dim_array res(ndim);
    res.set(0) = (*bconfig)(BCO_RES);
    res.set(1) = res(0);
    res.set(2) = res(0) - 1;

    // FIXME not sure if it's only about oddness...
    if(res(0) % 2 == 0 || res(2) % 2 != 0){
      std::cout << "New Resolution is invalid.  Must be odd (9,11,13,etc)" << std::endl;
      std::_Exit(EXIT_FAILURE);
    }

    std::cout << "Resolution of new space: "
      << res(0) << " (r), "
      << res(1) << " (theta), "
      << res(2) << " (phi)" << std::endl;

    // get the type of the colocation points
    int type_coloc = space->get_type_base();

    // Update config
    bconfig->set(BCO_PARAMS::RIN)  = 0.5 * r_min;
    bconfig->set(BCO_PARAMS::ROUT) = 1.5 * r_max;
    bconfig->set(BCO_PARAMS::RMID) = r_max;
    // end update config
    
    // setup radius bounds of the domains

    int ndom = 4 + (*bconfig)(BCO_PARAMS::NSHELLS);
    std::vector<double> bounds(ndom-1);
    Kadath::bco_utils::set_NS_bounds(bounds, *bconfig);
    
    Kadath::bco_utils::print_bounds("New bounds: ", bounds);

    // get origin of nucleus domain
    Point center = space->get_domain(0)->get_center();

    // initialize space with new resolution and domain decomposition
    Space_spheric_adapted new_space(type_coloc, center, res, bounds);
    Base_tensor basis(new_space, CARTESIAN_BASIS);

    // get adapted domains to update the radius
    const Domain_shell_outer_adapted* new_outer_adapted = dynamic_cast<const Domain_shell_outer_adapted*>(new_space.get_domain(1));
    const Domain_shell_inner_adapted* new_inner_adapted = dynamic_cast<const Domain_shell_inner_adapted*>(new_space.get_domain(2));

    // update adapted domain mapping
    Kadath::bco_utils::interp_adapted_mapping(new_outer_adapted, 1, old_space_radius);
    Kadath::bco_utils::interp_adapted_mapping(new_inner_adapted, 1, old_space_radius);

    // setup new fields
    // initialize to one or zero first
    Scalar new_conf(new_space);
    new_conf = 1.;
    new_conf.std_base();

    Scalar new_lapse(new_space);
    new_lapse = 1.;
    new_lapse.std_base();

    Vector new_shift(new_space, CON, basis);
    for (int i = 1; i <= 3; i++)
      new_shift.set(i).annule_hard();
    new_shift.std_base();

    Scalar new_logh(new_space);
    new_logh.annule_hard();
    new_logh.std_base();

    // end setup new fields
    
    // import data from fields in the old space
    new_conf.import(*conformal_factor);
    new_lapse.import(*lapse);
    new_logh.import(*logh);

    new_shift.set(1).import(shift->set(1));
    new_shift.set(2).import(shift->set(2));
    new_shift.set(3).import(shift->set(3));

    // end import old fields

    // enforce spectral decomposition compatible with the parities
    new_lapse.std_base();
    new_conf.std_base();
    new_logh.std_base();
    new_shift.std_base();
    
    // output data  
    bconfig->set_filename(outputfile);
    Kadath::bco_utils::save_to_file(new_space, *bconfig, new_conf, new_lapse, new_shift, new_logh);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Ensure all ranks have the same config file
    bconfig->set_filename(outputfile);
    bconfig->open_config();
    
    // Update stored fields and containers
    reset_all_ptrs();
    load_solution_from_file();
    initialize_support_containers();
  }

  void NS_XCTS_BASE::load_solution_from_file() {
    std::string spacein{bconfig->space_filename()};
    FILE* ff1 = fopen (spacein.c_str(), "r") ;

    space.reset(new base_space_t{ff1});
    conformal_factor.reset( new Scalar(*space, ff1)) ;
    lapse.reset( new Scalar(*space, ff1)) ;
    shift.reset( new Vector(*space, ff1)) ;
    logh.reset( new Scalar(*space, ff1)) ;    
    fclose(ff1);
        
    ndom = space->get_nbr_domains();
  }

  bool NS_XCTS_BASE::increment_resolution() {
    if(!(resolution->final() > resolution->init()) || last_stage_idx != solver_stage)
      return false;
    
    auto resolution_indices = resolution->get_indices();
    auto const & final_res = resolution->final();
    if((*bconfig)(resolution_indices) > final_res)
      return false;

    int next_res = bco_utils::next_resolution((*bconfig)(resolution_indices));

    // Greater would mean that the desired resolution may not be possible
    // (see bco_utils::next_resolution) so we go to the next available
    // resolution
    if(next_res >= final_res) {
      bconfig->set(resolution_indices) = next_res;
    } else {
      bconfig->set(resolution_indices) = next_res;
    }
    return true;
  }

  bool NS_XCTS_BASE::increment_seq() {
    if(!seq->is_set() || last_stage_idx != solver_stage)
      return false;
    
    auto sequence_var_indices = seq->get_indices();
    auto const & dx = seq->step_size();
    auto x = bconfig->set(sequence_var_indices) + dx;
    if(seq->loop_condition(x)) {
      bconfig->set(sequence_var_indices) = x;
      return true;
    }
    return false;
  }

  int NS_XCTS_BASE::do_newton() {
    int exit_status = EXIT_SUCCESS;
    // parameters for the solver loop
    bool endloop = false;
    int ite = 1;
    double conv;
  
    // solve until convergence is achieved
    while (!endloop) {  
      // do exactly one newton step, given the system above
      endloop = syst->do_newton(bconfig->seq_setting(SEQ_SETTINGS::PREC), conv);
  
      update_config_quantities();
      // output files at this iteration and print diagnostics
      std::stringstream ss;
      ss << stagename
         << "_ckpt_" << ite - 1;
      
      bconfig->set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(ite, conv);
        std::cout << std::endl;
        if(bconfig->control(CHECKPOINT))
          checkpoint();
      }
  
      // update all coordinate fields, in case the domain extents have changed
      update_fields_co(*cfields, *coord_vectors, {}, 0.);

      ite++;
      check_max_iter_exceeded(*this, ite, conv);
    }
  
    bconfig->set_filename(converged_filename(stagename));
    if (rank == 0) {
      checkpoint();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return exit_status;
  }
}