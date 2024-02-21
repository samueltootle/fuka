#include "Solvers/fuka_syst/fuka_syst_setup.hpp"
#include "utilities.hpp"
namespace Kadath::FUKA_Solvers {
  // NOROT Routines
  template<class eos_t>
  NS_XCTS_UNIFORM_ROT<eos_t>::NS_XCTS_UNIFORM_ROT(NS_XCTS_BASE::base_config_t& config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, int const rank_) :
      NS_XCTS_BASE(config_, seq_, res_, outputdir_, rank_) {
    
    solver_stage = ::Kadath::FUKA_Config::STAGES::NOROT_BC;
    if(!seq->is_set() && !bconfig->control(CONTROLS::SEQUENCES)) {
      initialize_config_from_fixing_values(*bconfig, *seq);
    }
    load_solution_from_file();
    initialize_EOS(*this);
    initialize_support_containers();

    if(rank == 0)
      cout << *seq << endl;
  }
  
  template<class eos_t>
  std::string NS_XCTS_UNIFORM_ROT<eos_t>::converged_filename(const std::string stage) const {
    // FIXME assumes a fix resolution for all domains
    auto res = space->get_domain(0)->get_nbr_points()(0);
    const std::string eosname{extract_eos_name(*bconfig)};
    std::stringstream ss;
    ss << "NS";
    if(stage != "") ss  << "_" << stage << ".";
    else ss << ".";
    ss << eosname << ".";
    
    // Add mass fixing parameter to filename
    auto default_idx = BCO_PARAMS::MADM;
    if(seq) {
      update_filename_from_mass_fixing(*bconfig, seq, ss);
    } else {
      auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
      ss << seq_key << "." << (*bconfig)(default_idx) << "."; 
    }  
    
    default_idx = BCO_PARAMS::CHI;
    if(seq) {
      update_filename_from_spin_fixing(*bconfig, seq, ss);
    } else {
      auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
      ss << seq_key << "." << (*bconfig)(default_idx) << "."; 
    }
    ss << (*bconfig)(BCO_PARAMS::NSHELLS) << "."
       <<std::setfill('0') << std::setw(2) << res;
    return ss.str();
  }

  template<class eos_t>
  void NS_XCTS_UNIFORM_ROT<eos_t>::setup_syst() {
    int exit_status = EXIT_SUCCESS;
    double loghc = std::log((*bconfig)(BCO_PARAMS::HC));

    // We use `config_filename()` vs `config_filename_abs()` since
    // `solution_exists` will probe the HOME_KADATH/COs directory
    // auto const current = bconfig.config_filename();
    // if(!bconfig.control(RESOLVE) && solution_exists(stagename)) {    
    //   if(rank == 0)
    //     std::cout << "Solved previously: " \
    //               << bconfig.config_filename_abs() << std::endl;
    //   return (current == bconfig.config_filename()) ? \
    //     EXIT_SUCCESS : RELOAD_FILE;
    // }

    update_fields_co(*cfields, *coord_vectors,{}, 0.);
    syst.reset(new System_of_eqs(*space));
    syst_init();
    
    std::string central_fixing_definition{"H - Hc"};
    std::string spin_fixing_definition{"integ(intJ) - chi * Madm * Madm = 0"};
    
    if(seq) {
      central_fixing_definition = ::Kadath::FUKA_Syst_tools::set_ns_mass_fixing(*syst, *bconfig, seq);
      spin_fixing_definition = ::Kadath::FUKA_Syst_tools::set_ns_spin_fixing(*syst, *bconfig, seq);
    } else {
      syst->add_var("Hc", loghc);
      syst->add_cst("chi" , (*bconfig)(BCO_PARAMS::CHI));
      syst->add_var("ome" , (*bconfig)(BCO_PARAMS::OMEGA));
      syst->add_cst("Madm", (*bconfig)(BCO_PARAMS::MADM));
    }
    syst->add_def("omega^i = bet^i + ome * mg^i");

    for (int d = 0; d < ndom; d++) {
      switch (d) {
      case 0:
      case 1:
        syst->add_def(d, "U^i = omega^i / N");
        syst->add_def(d, "Usquare = P^4 * U_i * U^i");
        syst->add_def(d, "Wsquare = 1. / (1. - Usquare)");
        syst->add_def(d, "W = sqrt(Wsquare)");

        syst->add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
        syst->add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare") ;
        syst->add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

        syst->add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde") ;
        syst->add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                              "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst->add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                              "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

        syst->add_def(d, "intMb = P^6 * rho(h) * W");
        syst->add_def(d, "firstint = H + log(N) - log(W)");

        break;
      default:
        syst->add_eq_full(d, "H = 0");

        syst->add_def(d, "eqP = D^i D_i P + A_ij * A^ij / P^7 / 8");
        syst->add_def(d, "eqNP = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
        syst->add_def(d, "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                        "A^ij * D_j Ntilde");
        break;
      }
    }

    if (rank == 0) {
      std::cout << "############################" << std::endl
                << "Uniformly Rotating NS Solver" << std::endl
                << "############################" << std::endl;
    }

    // add the constraint equations and demand continuity their normal derivative across domain boundaries
    space->add_eq(*syst, "eqNP= 0", "N", "dn(N)");
    space->add_eq(*syst, "eqP = 0", "P", "dn(P)");
    space->add_eq(*syst, "eqbet^i= 0", "bet^i", "dn(bet^i)");
    
    // boundary conditions at infinity
    syst->add_eq_bc(ndom - 1, OUTER_BC, "N=1");
    syst->add_eq_bc(ndom - 1, OUTER_BC, "P=1");
    syst->add_eq_bc(ndom - 1, OUTER_BC, "bet^i=0");

    // if the surface is resolved, define it to be where the matter vanishes
    syst->add_eq_bc(1, OUTER_BC, "H = 0");

    // first integral in the innermost domains with non-zero matter content
    // and condition on the central value, either fixed directly or by the
    // integral below
    syst->add_eq_first_integral(0, 1, "firstint", central_fixing_definition.c_str());
  
    // constrain stellar mass by these integrals and the central log enthalpy
    if(seq) {
      auto idx{seq->mass_idx()};
      bool add_Madm_int = true;
      switch(idx) {
        case BCO_PARAMS::MADM:
          space->add_eq_int_inf(*syst, "integ(intMadm) = Madm");
          add_Madm_int = false;
          break;
        case BCO_PARAMS::MB:
          space->add_eq_int_volume(*syst, 2, "integvolume(intMb) = Mb");
          break;
        default:
          break;
      }
      idx = seq->spin_idx();
      switch(idx) {
        case BCO_PARAMS::JADM:
          space->add_eq_int_inf(*syst, spin_fixing_definition.c_str());
          break;
        case BCO_PARAMS::CHI:
          // Since we need MADM to compute CHI, we need to ensure
          // that if it isn't a fixed quantity that it becomes a
          // variable in our system of equations and the appropriate
          // constraint equation is added
          if(add_Madm_int) {
            syst->add_var("Madm", (*bconfig)(BCO_PARAMS::MADM));
            space->add_eq_int_inf(*syst, "integ(intMadm) = Madm");
          }
          
          space->add_eq_int_inf(*syst, spin_fixing_definition.c_str());
          break;
        default:
          break;
      }
    } else {
    space->add_eq_int_inf(*syst, spin_fixing_definition.c_str());
    space->add_eq_int_inf(*syst, "integ(intMadm) = Madm");
    }
  }

  template<class eos_t>
  int NS_XCTS_UNIFORM_ROT<eos_t>::do_newton() {
    int exit_status = EXIT_SUCCESS;
    std::string stagename = "NOROT_BC";
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
      ss << "rot_3d_"
         << "_" << ite - 1;
      
      bconfig->set_filename(ss.str());
      if (rank == 0) {
        print_diagnostics(ite, conv);
        if(bconfig->control(CHECKPOINT))
          checkpoint();
      }
  
      // update all coordinate fields, in case the domain extents have changed
      update_fields_co(*cfields, *coord_vectors, {}, 0., &(*syst));

      ite++;
      check_max_iter_exceeded(*this, ite, conv);
    }
  
    bconfig->set_filename(converged_filename(stagename));
    if (rank == 0) {
      checkpoint();
    }
    return exit_status;
  }

  template<class eos_t>
  void NS_XCTS_UNIFORM_ROT<eos_t>::syst_init() {
    using namespace ::Kadath::Margherita;
    
    // call the (flat) conformal metric "f"
    fmet->set_system(*syst, "f");
  
    // define numerical constants
    syst->add_cst("4piG", (*bconfig)(BCO_QPIG));
    
    // include the coordinate fields
    syst->add_cst("mg"  , *(*coord_vectors)[GLOBAL_ROT]);
    syst->add_cst("sm"  , *(*coord_vectors)[S_BCO1]);
    syst->add_cst("einf", *(*coord_vectors)[S_INF]);
    
    // the basic fields, conformal factor, lapse and (log) enthalpy
    syst->add_var("P"   , *conformal_factor);
    syst->add_var("N"   , *lapse);
    syst->add_var("H"   , *logh);
    syst->add_var("bet" , *shift);
    
    // define common combinations of conformal factor and lapse
    syst->add_def("NP = P*N");
    syst->add_def("Ntilde = N / P^6");
    syst->add_def("A^ij = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / "
                  "2. / Ntilde");
  
    // define quantity to be integrated at infinity
    // two (in this case) equivalent definitions of ADM mass
    // as well as the Komar mass
    syst->add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
    syst->add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
    syst->add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");

    // ADM Angular momentum
    syst->add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");
    // Quasi-local spin angular momentum
    syst->add_def(2,"intS = A_ij * mg^i * sm^j / 2. / 4piG") ;
    
    // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
    syst->add_def("h = exp(H)");
  
    // define the EOS operators
    Param p;
    set_eos_ope<eos_t>(*syst, p);
  
    // define rest-mass density, internal energy and pressure through the enthalpy
    syst->add_def("rho = rho(h)");
    syst->add_def("eps = eps(h)");
    syst->add_def("press = press(h)");
    syst->add_def("dHdlnrho = dHdlnrho(h)");

    // definition to rescale the equations
    // delta = p / rho
    syst->add_def("delta = h - eps - 1.");
  }

  template<class eos_t>
  void NS_XCTS_UNIFORM_ROT<eos_t>::print_diagnostics(const int ite, const double conv) const {

    // compute the baryonic mass at volume integral from the given integrant
    double baryonic_mass =
        syst->give_val_def("intMb")()(0).integ_volume() +
        syst->give_val_def("intMb")()(1).integ_volume();

    // compute the ADM mass as surface integral at infinity  
    Val_domain integMadm(syst->give_val_def("intMadm")()(ndom - 1));
    double Madm = space->get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

    // compute the Komar mass as surface integral at infinity
    Val_domain integMk(syst->give_val_def("intMk")()(ndom - 1));
    double Mk = space->get_domain(ndom - 1)->integ(integMk, OUTER_BC);

    // get the maximum and minimum coordinate radius along the surface,
    // i.e. the adapted domain boundary
    auto rs = Kadath::bco_utils::get_rmin_rmax(*space, 1);

    // alternative, equivalent ADM mass integral
    Val_domain integMadmalt(syst->give_val_def("intMadmalt")()(ndom - 1));
    double Madmalt = space->get_domain(ndom - 1)->integ(integMadmalt, OUTER_BC);

    // compute the ADM angular momentum as surface integral at infinity  
    Val_domain integJ(syst->give_val_def("intJ")()(ndom - 1));
    double J = space->get_domain(ndom - 1)->integ(integJ, OUTER_BC);

    // output to standard output  
    #define FORMAT std::setw(13) << std::left << std::showpos
    std::ios_base::fmtflags f( std::cout.flags() );
    std::cout << "=======================================" << std::endl
              << FORMAT << "Iter: " << ite << std::endl
              << FORMAT << "Error: " << conv << std::endl
              << FORMAT << "Mb: " << baryonic_mass << std::endl
              << FORMAT << "Madm: " << Madm << std::endl
              << FORMAT << "Madm_ql: " << Madmalt 
              << " [" << std::abs(Madm - Madmalt) / Madm << "]" << std::endl
              << FORMAT << "Mk: " << Mk << " [" 
              << std::abs(Madm - Mk) / Madm << "]" << std::endl;
    std::cout << FORMAT << "R: " << rs[0] << " " << rs[1];
    std::cout << FORMAT << "Jadm: " << J << std::endl
              << FORMAT << "Chi: " << J / Madm / Madm << " [" << (*bconfig)(CHI) << "]\n"
              << FORMAT << "Omega: " << (*bconfig)(OMEGA) << std::endl;
    std::cout.flags(f);
    #undef FORMAT
    std::cout << "=======================================" << "\n\n";
  } // end print diagnostics norot

  template<class eos_t>
  void NS_XCTS_UNIFORM_ROT<eos_t>::update_config_quantities() {

    auto rs = bco_utils::get_rmin_rmax(*space, 1);
    bconfig->set(BCO_PARAMS::RMID) = rs[0];

    // compute the ADM mass as surface integral at infinity
    Val_domain integMadm(syst->give_val_def("intMadm")()(ndom - 1));
    double Madm = space->get_domain(ndom - 1)->integ(integMadm, OUTER_BC);

    // compute the baryonic mass at volume integral from the given integrant
    double baryonic_mass =
      syst->give_val_def("intMb")()(0).integ_volume() +
      syst->give_val_def("intMb")()(1).integ_volume();

    auto loghc = bco_utils::get_boundary_val(0, *logh, INNER_BC);

    // compute the ADM Angular Momentum as surface integral at infinity    
    Val_domain integJ(syst->give_val_def("intJ")()(ndom - 1));
    double const Jadm = space->get_domain(ndom - 1)->integ(integJ, OUTER_BC);
    double const chi = Jadm / Madm / Madm;

    if(seq) {
      auto idx{seq->mass_idx()};
      switch(idx) {
        case BCO_PARAMS::HC:
          bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
          bconfig->set(BCO_PARAMS::MADM) = Madm;
          bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
          break;
        case BCO_PARAMS::NC:
          bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
          bconfig->set(BCO_PARAMS::MADM) = Madm;
          bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
          break;
        case BCO_PARAMS::MB:
          bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
          bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
          bconfig->set(BCO_PARAMS::MADM) = Madm;
          break;
        default:
          bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
          bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
          bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
          break;
      }

      idx = seq->spin_idx();
      switch(idx) {
        case BCO_PARAMS::CHI:
          break;
        default:
          bconfig->set(BCO_PARAMS::CHI) = chi;
          break;
      }
    } else {
      bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
      bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
      bconfig->set(BCO_PARAMS::MB) = baryonic_mass;
      bconfig->set(BCO_PARAMS::CHI) = chi;
    }
    bconfig->set(BCO_PARAMS::QLMADM) = bconfig->set(BCO_PARAMS::MADM) ;
  }

  template<class eos_t>
  void NS_XCTS_UNIFORM_ROT<eos_t>::load_solution_from_file() {
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

  template<class eos_t>
  bool NS_XCTS_UNIFORM_ROT<eos_t>::increment_spin() {
    if(!spinup || !spinup->is_set())
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

  template<class eos_t>
  bool NS_XCTS_UNIFORM_ROT<eos_t>::increment_seq() {
    if(!seq->is_set() || last_stage_idx != ::Kadath::FUKA_Config::STAGES::NOROT_BC)
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

  template<class eos_t>
  bool NS_XCTS_UNIFORM_ROT<eos_t>::increment_resolution() {
    if(!(resolution->final() > resolution->init()) || last_stage_idx != ::Kadath::FUKA_Config::STAGES::NOROT_BC)
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
    
  template<class eos_t>
  void NS_XCTS_UNIFORM_ROT<eos_t>::regrid() {
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
}