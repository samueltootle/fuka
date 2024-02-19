#include "utilities.hpp"
namespace Kadath::FUKA_Solvers {
  // NOROT Routines
  template<class eos_t>
  std::string NS_XCTS_NOROT<eos_t>::converged_filename(const std::string stage) const {
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
    
    ss << "0.";
    ss << (*bconfig)(BCO_PARAMS::NSHELLS) << "."
       <<std::setfill('0') << std::setw(2) << res;
    return ss.str();
  }

  template<class eos_t>
  int NS_XCTS_NOROT<eos_t>::solve(bool fixed) {
    int exit_status = EXIT_SUCCESS;
    double loghc = std::log((*bconfig)(BCO_PARAMS::HC));
    std::string stagename = (fixed) ? "NOROT_FIXED" : "NOROT_BC";

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
    for (int d = 0; d < ndom; d++) {
    switch (d) {
    // in the star the constraint equations are sourced by the matter
    case 0:
    case 1:
      // sources
      syst->add_def(d, "Etilde = press * h - press * delta") ;
      syst->add_def(d, "Stilde = 3 * press * delta") ;
 
      // constraint equations
      syst->add_def(d, "eqP    = delta * D^i D_i P + 4piG / 2. * P^5 * Etilde") ;
      syst->add_def(d, "eqNP   = delta * D^i D_i NP - 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
 
      // definition for the baryonic mass integral
      syst->add_def(d, "intMb = P^6 * rho");
      // first integral of the euler equation for a static, non-rotating star, i.e. a TOV
      syst->add_def(d, "firstint = H + log(N)");
 
      break;
    // outside the matter is absent and the sources are zero
    default:
      syst->add_eq_full(d, "H = 0");
 
      syst->add_def(d, "eqP = D^i D_i P");
      syst->add_def(d, "eqNP = D^i D_i NP");
      break;
    }
  }

  }

  template<class eos_t>
  void NS_XCTS_NOROT<eos_t>::syst_init() {
    using namespace ::Kadath::Margherita;
    
    // call the (flat) conformal metric "f"
    fmet->set_system(*syst, "f");
    syst->add_var("H"   , logh);
  
    // define numerical constants
    syst->add_cst("4piG", (*bconfig)(BCO_QPIG));
    
    // include the coordinate fields
    syst->add_cst("einf", *(*coord_vectors)[S_INF]);
    
    // the basic fields, conformal factor, lapse and (log) enthalpy
    syst->add_var("P"   , *conformal_factor);
    syst->add_var("N"   , *lapse);
    
    // define common combinations of conformal factor and lapse
    syst->add_def("NP = P*N");
    syst->add_def("Ntilde = N / P^6");
  
    // define quantity to be integrated at infinity
    // two (in this case) equivalent definitions of ADM mass
    // as well as the Komar mass
    syst->add_def(ndom - 1, "intMadm = - einf^i * D_i P / 4piG * 2");
    syst->add_def(ndom - 1, "intMk = einf^i * D_i N / 4piG");
    syst->add_def(ndom - 1, "intMadmalt = -dr(P) * 2 / 4piG");
    
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
  void NS_XCTS_NOROT<eos_t>::print_diagnostics(const int ite, const double conv) const {

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
    std::cout << FORMAT << "R: " << rs[0] << " " << rs[1] << "\n";
    std::cout.flags(f);
    #undef FORMAT
    std::cout << "=======================================" << "\n\n";
  } // end print diagnostics norot

  template<class eos_t>
  void NS_XCTS_NOROT<eos_t>::update_config_quantities() {

    auto rs = bco_utils::get_rmin_rmax(*space, 1);
    bconfig->set(BCO_PARAMS::RMID) = rs[0];

    // compute the ADM mass as surface integral at infinity
    Val_domain integMadm(syst->give_val_def("intMadm")()(ndom - 1));
    double Madm = space->get_domain(ndom - 1)->integ(integMadm, OUTER_BC);
    bconfig->set(BCO_PARAMS::MADM) = Madm;

    // compute the baryonic mass at volume integral from the given integrant
    double baryonic_mass =
      syst->give_val_def("intMb")()(0).integ_volume() +
      syst->give_val_def("intMb")()(1).integ_volume();
    bconfig->set(BCO_PARAMS::MB) = baryonic_mass;

    auto loghc = bco_utils::get_boundary_val(0, *logh, INNER_BC);

    // compute the ADM Angular Momentum as surface integral at infinity
    double chi{0};

    if(seq) {
      auto idx{seq->mass_idx()};
      switch(idx) {
        case BCO_PARAMS::HC:
          bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
          break;
        case BCO_PARAMS::NC:
          bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
          break;
        default:
          bconfig->set(BCO_PARAMS::HC) = std::exp(loghc);
          bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
          break;
      }
    } else {
      bconfig->set(BCO_PARAMS::NC) = EOS<eos_t,DENSITY>::get(bconfig->set(BCO_PARAMS::HC));
      bconfig->set(BCO_PARAMS::CHI) = chi;
    }
    bconfig->set(BCO_PARAMS::QLMADM) = bconfig->set(BCO_PARAMS::MADM) ;
  }

  template<class eos_t>
  void NS_XCTS_NOROT<eos_t>::save_to_file() const {
    Kadath::bco_utils::save_to_file(*space, *bconfig, *conformal_factor, *lapse, *shift, *logh);
  }

  template<class eos_t>
  void NS_XCTS_NOROT<eos_t>::load_solution_from_file() {
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
  NS_XCTS_NOROT<eos_t>::NS_XCTS_NOROT(NS_XCTS_BASE::base_config_t& config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, int const rank_) :
      NS_XCTS_BASE(config_, seq_, res_, outputdir_, rank_) {
    
    solver_stage = ::Kadath::FUKA_Config::STAGES::NOROT_BC;
    if(!seq->is_set() && !bconfig->control(CONTROLS::SEQUENCES)) {
      initialize_config_from_fixing_values(*bconfig, *seq);
    }
    load_solution_from_file();
    initialize_EOS(*this);
    initialize_support_containers();
    cout << *seq << endl;
  }
    
  // template<class config_t>
  // inline auto initialize_NS_XCTS_NOROT_solver(config_t& bconfig) {
  //   // Initialize EOS
  //   auto h_cut    = bconfig.template eos<double>(Kadath::FUKA_Config::EOS_PARAMS::HCUT);
  //   auto eos_file = bconfig.template eos<std::string>(Kadath::FUKA_Config::EOS_PARAMS::EOSFILE);
  //   auto eos_type = bconfig.template eos<std::string>(Kadath::FUKA_Config::EOS_PARAMS::EOSTYPE);

  //   if(eos_type == "Cold_Table") {
  //     using namespace Kadath::Margherita;
  //     using eos_t = Kadath::Margherita::Cold_Table;

  //     const int interp_pts = (bconfig.template eos<int>(Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS) == 0) ? \
  //                             2000 : (*bconfig.template eos<int>(Kadath::FUKA_Config::EOS_PARAMS::INTERP_PTS);

  //     EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
  //     NS_XCTS_NOROT<eos_t> solver(bconfig);
  //     solver.set_h_cut() = h_cut;
  //     solver.set_eos_file() = eos_file;
  //     solver.set_eos_type() = eos_type;
  //     return solver;
  //   } else if(eos_type == "Cold_PWPoly") {
  //     using namespace Kadath::Margherita;
  //     using eos_t = Kadath::Margherita::Cold_PWPoly;
  //     EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
  //     NS_XCTS_NOROT<eos_t> solver(bconfig);
  //     solver.set_h_cut() = h_cut;
  //     solver.set_eos_file() = eos_file;
  //     solver.set_eos_type() = eos_type;
  //     return solver;
  //   } 
  //   throw std::invalid_argument("\nInvalid EOS Type\n)");
  // }
}