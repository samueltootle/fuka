#include "Solvers/fuka_syst/fuka_syst_setup.hpp"
#include "utilities.hpp"
namespace Kadath::FUKA_Solvers {
  // NOROT Routines
  template<class eos_t>
  NS_XCTS_DIFF_ROT<eos_t>::NS_XCTS_DIFF_ROT(NS_XCTS_BASE::base_config_t* config_, ns_sequence const & seq_, 
    Parameter_sequence<BCO_PARAMS> const & res_, std::string outputdir_, int const rank_) :
      NS_XCTS_BASE(config_, seq_, res_, outputdir_, rank_), spinup(nullptr) {
    
    stagename = "DIFF_ROT";
    solver_stage = ::Kadath::FUKA_Config::STAGES::DIFF_ROT;
    if(!seq->is_set() && !bconfig->control(CONTROLS::SEQUENCES)) {
      initialize_config_from_fixing_values(*bconfig, *seq);
    }
    load_solution_from_file();
    initialize_EOS(*this);
    initialize_support_containers();
    initialize_diffrot_params();
    initialize_spinup();

    if(rank == 0) {
      cout << *seq << endl;
      cout << *resolution << endl;
    }
  }

  template<class eos_t>
  std::string NS_XCTS_DIFF_ROT<eos_t>::converged_filename(const std::string stage) const {
    // FIXME assumes a fix resolution for all domains
    auto res = space->get_domain(0)->get_nbr_points()(0);
    const std::string eosname{extract_eos_name(*bconfig)};
    std::stringstream ss;
    ss << "NS";
    if(stage != "") ss  << "_" << stage << ".";
    else ss << ".";
    ss << eosname << "."
       << law << ".";
    
    // Add mass fixing parameter to filename
    auto default_idx = BCO_PARAMS::MADM;
    if(seq) {
      update_filename_from_mass_fixing(*bconfig, seq, ss);
    } else {
      auto [ seq_key, tidx ] = get_key_val_pair_from_val(MBCO_PARAMS, default_idx);
      ss << seq_key << "." << (*bconfig)(default_idx) << "."; 
    }  
    
    if(law == "keh") {
      ss << "Ar." << (*bconfig).template diffrot<double>(DIFFROT_PARAMS::DIFF_ARATIO) << "."
         << "Rr." << (*bconfig).template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO) << ".";
    }
    ss << (*bconfig)(BCO_PARAMS::NSHELLS) << "."
       <<std::setfill('0') << std::setw(2) << res;
    return ss.str();
  }

  template<class eos_t>
  void NS_XCTS_DIFF_ROT<eos_t>::setup_syst() {
    int exit_status = EXIT_SUCCESS;

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
    
    // Update vector fields
    update_fields_co(*cfields, *coord_vectors,{}, 0.);
    initialize_diffrot_params();

    // Initialize KEH parameters
    auto adpt_dom = space->get_domain(1);
    double R0 = adpt_dom->get_radius()(*pos_eq);
    double Rp = adpt_dom->get_radius()(*pos_pole);
    auto& diffAratio = diffrot_params[DIFFROT_PARAMS::DIFF_ARATIO];
    auto& diffRratio = diffrot_params[DIFFROT_PARAMS::DIFF_RRATIO];
    diffA = diffrot_params[DIFFROT_PARAMS::DIFF_ARATIO] * R0;
    bconfig->set(BCO_PARAMS::RMID) = R0;
    std::string firstint{"firstint = (H + log(N) - log(W)) - 0.5 * j^2 / diffA^2"};

    if (rank == 0) {
      std::cout << "###################################" << std::endl
                << "Differential Rotating models (KEH)"  << std::endl
                << firstint << std::endl
                << "Fixed A / R0: " << diffAratio << std::endl
                << "Fixed Rp / Re: " << diffRratio << std::endl
                << "Initial Rp/Re: " << Rp / R0 << "\n"
                << "Initial R0: " << R0 <<std::endl
                << "###################################" << "\n\n";
    }
    // Setup System of equations
    syst.reset(new System_of_eqs(*space));
    syst_init();

    // Setup mass fixing parameter
    std::string central_fixing_definition{"h - hc"};    
    if(seq) {
      central_fixing_definition = ::Kadath::FUKA_Syst_tools::set_ns_mass_fixing(*syst, *bconfig, seq);
    } else {
      syst->add_cst("hc", (*bconfig)(BCO_PARAMS::HC));
    }
    
    // Utility field
    syst->add_cst("one", *ones);
    
    // KEH Constants
    syst->add_cst("diffAratio", diffrot_params[DIFFROT_PARAMS::DIFF_ARATIO]);
    syst->add_cst("Rratio"    , diffrot_params[DIFFROT_PARAMS::DIFF_RRATIO]);

    // KEH vars
    syst->add_var("diffA", diffA);
    syst->add_var("omec", (*bconfig)(BCO_PARAMS::OMEGA));
    syst->add_var("R0"  , (*bconfig)(BCO_PARAMS::RMID));
    syst->add_var("Omega", *diff_omega);

    // KEH Definitions
    syst->add_def("diffAField = one * diffA");
    syst->add_def("r = multr(one)");

    syst->add_def("omega^i = bet^i + Omega * mg^i");
    syst->add_def("U^i = omega^i / N");
    syst->add_def("Usquare = P^4 * U_i * U^i");
    syst->add_def("Wsquare = 1. / (1. - Usquare)");
    syst->add_def("W = sqrt(Wsquare)");
    syst->add_def("j = P^4 * Wsquare * f_ij * U^i * mg^j / N");
    syst->add_def("omelaw = omec - j / diffA^2");
    
    for (int d = 0; d < ndom; d++) {
      switch (d) {
      case 0:
      case 1:
        syst->add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
        syst->add_def(d, "Stilde = 3 * press * delta + (Etilde + press * delta) * Usquare") ;
        syst->add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

        syst->add_def(d, "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * delta + 4piG / 2. * P^5 * Etilde") ;
        syst->add_def(d, "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * delta * A_ij *A^ij "
                              "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
        syst->add_def(d, "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
                              "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * ptilde^i");

        syst->add_def(d, "intMb = P^6 * rho(h) * W");
        syst->add_def(d, firstint.c_str());
        syst->add_eq_full(d, "Omega - omelaw = 0");
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
    // Ensure Omega field matches the interior solution, but is zero otherwise.
    // Much more reliable convergence
    syst->add_eq_matching(2, INNER_BC, "Omega");
    syst->add_eq_matching(2, INNER_BC, "dn(Omega)");
    syst->add_eq_inside(2,"Omega = 0");
    syst->add_eq_full(3, "Omega = 0");

    // Enforce Differential rotation parameters at the origin,
    // equator, and pole.
    syst->add_eq_val(0, "diffAField/R0 - diffAratio", *pos_origin);
    syst->add_eq_val(1, "r/R0 - 1", *pos_eq);
    syst->add_eq_val(1, "r/R0 - Rratio", *pos_pole);

    // add the constraint equations and demand continuity 
    // of their normal derivative across domain boundaries
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
    }
  }

  template<class eos_t>
  void NS_XCTS_DIFF_ROT<eos_t>::syst_init() {
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
  void NS_XCTS_DIFF_ROT<eos_t>::print_diagnostics(const int ite, const double conv) const {

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
              << std::abs(Madm - Mk) / Madm << "]" << std::endl
              << FORMAT << "R: " << rs[0] << " " << rs[1] 
                        << " [" << rs[0] / rs[1] << "]" << std::endl;
    std::cout << FORMAT << "Jadm: " << J << std::endl
              << FORMAT << "Chi: " << J / Madm / Madm << " [" << (*bconfig)(CHI) << "]\n"
              << FORMAT << "Omega: " << (*bconfig)(OMEGA) << std::endl;
    std::cout.flags(f);
    #undef FORMAT
    std::cout << "=======================================" << "\n\n";
  } // end print diagnostics norot

  template<class eos_t>
  void NS_XCTS_DIFF_ROT<eos_t>::update_config_quantities() {

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
  void NS_XCTS_DIFF_ROT<eos_t>::initialize_spinup() {
    
    auto adpt_dom = space->get_domain(1);
    double const R0 = adpt_dom->get_radius()(*pos_eq);
    double const Rp = adpt_dom->get_radius()(*pos_pole);
    axis_ratio = Rp / R0;
    double const diffRratio = diffrot_params[DIFFROT_PARAMS::DIFF_RRATIO];

    // no reason to spinup if the solution is already sufficiently rotating
    if( std::fabs(1. - axis_ratio / diffRratio) < 0.1){
      return;
    }
    spinup.reset(new Parameter_sequence<DIFFROT_PARAMS>("R_ratio",DIFFROT_PARAMS::DIFF_RRATIO));
    spinup->set(axis_ratio, axis_ratio, diffRratio);
    auto const spinidx = seq->spin_idx();
    const double dx = 0.1;
    const int N = int((axis_ratio - diffRratio) / dx);
    spinup->set_N(N);
    bconfig->set(spinidx) = 1.0;
    diffrot_params[spinidx] = 1.0;
  }

  template<class eos_t>
  bool NS_XCTS_DIFF_ROT<eos_t>::increment_spin() {
    if(!spinup || !spinup->is_set())
      return false;
    auto sequence_var_indices = spinup->get_indices();
    auto const & dx = spinup->step_size();
    double const diffRratio = (*bconfig).template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO);
    auto x = diffRratio + dx;
    if(spinup->loop_condition(x)) {
      x = (x > spinup->final()) ? spinup->final() : x;
      bconfig->set_diffrot(sequence_var_indices) = x;
      return true;
    }
    return false;
  }

  template<class eos_t>
  void NS_XCTS_DIFF_ROT<eos_t>::initialize_diffrot_params() {
    law = [&]() -> std::string {
      auto v = (*bconfig).template diffrot<std::string>(DIFFROT_PARAMS::DIFF_LAW);
      return str_tolower(v);
    }();
    if(law == "keh") {
      diffrot_params.fill(false);
      diffrot_params[DIFFROT_PARAMS::DIFF_ARATIO] = (*bconfig).template diffrot<double>(DIFFROT_PARAMS::DIFF_ARATIO);
      diffrot_params[DIFFROT_PARAMS::DIFF_RRATIO] = (*bconfig).template diffrot<double>(DIFFROT_PARAMS::DIFF_RRATIO);
    }
    pos_origin.reset(new Index(space->get_domain(0)->get_nbr_points()));
    auto npts = space->get_domain(1)->get_nbr_points();
    pos_eq.reset(new Index(npts));
    pos_eq->set(0) = npts(0) - 1; /// Set to outer radius
    pos_eq->set(1) = npts(1) - 1; /// Set theta to be on the xy plane.

    pos_pole.reset(new Index(npts));
    pos_pole->set(0) = npts(0) - 1; /// Set to outer radius
  
    ones.reset(new Scalar(*space));
    *ones = 1.;
    ones->std_base();
  }

  template<class eos_t>
  void NS_XCTS_DIFF_ROT<eos_t>::reset_all_ptrs() {
    // Fields
    conformal_factor.reset(nullptr) ;
    lapse.reset(nullptr);
    shift.reset(nullptr);
    logh.reset(nullptr);
    diff_omega.reset(nullptr);

    // Containers
    basis.reset(nullptr);
    fmet.reset(nullptr);
    syst.reset(nullptr);
    cfields.reset(nullptr);
    coord_vectors.reset(nullptr);
    ones.reset(nullptr);
    pos_origin.reset(nullptr);
    pos_eq.reset(nullptr);
    pos_pole.reset(nullptr);

    // Space
    space.reset(nullptr);
  }
}
