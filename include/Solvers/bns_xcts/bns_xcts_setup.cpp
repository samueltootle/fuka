/**
 * \addtogroup BNS_XCTS
 * \ingroup FUKA
 * @{*/
#include "EOS/FUKA_EOS_Utilities.hh"
namespace Kadath {
namespace FUKA_Solvers {
template <class config_t>
inline void bns_xcts_setup_headon_config(config_t& bconfig) {
  using namespace ::Kadath::bco_utils;
  check_dist(bconfig(BIN_PARAMS::DIST), bconfig(BCO_PARAMS::MADM, NODES::BCO1),
             bconfig(BCO_PARAMS::MADM, NODES::BCO2));

  // Binary Parameters
  bconfig.set(BIN_PARAMS::REXT) = 2 * bconfig(BIN_PARAMS::DIST);

  bconfig.set(BIN_PARAMS::Q) = bconfig(BCO_PARAMS::MADM, NODES::BCO2) /
                               bconfig(BCO_PARAMS::MADM, NODES::BCO1);

  // classical Newtonian estimate
  bconfig.set(BIN_PARAMS::COM) = com_estimate(
      bconfig(BIN_PARAMS::DIST), bconfig(BCO_PARAMS::MADM, NODES::BCO1),
      bconfig(BCO_PARAMS::MADM, NODES::BCO2));
}

template <class config_t>
inline void bns_xcts_setup_bin_config(config_t& bconfig) {
  using namespace ::Kadath::bco_utils;
  bns_xcts_setup_headon_config(bconfig);

  // obtain 3PN estimate for the global, orbital omega
  KadathPNOrbitalParams(bconfig, bconfig(BCO_PARAMS::MADM, NODES::BCO1),
                        bconfig(BCO_PARAMS::MADM, NODES::BCO2));

  // delete ADOT, this can always be recalculated during
  // the eccentricity reduction stage
  bconfig.reset(BIN_PARAMS::ADOT);
}

template <class config_t>
void bns_xcts_setup_space(config_t& bconfig) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<int, 2> bcos{NODES::BCO1, NODES::BCO2};
  std::array<std::string, 2> filenames;

  for (int i = 0; i < 2; ++i)
    filenames[i] = solve_NS_from_binary(bconfig, bcos[i]);

  // debugging only
  for (auto& f : filenames)
    if (rank == 0)
      std::cout << f << std::endl;

  if (rank == 0)
    bns_xcts_superimposed_import(bconfig, filenames);
  MPI_Barrier(MPI_COMM_WORLD);
  bconfig.control(CONTROLS::SEQUENCES) = false;
  bconfig.open_config();
}

template <typename eos_t>
struct bns_setup_boosted_3d {
inline void operator()(kadath_config_boost<BCO_NS_INFO>& NS1config,
                                 kadath_config_boost<BCO_NS_INFO>& NS2config,
                                 kadath_config_boost<BIN_INFO>& bconfig) {
  using namespace ::Kadath::bco_utils;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // read single NS configuration and data
  std::string nsspacein = NS1config.space_filename();

  FILE* ff1 = fopen(nsspacein.c_str(), "r");
  Space_spheric_adapted spacein1(ff1);
  Scalar confin1(spacein1, ff1);
  Scalar lapsein1(spacein1, ff1);
  Vector shiftin1(spacein1, ff1);
  Scalar loghin1(spacein1, ff1);
  Scalar phiin1(spacein1, ff1);
  fclose(ff1);
  int ndomin1 = spacein1.get_nbr_domains();

  // update NS1config quantities before updating binary configuration file
  NS1config.set(BCO_PARAMS::HC) =
      std::exp(get_boundary_val(0, loghin1, INNER_BC));
  NS1config.set(BCO_PARAMS::NC) =
      EOS<eos_t, DENSITY>::get(NS1config(BCO_PARAMS::HC));
  update_config_NS_radii(spacein1, NS1config, 1);

  // read single NS configuration and data
  nsspacein = NS2config.space_filename();

  FILE* ff2 = fopen(nsspacein.c_str(), "r");
  Space_spheric_adapted spacein2(ff1);
  Scalar confin2(spacein2, ff2);
  Scalar lapsein2(spacein2, ff2);
  Vector shiftin2(spacein2, ff2);
  Scalar loghin2(spacein2, ff2);
  Scalar phiin2(spacein2, ff2);
  fclose(ff2);
  int ndomin2 = spacein2.get_nbr_domains();

  // update NS2config quantities before updating binary configuration file
  NS2config.set(BCO_PARAMS::HC) =
      std::exp(get_boundary_val(0, loghin2, INNER_BC));
  NS2config.set(BCO_PARAMS::NC) =
      EOS<eos_t, DENSITY>::get(NS2config(BCO_PARAMS::HC));
  update_config_NS_radii(spacein2, NS2config, 1);

  // update NS parameters in binary config
  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    bconfig.set(i, NODES::BCO1) = NS1config.set(i);

  for (int i = 0; i < BCO_PARAMS::NUM_BCO_PARAMS; ++i)
    bconfig.set(i, NODES::BCO2) = NS2config.set(i);

  auto gen_radius_field = [&](auto& spacein, auto& old_space_radius,
                              const int ndomin) {
    // get the adapted domain of the single star
    const Domain_shell_outer_adapted* old_outer_adapted =
        dynamic_cast<const Domain_shell_outer_adapted*>(spacein.get_domain(1));
    for (int i = 0; i < ndomin; ++i)
      if (i != 1)
        old_space_radius.set_domain(i) = spacein.get_domain(i)->get_radius();
    old_space_radius.set_domain(1) = old_outer_adapted->get_outer_radius();
    old_space_radius.std_base();
  };

  auto interp_field = [&](auto& space, int outer_dom, auto& old_phi) {
    const int d = outer_dom;
    const Domain_shell_inner_adapted* old_inner =
        dynamic_cast<const Domain_shell_inner_adapted*>(
            space.get_domain(d + 1));
    // interpolate old_phi field outside of the star for import
    update_adapted_field(old_phi, d, d + 1, old_inner, INNER_BC);
  };
  // in case we used boosted TOVs, we need to import PHI
  if (NS1config.set_field(BCO_FIELDS::PHI) == true)
    interp_field(spacein1, 1, phiin1);
  if (NS2config.set_field(BCO_FIELDS::PHI) == true)
    interp_field(spacein2, 1, phiin2);

  // create scalar field representing the radius in the single star space
  Scalar old_space_radius1(spacein1);
  old_space_radius1.annule_hard();
  gen_radius_field(spacein1, old_space_radius1, ndomin1);

  // create scalar field representing the radius in the single star space
  Scalar old_space_radius2(spacein2);
  old_space_radius2.annule_hard();
  gen_radius_field(spacein2, old_space_radius2, ndomin2);

  // start Update config vars
  double r_max_tot = std::max(bconfig(BCO_PARAMS::RMID, NODES::BCO1),
                              bconfig(BCO_PARAMS::RMID, NODES::BCO2));
  const double rout_sep_est =
      (bconfig(BIN_PARAMS::DIST) / 2. - r_max_tot) / 3. + r_max_tot;
  bconfig.set(BCO_PARAMS::ROUT, NODES::BCO1) = rout_sep_est;
  bconfig.set(BCO_PARAMS::ROUT, NODES::BCO2) =
      bconfig(BCO_PARAMS::ROUT, NODES::BCO1);
  // end updating config vars

  // setup domain boundaries
  std::vector<double> out_bounds(1 + bconfig(BIN_PARAMS::OUTER_SHELLS));

  // scale outer shells by constant steps of 1/4 for the time being
  for (int e = 0; e < out_bounds.size(); ++e)
    out_bounds[e] = bconfig(BIN_PARAMS::REXT) * (1. + e * 0.25);

  // Determine optimal grid struction based on isolated solutions
  std::vector<double> NS1_bounds;
  {
    auto drPsi(compute_drPsi(spacein1, confin1,
                             Metric_flat(spacein1, shiftin1.get_basis()),
                             {0, 1}));
    NS1_bounds = bco_u::set_arb_boundsv3(bconfig, drPsi, 2, NODES::BCO1);
  }
  std::vector<double> NS2_bounds;
  {
    auto drPsi(compute_drPsi(spacein2, confin2,
                             Metric_flat(spacein2, shiftin2.get_basis()),
                             {0, 1}));
    NS2_bounds = bco_u::set_arb_boundsv3(bconfig, drPsi, 2, NODES::BCO2);
  }
  // end setup domain boundaries

  // output stellar domain radii
  if (rank == 0) {
    std::cout << "Bounds:" << std::endl;
    print_bounds("NS1-bounds", NS1_bounds);
    print_bounds("NS2-bounds", NS2_bounds);
    print_bounds("Outer", out_bounds);
  }

  // create a binary neutron star space
  int typer = CHEB_TYPE;
  Space_bin_ns space(typer, bconfig(BIN_PARAMS::DIST), NS1_bounds, NS2_bounds,
                     out_bounds, bconfig(BIN_PARAMS::BIN_RES));
  // with cartesian type basis
  Base_tensor basis(space, CARTESIAN_BASIS);

  // get the domains with inner respectively outer adapted boundary
  // for each of the stars
  std::array<const Domain_shell_outer_adapted*, 2> new_outer_adapted{
      dynamic_cast<const Domain_shell_outer_adapted*>(
          space.get_domain(space.ADAPTED1)),
      dynamic_cast<const Domain_shell_outer_adapted*>(
          space.get_domain(space.ADAPTED2))};

  std::array<const Domain_shell_inner_adapted*, 2> new_inner_adapted{
      dynamic_cast<const Domain_shell_inner_adapted*>(
          space.get_domain(space.ADAPTED1 + 1)),
      dynamic_cast<const Domain_shell_inner_adapted*>(
          space.get_domain(space.ADAPTED2 + 1))};

  // updated the radius mapping for both NS
  interp_adapted_mapping(new_inner_adapted[0], 1, old_space_radius1);
  interp_adapted_mapping(new_outer_adapted[0], 1, old_space_radius1);

  interp_adapted_mapping(new_inner_adapted[1], 1, old_space_radius2);
  interp_adapted_mapping(new_outer_adapted[1], 1, old_space_radius2);

  // get and print center of each star
  double xc1 = get_center(space, space.NS1);
  double xc2 = get_center(space, space.NS2);

  if (rank == 0)
    std::cout << "xc1: " << xc1 << std::endl << "xc2: " << xc2 << std::endl;

  // start to create the new fields
  // initialized to zero globally
  Scalar logh(space);
  logh.annule_hard();

  Scalar conf(space);
  conf.annule_hard();

  Scalar lapse(space);
  lapse.annule_hard();

  Vector shift(space, CON, basis);
  for (int i = 1; i <= 3; i++)
    shift.set(i).annule_hard();

  Scalar phi(space);
  phi.annule_hard();
  // end create new fields

  const double ns1_invw4 = set_decay(bconfig, NODES::BCO1);
  const double ns2_invw4 = set_decay(bconfig, NODES::BCO2);

  if (rank == 0)
    std::cout << "WeightNS1: " << bconfig(BCO_PARAMS::DECAY, NODES::BCO1)
              << ", "
              << "WeightNS2: " << bconfig(BCO_PARAMS::DECAY, NODES::BCO2)
              << std::endl;

  // start importing the fields from the single star
  int ndom = space.get_nbr_domains();
  for (int dom = 0; dom < ndom; dom++) {
    // get an index in each domain to iterate over all colocation points
    Index new_pos(space.get_domain(dom)->get_nbr_points());

    do {
      // get cartesian coordinates of the current colocation point
      double x = space.get_domain(dom)->get_cart(1)(new_pos);
      double y = space.get_domain(dom)->get_cart(2)(new_pos);
      double z = space.get_domain(dom)->get_cart(3)(new_pos);

      // define a point shifted suitably to the stellar centers in the binary
      Point absol1(3);
      absol1.set(1) = (x - xc1);
      absol1.set(2) = y;
      absol1.set(3) = z;
      double r2 = y * y + z * z;
      double r2_1 = (x - xc1) * (x - xc1) + r2;
      double r4_1 = r2_1 * r2_1;
      double r4_invw4_1 = r4_1 * ns1_invw4;
      double decay_1 = std::exp(-r4_invw4_1);

      Point absol2(3);
      absol2.set(1) = (x - xc2);
      absol2.set(2) = y;
      absol2.set(3) = z;
      double r2_2 = (x - xc2) * (x - xc2) + r2;
      double r4_2 = r2_2 * r2_2;
      double r4_invw4_2 = r4_2 * ns2_invw4;
      double decay_2 = std::exp(-r4_invw4_2);

      if (dom < ndom - 1) {
        conf.set_domain(dom).set(new_pos) =
            1. + decay_1 * (confin1.val_point(absol1) - 1.) +
            decay_2 * (confin2.val_point(absol2) - 1.);

        lapse.set_domain(dom).set(new_pos) =
            1. + decay_1 * (lapsein1.val_point(absol1) - 1.) +
            decay_2 * (lapsein2.val_point(absol2) - 1.);

        logh.set_domain(dom).set(new_pos) =
            0. + decay_1 * loghin1.val_point(absol1) +
            decay_2 * loghin2.val_point(absol2);

        phi.set_domain(dom).set(new_pos) = 0;
        if (NS1config.set_field(BCO_FIELDS::PHI) == true)
          phi.set_domain(dom).set(new_pos) +=
              decay_1 * phiin1.val_point(absol1);
        if (NS2config.set_field(BCO_FIELDS::PHI) == true)
          phi.set_domain(dom).set(new_pos) +=
              decay_2 * phiin2.val_point(absol2);

        for (int i = 1; i <= 3; i++)
          shift.set(i).set_domain(dom).set(new_pos) =
              0. + decay_1 * shiftin1(i).val_point(absol1) +
              decay_2 * shiftin2(i).val_point(absol2);

      } else {
        // We have to set the compactified domain manually
        // since the outer collocation point is always
        // at inf which is undefined numerically
        conf.set_domain(dom).set(new_pos) = 1.;
        lapse.set_domain(dom).set(new_pos) = 1.;
      }
      // loop over all colocation points
    } while (new_pos.inc());
  }  // end importing fields

  // set logh and phi to zero outside of stars explicitly
  // logh.set_domain(space.ADAPTED1+1).annule_hard();
  for (auto i : {space.ADAPTED1 + 1, space.ADAPTED2 + 1}) {
    phi.set_domain(i).annule_hard();
    logh.set_domain(i).annule_hard();
  }
  for (int i = space.OUTER; i < ndom; ++i) {
    phi.set_domain(i).annule_hard();
    logh.set_domain(i).annule_hard();
  }

  // employ standard spectral expansion, compatible with the given paraties
  conf.std_base();
  lapse.std_base();
  logh.std_base();
  shift.std_base();
  phi.std_base();

  // save everything to a binary file
  save_to_file(space, bconfig, conf, lapse, shift, logh, phi);
}
};

template <class config_t>
void bns_xcts_superimposed_import(config_t& bconfig,
                                  std::array<std::string, 2> NSfilenames) {
  using namespace Kadath::FUKA_EOS;

  // load single NS configuration
  std::string ns1filename{NSfilenames[0]};
  kadath_config_boost<BCO_NS_INFO> NS1config(ns1filename);

  std::string ns2filename{NSfilenames[1]};
  kadath_config_boost<BCO_NS_INFO> NS2config(ns2filename);

  // setup eos and update central density
  const std::string eos_type = NS1config.eos<std::string>(EOS_PARAMS::EOSTYPE);
  EOS_Function_Dispatcher::dispatch<bns_setup_boosted_3d>(bconfig, eos_type, NS1config, NS2config,
                                                          bconfig);
}
/** @}*/
}  // namespace FUKA_Solvers
}  // namespace Kadath