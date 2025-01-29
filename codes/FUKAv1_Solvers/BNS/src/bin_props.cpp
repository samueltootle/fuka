#include "Configurator/config_binary.hpp"
#include "EOS/EOS.hh"
#include "coord_fields.hpp"
#include "kadath_bin_ns.hpp"
#include "python_reader.hpp"

using namespace Kadath::Margherita;

template <eos_var_t var>
using eos = EOS<Cold_Table, var>;

// the space type
typedef Kadath::Space_bin_ns space_t;

// specialized quantities for a BNS system
struct bns_vars_t : public Kadath::vars_base_t<bns_vars_t> {};

// define the actual quantities and their order in the file!
template <>
Kadath::var_vector Kadath::vars_base_t<bns_vars_t>::vars = {
    {"conf", SCALAR}, {"lapse", SCALAR}, {"shift", VECTOR},
    {"logh", SCALAR}, {"phi", SCALAR},
};

class bns_reader_t : public Kadath::python_reader_t<space_t, bns_vars_t> {
  std::string config_filename;
  kadath_config_boost<BIN_INFO> bconfig;

 public:
  bns_reader_t(std::string const filename)
      : Kadath::python_reader_t<space_t, bns_vars_t>(filename),
        config_filename(filename.substr(0, filename.size() - 3) + "info"),
        bconfig(config_filename) {
    this->compute_defs();
  }

  void compute_defs() {
    Kadath::Scalar const& conf = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const& lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const& shift = extractField<Kadath::Vector>("shift");
    Kadath::Scalar const& logh = extractField<Kadath::Scalar>("logh");
    Kadath::Scalar const& phi = extractField<Kadath::Scalar>("phi");

    double& Mb1 = bconfig(MB, BCO1);
    double& Mb2 = bconfig(MB, BCO2);
    double Mbtot = Mb1 + Mb2;

    double& Madm1 = bconfig(MADM, BCO1);
    double& Madm2 = bconfig(MADM, BCO2);
    double Minf = Madm1 + Madm2;

    double& q = bconfig(Q);

    double& chi1 = bconfig(CHI, BCO1);
    double& chi2 = bconfig(CHI, BCO2);

    Base_tensor basis(shift.get_basis());
    int ndom = space.get_nbr_domains();

    Index center_pos(space.get_domain(space.NS1)->get_nbr_points());
    double xc1 = space.get_domain(space.NS1)->get_cart(1)(center_pos);
    double xc2 = space.get_domain(space.NS2)->get_cart(1)(center_pos);
    double xo = space.get_domain(ndom - 1)->get_cart(1)(center_pos);

    CoordFields<Space_bin_ns> cfields(space);

    std::array<Vector, NUM_VECTORS> coord_vectors{
        Vector(space, CON, basis), Vector(space, CON, basis),
        Vector(space, CON, basis), Vector(space, CON, basis),
        Vector(space, CON, basis), Vector(space, CON, basis),
        Vector(space, CON, basis), Vector(space, CON, basis),
        Vector(space, CON, basis)};
    std::array<Scalar, NUM_SCALARS> coord_scalars{Scalar(space), Scalar(space)};
    update_fields(cfields, coord_vectors, coord_scalars, xo, xc1, xc2);

    Metric_flat fmet(space, basis);

    {
      System_of_eqs syst(space, 0, ndom - 1);

      fmet.set_system(syst, "f");

      if (bconfig.template eos<std::string>(EOSTYPE, BCO1) == "Cold_PWPoly") {
        using eos_t = Kadath::Margherita::Cold_PWPoly;

        EOS<eos_t, PRESSURE>::init();

        Param p;
        syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
        syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
        syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
        syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

      } else if (bconfig.template eos<std::string>(EOSTYPE, BCO1) ==
                 "Cold_Table") {
        using eos_t = Kadath::Margherita::Cold_Table;

        EOS<eos_t, PRESSURE>::init(
            bconfig.template eos<std::string>(EOSFILE, BCO1));

        Param p;
        syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
        syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
        syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
        syst.add_ope("dHdlnrho", &EOS<eos_t, DHDRHO>::action, &p);

      } else {
        abort();
      }

      Index p1(space.get_domain(0)->get_nbr_points());
      double loghc1 = logh(0)(p1);
      double rhoc1 = eos<DENSITY>::get(std::exp(loghc1));
      double pressc1 = eos<PRESSURE>::get(std::exp(loghc1));

      Index p2(space.get_domain(3)->get_nbr_points());
      double loghc2 = logh(3)(p2);
      double rhoc2 = eos<DENSITY>::get(std::exp(loghc2));
      double pressc2 = eos<PRESSURE>::get(std::exp(loghc2));

      syst.add_cst("4piG", bconfig(QPIG));
      syst.add_cst("PI", M_PI);

      syst.add_cst("Mb1", bconfig(MB, BCO1));
      syst.add_cst("chi1", bconfig(CHI, BCO1));

      syst.add_cst("Mb2", bconfig(MB, BCO2));
      syst.add_cst("chi2", bconfig(CHI, BCO2));

      syst.add_cst("mg", coord_vectors[GLOBAL_ROT]);
      syst.add_cst("mm", coord_vectors[BCO1_ROT]);
      syst.add_cst("mp", coord_vectors[BCO2_ROT]);

      syst.add_cst("ex", coord_vectors[EX]);
      syst.add_cst("ey", coord_vectors[EY]);
      syst.add_cst("ez", coord_vectors[EZ]);

      syst.add_cst("sm", coord_vectors[S_BCO1]);
      syst.add_cst("sp", coord_vectors[S_BCO2]);
      syst.add_cst("einf", coord_vectors[S_INF]);

      syst.add_cst("rm", coord_scalars[R_BCO1]);
      syst.add_cst("rp", coord_scalars[R_BCO2]);

      syst.add_cst("Madm1", bconfig(MADM, BCO1));
      syst.add_cst("Madm2", bconfig(MADM, BCO2));

      syst.add_cst("qlMadm1", bconfig(QLMADM, BCO1));
      syst.add_cst("qlMadm2", bconfig(QLMADM, BCO2));

      syst.add_cst("Hc1", loghc1);
      syst.add_cst("Hc2", loghc2);

      syst.add_cst("omes1", bconfig(OMEGA, BCO1));
      syst.add_cst("omes2", bconfig(OMEGA, BCO2));

      syst.add_cst("xaxis", bconfig(COM));
      syst.add_cst("ome", bconfig(GOMEGA));

      syst.add_cst("P", conf);
      syst.add_cst("N", lapse);
      syst.add_cst("bet", shift);
      syst.add_cst("H", logh);
      syst.add_cst("phi", phi);

      syst.add_def("h      = exp(H)");

      syst.add_def("press = press(h)");
      syst.add_def("eps = eps(h)");
      syst.add_def("rho = rho(h)");
      syst.add_def("dHdlnrho = dHdlnrho(h)");

      syst.add_def("delta = h - eps - 1.");

      syst.add_def("NP     = P*N");
      syst.add_def("Ntilde = N / P^6");

      syst.add_def("Morb^i = mg^i + xaxis * ey^i");
      syst.add_def("omega^i= bet^i + ome * Morb^i");

      for (int d = space.NS1; d <= space.ADAPTED1; ++d) {
        syst.add_def(d, "s^i  = omes1 * mm^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }
      for (int d = space.NS2; d <= space.ADAPTED2; ++d) {
        syst.add_def(d, "s^i  = omes2 * mp^i");
        syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
      }

      syst.add_def(
          "A^ij   = (D^i bet^j + D^j bet^i - 2. / 3.* D_k bet^k * f^ij) / 2. / "
          "Ntilde");

      syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i");
      syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i");
      syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i");

      syst.add_def(space.ADAPTED1 + 1,
                   "intS1 = A_ij * mm^i * sm^j / 2. / 4piG");
      syst.add_def(space.ADAPTED2 + 1,
                   "intS2 = A_ij * mp^i * sp^j / 2. / 4piG");

      for (int d = 0; d < ndom; d++) {
        if ((d >= space.ADAPTED2 + 1) || d == space.ADAPTED1 + 1) {
          syst.add_def(d, "eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8");
          syst.add_def(
              d, "eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
          syst.add_def(d,
                       "eqbet^i = D_j D^j bet^i + D^i D_j bet^j / 3. - 2. * "
                       "A^ij * D_j Ntilde");
        } else {
          syst.add_def(d, "Ucor^i    = omega^i / N");
          syst.add_def(d, "Ucorsquare = P^4 * Ucor_i * Ucor^i");
          syst.add_def(d, "Wcorsquare = 1 / (1 - Ucorsquare)");
          syst.add_def(d, "Wcor = sqrt(Wcorsquare)");

          syst.add_def(d, "Wsquare= eta^i * eta_i / h^2 / P^4 + 1.");
          syst.add_def(d, "W      = sqrt(Wsquare)");

          syst.add_def(d, "U^i    = eta^i / P^4 / h / W");
          syst.add_def(d, "V^i    = N * U^i - omega^i");
          syst.add_def(d, "Usquare= P^4 * U_i * U^i");

          syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)");

          syst.add_def(d,
                       "eqphi  = P^6 * W * V^i * D_i H + dHdlnrho * D_i (P^6 * "
                       "W * V^i)");

          syst.add_def(d, "Etilde = press * h * Wsquare - press * delta");
          syst.add_def(d,
                       "Stilde = 3 * press * delta + (Etilde + press * delta) "
                       "* Usquare");
          syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i");

          syst.add_def(d,
                       "eqP    = delta * D^i D_i P + A_ij * A^ij / P^7 / 8 * "
                       "delta + 4piG / 2. * P^5 * Etilde");
          syst.add_def(d,
                       "eqNP   = delta * D^i D_i NP - 7. / 8. * NP / P^8 * "
                       "delta * A_ij *A^ij "
                       "- 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
          syst.add_def(
              d,
              "eqbet^i= delta * D_j D^j bet^i + delta * D^i D_j bet^j / 3. "
              "- 2. * delta * A^ij * D_j Ntilde - 4. * 4piG * N * P^4 * "
              "ptilde^i");

          syst.add_def(d, "intMb  = P^6 * rho * W");
          syst.add_def(d, "intM   = - D_i D^i P * 2. / 4piG");
        }
      }

      syst.add_def(ndom - 1, "Madm = -dr(P) * 2 / 4piG");

      vars["dist"] = bconfig(DIST);
      vars["q"] = q;

      vars["chi1"] = chi1;
      vars["chi2"] = chi2;

      vars["omega"] = bconfig(GOMEGA);
      vars["mOmega"] = Minf * bconfig(GOMEGA);

      vars["Madm"] = space.get_domain(ndom - 1)->integ(syst.give_val_def(
                                                           "Madm")()(ndom - 1),
                                                       OUTER_BC);
      vars["Minf"] = Minf;

      vars["Madm1"] = Madm1;
      vars["Madm2"] = Madm2;

      vars["Eb"] = vars["Madm"] - vars["Minf"];
      vars["EboM"] = vars["Eb"] / vars["Minf"];
    }
  }
};

BOOST_PYTHON_MODULE(bin_props) {
  // initialize python types
  Kadath::initPythonBinding<space_t>();
  Kadath::constructPythonReader<bns_reader_t>("bns_reader");
}
