/*
 * Copyright 2021
 * This file is part of the KADATH library and published under
 * https://arxiv.org/abs/2103.09911
 *
 * Author: 
 * Samuel D. Tootle <tootle@itp.uni-frankfurt.de>
 * L. Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "kadath_adapted_polar.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "python_reader.hpp"
#include "include/fuka_py.hpp"
using namespace Kadath::Margherita ;


// the space type
typedef Kadath::Space_polar_adapted space_t;

// specialized quantities for a BNS system
struct ns_isotropic_vars_t : public Kadath::vars_base_t<ns_isotropic_vars_t> {};
// define the actual quantities and their order in the file!
template<> Kadath::var_vector Kadath::vars_base_t<ns_isotropic_vars_t>::vars = {
  {"lap_Aterm", SCALAR},
  {"nu", SCALAR},
  {"logh", SCALAR},
  {"lap_Bterm", SCALAR},
  {"lap_wterm", SCALAR},
  {"Omega", SCALAR},
};

class ns_isotropic_reader_t : public Kadath::python_reader_t<space_t, ns_isotropic_vars_t> {
  std::string config_filename;
  kadath_config_boost<BCO_ISO_NS_INFO> bconfig;

  public:
  ns_isotropic_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, ns_isotropic_vars_t>(filename),
                                             config_filename(filename.substr(0,filename.size()-3)+"info"),
                                             bconfig(config_filename) {
    // setup eos to before calling solver
    const double h_cut = bconfig.eos<double>(HCUT);
    const std::string eos_file = bconfig.eos<std::string>(EOSFILE);
    const std::string eos_type = bconfig.eos<std::string>(EOSTYPE);

    if(eos_type == "Cold_PWPoly") {
      using eos_t = Kadath::Margherita::Cold_PWPoly;

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut);
      this->compute_defs<eos_t>();
    } else if(eos_type == "Cold_Table") {
      using eos_t = Kadath::Margherita::Cold_Table;

      const int interp_pts = (bconfig.eos<int>(INTERP_PTS) == 0) ? \
                              2000 : bconfig.eos<int>(INTERP_PTS);

      EOS<eos_t,PRESSURE>::init(eos_file, h_cut, interp_pts);
      this->compute_defs<eos_t>();
    }
    else { 
      std::cerr << eos_type << " is not recognized.\n";
      std::_Exit(EXIT_FAILURE);
    }
    ns_Configurator_reader_t pybconfig(config_filename);
    config = pybconfig.config;
    // end eos setup and solver
   
  }

  template <typename eos_t>
  void compute_defs() {
    Kadath::Scalar const & lap_Aterm = extractField<Kadath::Scalar>("lap_Aterm");
    Kadath::Scalar const & nu = extractField<Kadath::Scalar>("nu");
    Kadath::Scalar const & lap_Bterm = extractField<Kadath::Scalar>("lap_Bterm");
    Kadath::Scalar const & lap_wterm = extractField<Kadath::Scalar>("lap_wterm");
    Kadath::Scalar const & Omega = extractField<Kadath::Scalar>("Omega");
    Kadath::Scalar const & logh = extractField<Kadath::Scalar>("logh");

  	int ndom = space.get_nbr_domains() ;

  	double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);

  	// Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;
    // define numerical constants
    syst.add_cst("4piG", bconfig(BCO_PARAMS::BCO_QPIG));

    // Fields - must be initialized before common setup
    syst.add_cst("H", logh);
    syst.add_cst("nu", nu);
    syst.add_cst("lapAterm", lap_Aterm);
    syst.add_cst("lapBterm", lap_Bterm);
    syst.add_cst("wrsint", lap_wterm);
    syst.add_cst("ome", Omega);

    syst.add_cst("omec", bconfig(OMEGA));

    syst.add_def("N = exp(nu)");
    syst.add_def("A = exp(lapAterm - nu)");
    syst.add_def("B = (divrsint(lapBterm) + 1) / N");
    syst.add_def("w = divrsint(wrsint)");
    syst.add_def("Fomega = B^2 * multrsint(multrsint(ome - w)) "
                        "/ (N^2 - multrsint(B * (ome - w))^2)");

    // define quantity to be integrated at infinity
    // two (in this case) equivalent definitions of ADM mass
    // as well as the Komar mass
    syst.add_def(ndom - 1, "intMadm = - (dr(A^2 + B^2) + divr(B^2 - A^2))  / 4 / 4piG ");
    syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");
    syst.add_def(ndom - 1, "intJ = -multrsint(multrsint(dr(w))) / 4 / 4piG");
    
    // enthalpy from the logarithmic enthalpy, the latter is the actual variable in this system
    syst.add_def("h = exp(H)");
  
    // define the EOS operators
    Param p;
    syst.add_ope("eps", &EOS<eos_t,EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t,PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t,DENSITY>::action, &p);
  
    // define rest-mass density, internal energy and pressure through the enthalpy
    syst.add_def("rho = rho(h)");
    syst.add_def("eps = eps(h)");
    syst.add_def("press = press(h)");

    // definition to rescale the equations
    // delta = p / rho
    syst.add_def("delta = h - eps - 1.");

    syst.add_def("U = multrsint(B / N * (ome - w))");
    syst.add_def("Usq = U*U");
    syst.add_def("Wsq = 1 / (1 - Usq)");
    syst.add_def("W = sqrt(Wsq)");

    // Avoid excision region (d=0,1) and compactified (d=ndom-1)
    // for(auto d = 2; d < ndom-1; ++d) {
    //   syst.add_def(d, "drP = dr(P)");
    //   syst.add_def(d, "ddrP = dr(drP)");
    // }

    for (int d = 0; d < ndom; d++) {
      
      switch (d) {
      // in the star the constraint equations are sourced by the matter
      case 0:
      case 1:

        // sources
        syst.add_def(d, "E = Wsq * press * h - press * delta");
        syst.add_def(d, "Srrtt = press * delta");
        syst.add_def(d, "pphi = B * (E + Srrtt) * U");
        syst.add_def(d, "Spp = delta * press * (1 + Usq) + E * Usq");
        syst.add_def(d, "S = 2 * Srrtt + Spp");
  
        // constraint equations
        syst.add_def(d, "eqnu  = delta * lap(nu) + delta * scal(grad(nu), grad(nu + log(B))) "
                              "- delta * multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w)) "
                              "- 4piG * A^2 * (E + S)");
        syst.add_def(d, "eqAterm = delta * lap2(lapAterm) + delta * scal(grad(nu), grad(nu))"
                        "- 3 * delta * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))"
                        "- 2 * 4piG * A^2 * Spp");
        syst.add_def(d, "eqBterm = delta * lap2(lapBterm) - 2 * 4piG * N * A^2 * multrsint(B) * (2 * Srrtt)");
        syst.add_def(d, "eqwrsint = delta * lap(wrsint) - delta * multrsint(scal(grad(w), grad(nu - 3 * log(B))))"
                            "+ 4 * 4piG * N * A^2 / B * pphi");
  
        // definition for the baryonic mass integral
        syst.add_def(d, "intMb = W * rho * A^2 * B * 4piG / 2");
        syst.add_def(d, "intH  = W * H * A^2 * B * 4piG / 2") ;
              break;
        // outside the matter is absent and the sources are zero
        default:

          syst.add_def(d, "eqnu  = lap(nu) + scal(grad(nu), grad(nu + log(B))) "
                          "- multrsint(multrsint(B^2)) / 2 / N^2 * scal(grad(w), grad(w))");
          syst.add_def(d, "eqAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))"
                    "- 3 * multrsint(multrsint(B^2)) / 4 / N^2 * scal(grad(w), grad(w))");
          syst.add_def(d, "eqBterm = lap2(lapBterm)");
          syst.add_def(d, "eqwrsint = lap(wrsint) - multrsint(scal(grad(w), grad(nu - 3 * log(B))))");
          break;
      }
    }
    auto add_surf_integ = [&](auto varstr, auto defstr, auto dom, auto bc) {
      vars[varstr]  = syst.get_space().get_domain(dom)->integ(
        syst.give_val_def(defstr)()(dom), bc
      );
    };
    add_surf_integ("Jadm", "intJ"   , ndom-1, OUTER_BC);
    add_surf_integ("Madm", "intMadm"   , ndom-1, OUTER_BC);
    add_surf_integ("Mk", "intMk"   , ndom-1, OUTER_BC);

    // Populate vars dictionary
    // FUKA_Syst_tools::syst_vars(vars, syst);
    // FUKA_Syst_tools::syst_vars_hydro(vars, syst);
    // FUKA_Syst_tools::export_radii(space, vars, 0, ndom-1, "NS_R");
    // FUKA_Syst_tools::dict_add_vector_cmp(
    //   syst, vars, "shift", Tensor(shift)
    // );

    // double Madm = boost::python::extract<double>(vars["Madm"]);
    FUKA_Syst_tools::syst_vars_NS_isotropic(vars, syst, 2);
    FUKA_Syst_tools::syst_add_resolution_list(space, vars);
    vars["nc"] = EOS<eos_t,DENSITY>::get(bconfig(BCO_PARAMS::HC));
    vars["hc"] = bconfig(BCO_PARAMS::HC);
  }
};

BOOST_PYTHON_MODULE(_ns_isotropic_reader)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    Kadath::constructPythonReader<ns_isotropic_reader_t>("ns_isotropic_reader");
}