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
#include "Solvers/ns_isotropic/ns_isotropic_exporter.hpp"
using namespace Kadath::Margherita ;


// the space type
typedef Kadath::Space_polar_adapted space_t;

// specialized quantities for a non-rotating isotropic NS
struct ns_isotropic_vars_t : public Kadath::vars_base_t<ns_isotropic_vars_t> {};
// define the actual quantities and their order in the file!
template<> Kadath::var_vector Kadath::vars_base_t<ns_isotropic_vars_t>::vars = {
  {"lap_Aterm", SCALAR},
  {"nu", SCALAR},
  {"logh", SCALAR},
  {"lap_Bterm", SCALAR},
};

class ns_isotropic_reader_t : public Kadath::python_reader_t<space_t, ns_isotropic_vars_t> {
  using exporter_t = Kadath::FUKA_Solvers::CFMS_NS_ISO_Exporter;
  std::string config_filename;
  kadath_config_boost<BCO_ISO_NS_INFO> bconfig;
  exporter_t exporter;

  public:
  ns_isotropic_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, ns_isotropic_vars_t>(filename),
                                             config_filename(filename.substr(0,filename.size()-3)+"info"),
                                             bconfig(config_filename), exporter(config_filename) {
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
      std::string msg{"Invalid eos type: "+eos_type};
      throw std::invalid_argument(msg.c_str());
    }
    ns_Configurator_reader_t pybconfig(config_filename);
    config = pybconfig.config;
    // end eos setup and solver
   
  }

  template <typename eos_t>
  void compute_defs() {
    Kadath::Scalar const & lap_Aterm = extractField<Kadath::Scalar>("lap_Aterm");
    Kadath::Scalar const & nu = extractField<Kadath::Scalar>("nu");
    Kadath::Scalar const & logh = extractField<Kadath::Scalar>("logh");
    Kadath::Scalar const & lap_Bterm = extractField<Kadath::Scalar>("lap_Bterm");

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

    syst.add_def("N = exp(nu)");
    syst.add_def("A = exp(lapAterm - nu)");

    // define quantity to be integrated at infinity
    // two (in this case) equivalent definitions of ADM mass
    // as well as the Komar mass
    syst.add_def(ndom - 1, "intMadm = -dr(A)  / 4piG ");
    syst.add_def(ndom - 1, "intMk = dr(N)  / 4piG");
    
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
        syst.add_def(d, "E = press * (1 + eps)");
        syst.add_def(d, "S = delta * 3 * press");
        syst.add_def(d, "Spp = press * delta");
  
        // constraint equations
        syst.add_def(d, "eqnu = delta * ( lap(nu) + dr(nu) * dr(lapAterm) ) - 4piG * A^2 * (E + S)") ;
        syst.add_def(d, "eqAterm = delta * ( lap2(lapAterm) + dr(nu) * dr(nu) ) - 2 * 4piG * A^2 * Spp") ;
  
        // definition for the baryonic mass integral
        syst.add_def(d, "intMb = rho * A^3 * 4piG / 2");
        syst.add_def(d, "intH  = H * A^3 * 4piG / 2") ;
              break;
        // outside the matter is absent and the sources are zero
        default:

          syst.add_def(d, "eqnu = lap(nu) + scal(grad(nu), grad(lapAterm))") ;
          syst.add_def(d, "eqAterm = lap2(lapAterm) + scal(grad(nu), grad(nu))") ;
          break;
      }
    }
    auto add_surf_integ = [&](auto varstr, auto defstr, auto dom, auto bc) {
      vars[varstr]  = syst.get_space().get_domain(dom)->integ(
        syst.give_val_def(defstr)()(dom), bc
      );
    };
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
  
    auto add_from_def = [&](std::string in_str, std::string out_str="") {
      if(out_str == "") out_str = in_str;
      Scalar tmp(syst.give_val_def(in_str.c_str()));
      tmp.coef_i();
      tmp.std_base();
      vars[out_str.c_str()] = tmp;
    };
    add_from_def("A");
    add_from_def("N");
    add_from_def("eqnu", "cLapse");
    add_from_def("eqAterm", "cA");

    auto npts = space.get_domain(1)->get_nbr_points();
    Index pos_eq (npts);
    pos_eq.set(0) = npts(0) - 1; /// Set to outer radius
    pos_eq.set(1) = npts(1) - 1; /// Set theta to be on the xy plane.
    auto B(syst.give_val_def("A")()(1));
    auto r(space.get_domain(1)->get_radius());
    double CR = B(pos_eq) * r(pos_eq);
    vars["CR"] = CR;
  }

  boost::python::list getExporterFieldValues__cartesian(std::string const & fieldname, boost::python::list const & coord_list) {
    // list of values to return
    boost::python::list values;

    
    // From CPPReference
    auto str_tolower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), 
                      [](unsigned char c){ return std::tolower(c); } // correct
                      );
        return s;
    };
    auto key{str_tolower(fieldname)};
    size_t idx;
    if (auto search = exporter.output_var_map.find(key); search == exporter.output_var_map.end()) {
      std::string msg{"Invalid output fieldname pass: "+fieldname};
      throw std::invalid_argument(msg.c_str());
    } else {
      idx = exporter.output_var_map[key];
      cout << key << ": " << idx << endl;
    }
    
    // loop through all given coords
    for(int i = 0; i < boost::python::len(coord_list); ++i) {
      // extract coords
      boost::python::list coords = boost::python::extract<boost::python::list>(coord_list[i]);
      auto output_vars = exporter.export_pointwise(
        boost::python::extract<double>(coords[0]), 
        boost::python::extract<double>(coords[1]), 
        boost::python::extract<double>(coords[2])
      );
      values.append(output_vars[idx]);
    }
    return values;
  }
    boost::python::list getExporterFieldValues__spherical(std::string const & fieldname, boost::python::list const & coord_list) {
    // list of values to return
    boost::python::list values;

    
    // From CPPReference
    auto str_tolower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), 
                      [](unsigned char c){ return std::tolower(c); } // correct
                      );
        return s;
    };
    auto key{str_tolower(fieldname)};
    size_t idx;
    if (auto search = exporter.output_var_map.find(key); search == exporter.output_var_map.end()) {
      std::string msg{"Invalid output fieldname pass: "+fieldname};
      throw std::invalid_argument(msg.c_str());
    } else {
      idx = exporter.output_var_map[key];
      cout << key << ": " << idx << endl;
    }
    
    // loop through all given coords
    for(int i = 0; i < boost::python::len(coord_list); ++i) {
      // extract coords
      boost::python::list coords = boost::python::extract<boost::python::list>(coord_list[i]);
      auto output_vars = exporter.export_pointwise(
        boost::python::extract<double>(coords[0]), 
        boost::python::extract<double>(coords[1]), 
        boost::python::extract<double>(coords[2])
      );
      values.append(output_vars[idx]);
    }
    return values;
  }
  boost::python::list getExporterKeys() {
    boost::python::list values;
    for(auto t : exporter.output_var_map) {
      values.append(t.first);
    }
    return values;
  }
  boost::python::dict getallExporterFieldValues__spherical_pointwise(boost::python::list const & coord) {
    // list of values to return
    boost::python::dict values;

    if (boost::python::len(coord) > 3) {
      std::string msg{"getallExporterFieldValues_pointwise accepts a single coordinate only!"};
      throw std::invalid_argument(msg.c_str());
    }

    auto output_vars = exporter.export_pointwise__spherical(
      boost::python::extract<double>(coord[0]), 
      boost::python::extract<double>(coord[1]), 
      boost::python::extract<double>(coord[2])
    );
    
    // loop through all given coords
    for(auto& kvp : exporter.output_var_map) {
      auto k = kvp.first;
      auto idx = kvp.second;

      // Create dictionary
      values[k] = output_vars[idx];
    }
    return values;
  }
  boost::python::dict getallExporterFieldValues__cartesian_pointwise(boost::python::list const & coord) {
    // list of values to return
    boost::python::dict values;

    if (boost::python::len(coord) > 3) {
      std::string msg{"getallExporterFieldValues_pointwise accepts a single coordinate only!"};
      throw std::invalid_argument(msg.c_str());
    }

    auto output_vars = exporter.export_pointwise(
      boost::python::extract<double>(coord[0]), 
      boost::python::extract<double>(coord[1]), 
      boost::python::extract<double>(coord[2])
    );
    
    // loop through all given coords
    for(auto& kvp : exporter.output_var_map) {
      auto k = kvp.first;
      auto idx = kvp.second;

      // Create dictionary
      values[k] = output_vars[idx];
    }
    return values;
  }
};

// dummy constructor function, defining readers through boost python
template<typename reader_t>
void constructPythonReader_here(std::string reader_name) {
  using namespace boost::python;

  auto reader = class_<reader_t>(reader_name.c_str(), init<std::string>());
  reader.def("getFieldValues", &reader_t::getFieldValues);
  reader.def("getEOSValues", &reader_t::getEOSValues);
  reader.def("getExporterFieldValues__cartesian", &reader_t::getExporterFieldValues__cartesian);
  reader.def("getExporterFieldValues__spherical", &reader_t::getExporterFieldValues__spherical);
  reader.def("getExporterKeys", &reader_t::getExporterKeys);
  reader.def("getallExporterFieldValues__cartesian_pointwise", &reader_t::getallExporterFieldValues__cartesian_pointwise);
  reader.def("getallExporterFieldValues__spherical_pointwise", &reader_t::getallExporterFieldValues__spherical_pointwise);
}

BOOST_PYTHON_MODULE(_ns_isotropic_norot_reader)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    constructPythonReader_here<ns_isotropic_reader_t>("ns_isotropic_norot_reader");
    
}