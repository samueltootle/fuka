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
#include "kadath_adapted.hpp"
#include "Configurator/config_bco.hpp"
#include "bco_utilities.hpp"
#include "EOS/EOS.hh"
#include "coord_fields.hpp"
#include "python_reader.hpp"
#include "include/fuka_py.hpp"
#include "Solvers/ns_3d_xcts/ns_exporter.hpp"
using namespace Kadath::Margherita ;


// the space type
typedef Kadath::Space_spheric_adapted space_t;

// specialized quantities for a NS system
struct ns_vars_t : public Kadath::vars_base_t<ns_vars_t> {};
// define the actual quantities and their order in the file!
template<> Kadath::var_vector Kadath::vars_base_t<ns_vars_t>::vars = {
  {"conf", SCALAR},
  {"lapse", SCALAR},
  {"shift", VECTOR},
  {"logh", SCALAR},
};

class ns_reader_t : public Kadath::python_reader_t<space_t, ns_vars_t> {
  std::string config_filename;
  kadath_config_boost<BCO_NS_INFO> bconfig;
  using exporter_t = Kadath::FUKA_Solvers::CFMS_NS_Exporter;
  exporter_t exporter;

  public:
  ns_reader_t(std::string const filename) : Kadath::python_reader_t<space_t, ns_vars_t>(filename),
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
      std::cerr << eos_type << " is not recognized.\n";
      std::_Exit(EXIT_FAILURE);
    }
    ns_Configurator_reader_t pybconfig(config_filename);
    config = pybconfig.config;
    // end eos setup and solver
   
  }

  template <typename eos_t>
  void compute_defs() {
    Kadath::Scalar const & conf = extractField<Kadath::Scalar>("conf");
    Kadath::Scalar const & lapse = extractField<Kadath::Scalar>("lapse");
    Kadath::Vector const & shift = extractField<Kadath::Vector>("shift");
    Kadath::Scalar const & logh = extractField<Kadath::Scalar>("logh");

    Base_tensor basis(shift.get_basis());
  	Metric_flat fmet (space, basis) ;
  	int ndom = space.get_nbr_domains() ;

    double xo  = bco_utils::get_center(space, ndom-1);
    std::vector<int> matter_Domains{0,1};
    std::vector<int> vac_Domains{
      FUKA_Syst_tools::vector_of_domains(2, ndom)
    };

  	// setup coordinate vector fields for System_of_eqs
    CoordFields<Space_spheric_adapted> cfields(space);
    vec_ary_t coord_vectors {default_co_vector_ary(space)};
    update_fields_co(cfields, coord_vectors, {}, xo);
    // end coordinate field setup

  	double loghc = bco_utils::get_boundary_val(0, logh, INNER_BC);

  	// Setup system of equations and definitions
    System_of_eqs syst (space, 0, ndom-1) ;
    fmet.set_system(syst, "f") ;

    // Fields - must be initialized before common setup
    syst.add_cst("P"  , conf) ;
    syst.add_cst("N"  , lapse) ;
    syst.add_cst("bet", shift) ;
    syst.add_cst("H"  , logh) ;

    syst.add_cst("ome", bconfig(OMEGA));
   
    // Avoid excision region (d=0,1) and compactified (d=ndom-1)
    for(auto d = 2; d < ndom-1; ++d) {
      syst.add_def(d, "drP = dr(P)");
      syst.add_def(d, "ddrP = dr(drP)");
    }

    // Add interface to the EOS framework as a Syst operator
    FUKA_Syst_tools::syst_init_defs_hydro<eos_t>(syst);

    // Initialize constants and constant fields
    FUKA_Syst_tools::syst_init_co(syst, coord_vectors, bconfig);
    FUKA_Syst_tools::syst_init_eqdefs_hydro_corot(syst, matter_Domains);
    FUKA_Syst_tools::syst_init_quasi_local_defs_hydro(syst, matter_Domains);
    // Define equations in relevant domains
    FUKA_Syst_tools::syst_init_quasi_local_defs(
      syst, vac_Domains,"mg","sm"
    );
    FUKA_Syst_tools::syst_init_eqdefs_vac(syst, vac_Domains);
    FUKA_Syst_tools::syst_init_contraction_defs_vac(syst);
    FUKA_Syst_tools::syst_init_contraction_defs_hydro(syst);

    // Populate vars dictionary
    FUKA_Syst_tools::syst_vars(vars, syst);
    FUKA_Syst_tools::syst_vars_hydro(vars, syst);
    FUKA_Syst_tools::export_radii(space, vars, 0, ndom-1, "NS_R");
    FUKA_Syst_tools::dict_add_vector_cmp(
      syst, vars, "shift", Tensor(shift)
    );

    double Madm = boost::python::extract<double>(vars["Madm"]);
    FUKA_Syst_tools::syst_vars_NS(vars, syst, 2, Madm, matter_Domains);
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

  boost::python::list getExporterKeys() {
    boost::python::list values;
    for(auto t : exporter.output_var_map) {
      values.append(t.first);
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
  reader.def("getExporterKeys", &reader_t::getExporterKeys);
  reader.def("getallExporterFieldValues__cartesian_pointwise", &reader_t::getallExporterFieldValues__cartesian_pointwise);
  reader.def_readonly("vars", &reader_t::vars);
  reader.def_readonly("config", &reader_t::config);  
}

BOOST_PYTHON_MODULE(_ns_reader)
{
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    constructPythonReader_here<ns_reader_t>("ns_reader");
}