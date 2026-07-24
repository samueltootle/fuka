/*
 * Copyright 2022
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
#include <string>
#include "Configurator/config_binary.hpp"
#include "Solvers/bbh_xcts/bbh_exporter.hpp"
#include "bco_utilities.hpp"
#include "coord_fields.hpp"
#include "include/fuka_py.hpp"
#include "kadath_bin_bh.hpp"
#include "python_reader.hpp"
typedef Kadath::Space_bin_bh space_t;

struct bbh_vars_t : public Kadath::vars_base_t<bbh_vars_t> {};

template <>
Kadath::var_vector Kadath::vars_base_t<bbh_vars_t>::vars = {
    {"conf", SCALAR},
    {"lapse", SCALAR},
    {"shift", VECTOR},
};

class bbh_reader_t : public Kadath::python_reader_t<space_t, bbh_vars_t> {
    std::string config_filename;
    kadath_config_boost<BIN_INFO> bconfig;
    using exporter_t = Kadath::FUKA_Solvers::CFMS_BBH_Exporter;
    exporter_t exporter;

   public:
    bbh_reader_t(std::string const filename)
        : Kadath::python_reader_t<space_t, bbh_vars_t>(filename),
          config_filename(filename.substr(0, filename.size() - 3) + "info"),
          bconfig(config_filename),
          exporter(config_filename) {
        this->compute_defs();
        bin_Configurator_reader_t pybconfig(config_filename);
        config = pybconfig.config;
    }

    void compute_defs() {
        using namespace Kadath;

        Kadath::Scalar const& conf = extractField<Kadath::Scalar>("conf");
        Kadath::Scalar const& lapse = extractField<Kadath::Scalar>("lapse");
        Kadath::Vector const& shift = extractField<Kadath::Vector>("shift");

        Base_tensor basis(shift.get_basis());
        Metric_flat fmet(space, basis);
        int ndom = space.get_nbr_domains();

        double xc1 = bco_utils::get_center(space, space.BH1);
        double xc2 = bco_utils::get_center(space, space.BH2);
        double xo = bco_utils::get_center(space, ndom - 1);

        std::vector<int> BH1_Domains{
            FUKA_Syst_tools::vector_of_domains(space.BH1 + 2, space.BH2)};
        std::vector<int> BH2_Domains{
            FUKA_Syst_tools::vector_of_domains(space.BH2 + 2, space.OUTER)};
        std::vector<int> vac_Domains{
            FUKA_Syst_tools::vector_of_domains(0, ndom)};

        // setup coordinate vector fields for System_of_eqs
        CoordFields<Space_bin_bh> cf_generator(space);
        vec_ary_t coord_vectors = default_binary_vector_ary(space);
        update_fields(cf_generator, coord_vectors, {}, xo, xc1, xc2);

        Vector CART(space, CON, basis);
        CART = cf_generator.cart();
        // end coordinate field setup

        // Setup system of equations and definitions
        System_of_eqs syst(space, 0, ndom - 1);
        fmet.set_system(syst, "f");

        // some constants
        syst.add_cst("4piG", bconfig(QPIG));

        // Fields
        syst.add_cst("P", conf);
        syst.add_cst("N", lapse);
        syst.add_cst("bet", shift);

        // Binary Quantities
        syst.add_cst("ome", bconfig(GOMEGA));
        syst.add_cst("xaxis", bconfig(COM));
        syst.add_cst("yaxis", bconfig(COMY));

        // Component quantities
        syst.add_cst("omesm", bconfig(OMEGA, BCO1));
        syst.add_cst("omesp", bconfig(OMEGA, BCO2));

        for (auto d = 0; d < ndom; ++d) {
            if ((d != space.BH1 && d != space.BH1 + 1) &&
                (d != space.BH2 && d != space.BH2 + 1) && d < space.OUTER) {
                syst.add_def(d, "drP = dr(P)");
                syst.add_def(d, "ddrP = dr(drP)");
            }
        }
        // Initialize syst with constants, fields, etc...
        FUKA_Syst_tools::syst_init_binary(syst, coord_vectors, bconfig);
        FUKA_Syst_tools::syst_init_inspiral(syst, bconfig, CART);

        // Quasi-local spin definitions
        FUKA_Syst_tools::syst_init_quasi_local_defs(syst,
                                                    BH1_Domains,
                                                    "mm",
                                                    "sm");
        FUKA_Syst_tools::syst_init_quasi_local_defs(syst,
                                                    BH2_Domains,
                                                    "mp",
                                                    "sp");
        // Add vacuum constraint definitions
        FUKA_Syst_tools::syst_init_eqdefs_vac(syst, vac_Domains);
        FUKA_Syst_tools::syst_init_contraction_defs_vac(syst);

        // Populate vars dictionary
        FUKA_Syst_tools::syst_vars(vars, syst);
        FUKA_Syst_tools::syst_vars_BH(vars, syst, space.BH1 + 2, "BH1_");
        FUKA_Syst_tools::syst_vars_BH(vars, syst, space.BH2 + 2, "BH2_");
        FUKA_Syst_tools::export_radii(space,
                                      vars,
                                      space.BH1 + 1,
                                      space.BH2,
                                      "BH1_R");
        FUKA_Syst_tools::export_radii(space,
                                      vars,
                                      space.BH2 + 1,
                                      space.OUTER,
                                      "BH2_R");

        // BCO coordinate centers
        vars["xc1"] = xc1;
        vars["xc2"] = xc2;
        vars["xo"] = xo;

        FUKA_Syst_tools::dict_add_vector_cmp(syst,
                                             vars,
                                             "shift",
                                             Tensor(shift));

        //create old radius scalar field
        Scalar space_radius(space);
        space_radius.annule_hard();
        for (int d = 0; d < ndom; ++d) {
            space_radius.set_domain(d) = space.get_domain(d)->get_radius();
        }
        space_radius.std_base();

        auto find_rmin = [&](int d) {
            Index pos(space.get_domain(d)->get_nbr_points());
            pos.set_start();

            double rmin = space_radius(d)(pos);
            do {
                pos.set(0) = space.get_domain(d)->get_nbr_points()(0) - 1;
                double r = space_radius(d)(pos);

                if (r < rmin)
                    rmin = r;
            } while (pos.inc());
            return rmin;
        };

        // create a list of domain radii
        boost::python::list d_outer_r;
        for (int d = 0; d < ndom; ++d)
            d_outer_r.append(find_rmin(d));

        vars["dom_r"] = d_outer_r;

        // this was used to do a colored plot of the domains
        // like kids coloring books - color the numbered areas.
        Scalar dom_colors(conf);
        dom_colors.annule_hard();
        int c = 1;
        for (int d = 0; d < ndom; ++d) {
            //set BH interiors to the same color
            if (d == space.BH1 || d == space.BH1 + 1 || d == space.BH2 ||
                d == space.BH2 + 1) {
                dom_colors.set_domain(d) = 0;
            }
            //set inner_adapted domains to the same color
            else if (d == space.BH1 + 2 || d == space.BH2 + 2) {
                dom_colors.set_domain(d) = 1;
            }
            //set chi_first domains to the same color
            else if (d == space.OUTER || d == space.OUTER + 4) {
                dom_colors.set_domain(d) = 2;
            }
            //rect domains
            else if (d == space.OUTER + 1 || d == space.OUTER + 3) {
                dom_colors.set_domain(d) = 3;
            }
            //eta + shells + compactified domains
            else {
                dom_colors.set_domain(d) = 3 + c;
                c++;
            }
        }
        dom_colors.std_base();
        vars["dom_color_chart"] = dom_colors;
        vars["rfield"] = space_radius;
    }

    boost::python::list getExporterFieldValues__cartesian(
        std::string const& fieldname,
        boost::python::list const& coord_list) {
        // list of values to return
        boost::python::list values;

        // From CPPReference
        auto str_tolower = [](std::string s) {
            std::transform(s.begin(),
                           s.end(),
                           s.begin(),
                           [](unsigned char c) {
                               return std::tolower(c);
                           }  // correct
            );
            return s;
        };
        auto key{str_tolower(fieldname)};
        size_t idx;
        if (auto search = exporter.output_var_map.find(key);
            search == exporter.output_var_map.end()) {
            std::string msg{"Invalid output fieldname pass: " + fieldname};
            throw std::invalid_argument(msg.c_str());
        } else {
            idx = exporter.output_var_map[key];
            cout << key << ": " << idx << endl;
        }

        // loop through all given coords
        for (int i = 0; i < boost::python::len(coord_list); ++i) {
            // extract coords
            boost::python::list coords =
                boost::python::extract<boost::python::list>(coord_list[i]);
            auto output_vars = exporter.export_pointwise(
                boost::python::extract<double>(coords[0]),
                boost::python::extract<double>(coords[1]),
                boost::python::extract<double>(coords[2]));
            values.append(output_vars[idx]);
        }
        return values;
    }

    boost::python::list getExporterKeys() {
        boost::python::list values;
        for (auto t : exporter.output_var_map) {
            values.append(t.first);
        }
        return values;
    }

    boost::python::dict getallExporterFieldValues__cartesian_pointwise(
        boost::python::list const& coord) {
        // list of values to return
        boost::python::dict values;

        if (boost::python::len(coord) > 3) {
            std::string msg{
                "getallExporterFieldValues_pointwise accepts a single "
                "coordinate "
                "only!"};
            throw std::invalid_argument(msg.c_str());
        }

        auto output_vars =
            exporter.export_pointwise(boost::python::extract<double>(coord[0]),
                                      boost::python::extract<double>(coord[1]),
                                      boost::python::extract<double>(coord[2]));

        // loop through all given coords
        for (auto& kvp : exporter.output_var_map) {
            auto k = kvp.first;
            auto idx = kvp.second;

            // Create dictionary
            values[k] = output_vars[idx];
        }
        return values;
    }
};

// dummy constructor function, defining readers through boost python
template <typename reader_t>
void constructPythonReader_here(std::string reader_name) {
    using namespace boost::python;

    auto reader = class_<reader_t>(reader_name.c_str(), init<std::string>());
    reader.def("getFieldValues", &reader_t::getFieldValues);
    reader.def("getExporterFieldValues__cartesian",
               &reader_t::getExporterFieldValues__cartesian);
    reader.def("getExporterKeys", &reader_t::getExporterKeys);
    reader.def("getallExporterFieldValues__cartesian_pointwise",
               &reader_t::getallExporterFieldValues__cartesian_pointwise);
    reader.def_readonly("vars", &reader_t::vars);
    reader.def_readonly("config", &reader_t::config);
}

BOOST_PYTHON_MODULE(_bbh_reader) {
    // initialize python types
    Kadath::initPythonBinding<space_t>();
    constructPythonReader_here<bbh_reader_t>("bbh_reader");
}
