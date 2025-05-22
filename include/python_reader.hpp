/*
    Copyright 2018 Ludwig Jens Papenfort

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PYTHON_READER_HH
#define PYTHON_READER_HH

#include <iostream>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <string>
#include <vector>
#include <utility>

#include "EOS/EOS.hh"
#include "utilities.hpp"

namespace Kadath {
/*
  These classes make it possible to create and use
  Kadath scalar and tensor fields directly from
  python, for analysis purposes.

  See the codes directory for some examples.
*/

enum quantity_type {CONSTANT, SCALAR, VECTOR};

typedef std::pair<std::string, quantity_type> var_t;
typedef std::vector<var_t> var_vector;

// base class for any collection of quantities defined on some space
template<typename T>
struct vars_base_t {
  static var_vector vars;
};

// in any case, get storage for them
template<typename T> var_vector vars_base_t<T>::vars;

// general reader class to be used within python
// fills the internal dictionaries using a quantities type
template<typename space_t, typename vars_t>
class python_reader_t {
  protected:
  // the file from which gets read
  FILE * file;

  // the space on which everything is defined
  space_t space;

  public:
  boost::python::dict vars;
  boost::python::dict config;

  // general construcor, filling the dicts
  python_reader_t(std::string const & filename)
    : file(fopen(filename.c_str(), "r")), space(file) {
    for(auto q : vars_t::vars) {
      switch(q.second) {
        case CONSTANT:
          // read in constant, i.e. single scalars
          double tmp;
          fread_be(&tmp, sizeof(double), 1, file);
          vars[q.first] = tmp;
          break;
        case SCALAR:
          // read in scalar field
          vars[q.first] = Scalar(space, file);
          break;
        case VECTOR:
          // read in vector field
          vars[q.first] = Vector(space, file);
          break;
      }
    }
    // close the file again, we don't need it anymore
    fclose(file);
  }

  template<typename T>
  T const & extractField(std::string const & fieldname) {
    return boost::python::extract<T const &>(vars[fieldname]);
  }

  template<typename T>
  T extractConstant(std::string const & cname) {
    return boost::python::extract<T>(vars[cname]);
  }

  // point-wise tensor field exporter interface for python
  boost::python::list getFieldValues(std::string const & fieldname, boost::python::list const & coord_list, int bound = -1) {
    // list of values to return
    boost::python::list values;

    // extract field as a general tensor object
    Tensor const & tensor = extractField<Tensor const &>(fieldname);

    // get number of dimensions and construct index
    int ndim = tensor.get_ndim();
    auto ind = Index(Dim_array(ndim));

    // loop through all given coords
    for(int i = 0; i < boost::python::len(coord_list); ++i) {
      // extract coords
      boost::python::list coords = boost::python::extract<boost::python::list>(coord_list[i]);

      // check if we got the expected number of coords
      if(boost::python::len(coords) != ndim) {
        std::cerr << "ERROR: Got wrong number of coordinates, expecting exactly " << ndim << std::endl;
        return values;
      }

      // construct Kadath point
      Point abs_coords(ndim);
      for(int j = 0; j < ndim; ++j) {
        abs_coords.set(j+1) = boost::python::extract<double>(coords[j]);
      }

      // add values to output list
      if(bound == -1) {
        if(tensor.get_valence() == 0)
          values.append(tensor(ind).val_point(abs_coords));
        else {
          boost::python::list comps;
          do {
            comps.append(tensor(ind).val_point(abs_coords));
          } while(ind.inc());
          values.append(comps);
        }
      }
      else {
        if(tensor.get_valence() == 0)
          values.append(tensor(ind).val_point_bound(abs_coords,bound));
        else {
          boost::python::list comps;
          do {
            comps.append(tensor(ind).val_point_bound(abs_coords,bound));
          } while(ind.inc());
          values.append(comps);
        }
      }
    }
    return values;
  }

};

template<class eos_t>
struct PygetEOSValues {
  // point-wise EOS exporter interface for python
  // boost::python::list getEOSValues(boost::python::list const & coord_list) {
  template<class reader_t>
  boost::python::list operator()(boost::python::list const & coord_list, reader_t * this_reader) {
    auto logh_values = this_reader->getFieldValues("logh", coord_list);

    // list of EOS values to return
    boost::python::list values;

    for(int i = 0; i < len(coord_list); ++i) {
      double h = std::exp(boost::python::extract<double>(logh_values[i]));

      boost::python::list eos_tuple;

      eos_tuple.append(EOS<eos_t, DENSITY>::get(h));
      eos_tuple.append(EOS<eos_t, EPSILON>::get(h));
      eos_tuple.append(EOS<eos_t, PRESSURE>::get(h));

      values.append(eos_tuple);
    }

    return values;
  }
};

template<typename space_t>
void initPythonBinding() {
  using namespace boost::python;

  // this is needed, so that python lets us read in the fields
  class_<Kadath::var_t>("KadathVarType");
  // this is needed, so that python understands what the field types are
  class_<Kadath::Tensor>("KadathTensorType", init<Kadath::Tensor const &, bool>());
  class_<Kadath::Scalar, bases<Kadath::Tensor>>("KadathScalarType", init<space_t const &, FILE*>());
  class_<Kadath::Vector, bases<Kadath::Tensor>>("KadathVectorType", init<space_t const &, FILE*>());
}

// dummy constructor function, defining readers through boost python
template<typename reader_t>
void constructPythonReader(std::string reader_name) {
  using namespace boost::python;

  auto reader = class_<reader_t>(reader_name.c_str(), init<std::string>());
  reader.def("getFieldValues", &reader_t::getFieldValues);
  reader.def("getEOSValues", &reader_t::getEOSValues);
  reader.def_readonly("vars", &reader_t::vars);
  reader.def_readonly("config", &reader_t::config);
}


} // namespace
#endif
