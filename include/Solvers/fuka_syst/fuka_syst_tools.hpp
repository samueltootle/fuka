#pragma once
#include <string>
#include "exporter_utilities.hpp"

/**
 * \addtogroup Syst_tools
 * \ingroup FUKA
 * Centralized tools to populate and manipulate System_of_eqs
 * @{*/

namespace Kadath {
namespace FUKA_Syst_tools {

/**
 * @brief Generate a scalar field of domain outer radii for a range
 * of continuous domain indicies \in [dom_min, dom_max]
 *
 * @tparam space_t Numerical space type
 * @tparam dict_t Python dictionary type
 * @param space numerical space
 * @param vars python dictionary to add to
 * @param dom_min Minimum domain index
 * @param dom_max Maximum domain index
 * @param forestr String to append at the start
 */
template <class space_t, class dict_t>
void export_radii(space_t& space,
                  dict_t& vars,
                  int const dom_min,
                  int const dom_max,
                  std::string forestr) {
    int const ndom = space.get_nbr_domains();
    Scalar space_radius(space);
    space_radius.annule_hard();
    for (int d = 0; d < ndom; ++d) {
        space_radius.set_domain(d) = space.get_domain(d)->get_radius();
    }
    space_radius.std_base();

    unsigned int cnt = 1;
    for (int i = dom_min; i < dom_max; ++i) {
        Index pos(
            space.get_domain(i)->get_radius().get_conf().get_dimensions());
        pos.set(0) = space.get_domain(i)->get_nbr_points()(0) - 1;
        vars[forestr + std::to_string(cnt)] = space_radius(i)(pos);
        cnt++;
    }
}

/**
 * @brief Add rank 2 tensor data to python dictionary componentwise
 * along with extracting the trace from System of equations
 *
 * @tparam dict_t Python dictionary type
 * @param syst System of equations
 * @param vars python dictionary to add fields to
 * @param var System of equations string definition of tensor
 * @param field Tensor field to extract scalar fields from
 */
template <class dict_t>
inline void dict_add_tensor_cmp(System_of_eqs& syst,
                                dict_t& vars,
                                std::string var,
                                Tensor field) {
    int c = 0;
    vars[var.c_str()] = syst.give_val_def(("Trace" + var).c_str());
    for (std::string coord : {"_11", "_12", "_13", "_22", "_23", "_33"}) {
        auto tidx = export_utils::R2TensorSymmetricIndices[c];
        Array<int> ind(field.indices(tidx));
        vars[(var + coord).c_str()] = field(ind);
        c++;
    }
};

/**
 * @brief Add rank 1 tensor data to python dictionary componentwise
 *
 * @tparam dict_t Python dictionary type
 * @param syst System of equations
 * @param vars python dictionary to add fields to
 * @param var System of equations string definition of tensor
 * @param field Tensor field to extract scalar fields from
 */
template <class dict_t>
inline void dict_add_vector_cmp(System_of_eqs& syst,
                                dict_t& vars,
                                std::string var,
                                Tensor field) {
    int c = 0;
    for (std::string coord : {"_1", "_2", "_3"}) {
        Array<int> ind(field.indices(c));
        vars[(var + coord).c_str()] = field(ind);
        c++;
    }
};

/**
 * @brief Helper function to get a vector of continuous domain indicies
 * up to, but not including dom_max
 *
 * @param dom_min minimum domain index
 * @param dom_max maximum domain index
 * @return std::vector<int> vector of domains
 */
inline std::vector<int> vector_of_domains(int const dom_min,
                                          int const dom_max) {
#include <numeric>
    assert(dom_min <= dom_max);
    // dom_max is not added to the list
    std::vector<int> doms(dom_max - dom_min);
    std::iota(doms.begin(), doms.end(), dom_min);
    return doms;
}

/** @}*/
}  // namespace FUKA_Syst_tools
}  // namespace Kadath