#pragma once
#include "Solvers/fuka_syst/fuka_syst_tools.hpp"

namespace Kadath {
namespace FUKA_Syst_tools {

template <class dict_t>
void syst_vars(dict_t& vars, System_of_eqs& syst) {
#ifdef DEBUG
    std::cout << "Loading standard vars into dictionary\n";
#endif
    auto const ndom = syst.get_space().get_nbr_domains();
    // add vars to boost dictionary
    vars["drPsi"] = syst.give_val_def("drP");
    vars["ddrPsi"] = syst.give_val_def("ddrP");
    vars["cPsi"] = syst.give_val_def("eqP");
    vars["cLapsePsi"] = syst.give_val_def("eqNP");

    Kadath::Vector beta_res(syst.give_val_def("eqbet"));

    vars["cShift_1"] = beta_res(1);
    vars["cShift_2"] = beta_res(2);
    vars["cShift_3"] = beta_res(3);

    auto add_surf_integ = [&](auto varstr, auto defstr, auto dom, auto bc) {
        vars[varstr] =
            syst.get_space().get_domain(dom)->integ(syst.give_val_def(defstr)()(
                                                        dom),
                                                    bc);
    };
    add_surf_integ("Jadm", "intJ", ndom - 1, OUTER_BC);
    add_surf_integ("Madm", "intMadm", ndom - 1, OUTER_BC);
    add_surf_integ("Mk", "intMk", ndom - 1, OUTER_BC);
    add_surf_integ("Px", "intPx", ndom - 1, OUTER_BC);
    add_surf_integ("Py", "intPy", ndom - 1, OUTER_BC);
    add_surf_integ("Pz", "intPz", ndom - 1, OUTER_BC);

    FUKA_Syst_tools::dict_add_tensor_cmp(syst,
                                         vars,
                                         "A",
                                         syst.give_val_def("A"));
    FUKA_Syst_tools::dict_add_vector_cmp(syst,
                                         vars,
                                         "B",
                                         syst.give_val_def("B"));

    // Conformal Aij contractions
    vars["Axx"] = syst.give_val_def("Axx");
    vars["Ayy"] = syst.give_val_def("Ayy");
    vars["Azz"] = syst.give_val_def("Azz");
    vars["Axy"] = syst.give_val_def("Axy");
    vars["Ayz"] = syst.give_val_def("Ayz");
    vars["Axz"] = syst.give_val_def("Axz");
    vars["A"] = syst.give_val_def("TraceA");

    // inertial shift contractions
    vars["betx"] = syst.give_val_def("betx");
    vars["bety"] = syst.give_val_def("bety");
    vars["betz"] = syst.give_val_def("betz");

    // total shift contractions
    vars["Bx"] = syst.give_val_def("Bx");
    vars["By"] = syst.give_val_def("By");
    vars["Bz"] = syst.give_val_def("Bz");
}

template <class dict_t>
void syst_vars_BH(dict_t& vars,
                  System_of_eqs& syst,
                  int const adapt_d,
                  std::string iden = {}) {
#ifdef DEBUG
    std::cout << "Loading BH vars into dictionary.\n";
#endif
    auto& space = syst.get_space();

    // Local angular momentum on the AH
    double S =
        space.get_domain(adapt_d)->integ(syst.give_val_def("intS")()(adapt_d),
                                         INNER_BC);
    vars[std::string{iden + "S"}.c_str()] = S;

    // Irreducible mass
    double Mirrsq =
        space.get_domain(adapt_d)->integ(syst.give_val_def("intMirrsq")()(
                                             adapt_d),
                                         INNER_BC);
    double Mirr = sqrt(Mirrsq);
    vars[std::string{iden + "Mirr"}.c_str()] = Mirr;

    // Area radius
    vars[std::string{iden + "ArealR"}.c_str()] = 2. * Mirr;

    double Mch = std::sqrt(Mirrsq + S * S / 4. / Mirrsq);
    vars[std::string{iden + "Mch"}.c_str()] = Mch;

    double Chi = S / Mch / Mch;
    vars[std::string{iden + "Chi"}.c_str()] = Chi;
}

template <class space_t, class dict_t>
void syst_add_resolution_list(space_t& space, dict_t& vars) {
#include <boost/python.hpp>
    boost::python::list all_res;
    auto ndom = space.get_nbr_domains();
    for (auto d = 0; d < ndom; ++d) {
        boost::python::list dom_res;
        auto npts(space.get_domain(d)->get_nbr_points());
        auto ndim = npts.get_ndim();
        for (int i = 0; i < ndim; ++i) {
            dom_res.append(npts(i));
        }
        all_res.append(dom_res);
    }
    vars["domain_resolutions"] = all_res;
}
}  // namespace FUKA_Syst_tools
}  // namespace Kadath

#include "fuka_syst_vars_hydro.hpp"
