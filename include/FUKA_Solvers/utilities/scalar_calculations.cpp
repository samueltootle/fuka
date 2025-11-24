#pragma once
namespace Kadath::FUKA_Solvers {
    template <class space_t, class metric_t>
    Scalar compute_drPsi(space_t& space,
                        Scalar& conf,
                        metric_t metric,
                        std::vector<int> excluded_doms,
                        int bound_dom = -1) {
    std::string const def_drP{"drP = dr(P)"};
    auto ndom{space.get_nbr_domains()};
    System_of_eqs syst(space);
    metric.set_system(syst, "f");
    syst.add_cst("P", conf);

    if (excluded_doms.size() == 0) {
        syst.add_def(def_drP.c_str());
    } else {
        auto last_dom = (bound_dom == -1) ? ndom : bound_dom;
        for (auto dom = 0; dom < last_dom; ++dom) {
        auto res = std::find(excluded_doms.begin(), excluded_doms.end(), dom);
        if (res == std::end(excluded_doms)) {
            syst.add_def(dom, def_drP.c_str());
        }
        }
    }
    Scalar field(syst.give_val_def("drP"));
    field.std_base();
    return field;
    }



    template <class space_t, class metric_t>
    Scalar compute_ddrPsi(space_t& space,
                        Scalar& conf,
                        metric_t metric,
                        std::vector<int> excluded_doms,
                        int bound_dom = -1) {
        std::string const def_drP{"drP = dr(P)"};
        std::string const def_drdrP{"ddrP = dr(drP)"};
        auto ndom{space.get_nbr_domains()};
        System_of_eqs syst(space);
        metric.set_system(syst, "f");
        syst.add_cst("P", conf);

        if (excluded_doms.size() == 0) {
            syst.add_def(def_drP.c_str());
            syst.add_def(def_drdrP.c_str());
        } else {
            auto last_dom = (bound_dom == -1) ? ndom : bound_dom;
            for (auto dom = 0; dom < last_dom; ++dom) {
            auto res = std::find(excluded_doms.begin(), excluded_doms.end(), dom);
            if (res == std::end(excluded_doms)) {
                syst.add_def(dom, def_drP.c_str());
                syst.add_def(dom, def_drdrP.c_str());
            }
            }
        }
        Scalar field(syst.give_val_def("ddrP"));
        field.std_base();
        return field;
    }
}