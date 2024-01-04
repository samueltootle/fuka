#include <memory>
namespace Kadath::FUKA_Solvers {
template<class FUKA_READER>
struct Interpolator {
    protected:
    std::unique_ptr<FUKA_READER> solution{nullptr};

    public:
    Interpolator(FUKA_READER const & _sol) {
        solution.reset(new FUKA_READER{_sol});
    }

    // DEBUG ONLY
    void print_field_coefs(size_t FIELD_IDX) {
        auto& f = solution->get_field_vector(FIELD_IDX);
        
        uint cnt = 0;
        for(const auto & e : f[2]) {
            std::cout << e << ", ";
            if(cnt > solution->get_ncoefs_i(2, 0)) {
                cout << '\n';
                cnt = 0;
            }
            else cnt++;
        }
        cout << '\n';
    }
};
}