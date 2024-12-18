#include<EOS/FUKA_EOS_Utilities.hh>
#include<kadath.hpp>
using namespace Kadath::FUKA_EOS;
/**
 * @brief Setup the EOS operators in System_of_eqs
 *
 * @tparam eos_t FUKA_EOS_Wrapper type
 * @param syst System of equations
 * @param p Parameter container (unused, but required by add_ope)
 */
template <class eos_t>
void set_eos_ope_struct<eos_t>::operator()(System_of_eqs& syst, Param& p) {
    syst.add_ope("eps", &EOS<eos_t, EPSILON>::action, &p);
    syst.add_ope("press", &EOS<eos_t, PRESSURE>::action, &p);
    syst.add_ope("rho", &EOS<eos_t, DENSITY>::action, &p);
}