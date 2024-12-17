/*
 * =====================================================================================
 *
 *       Filename:  1dtov.cpp
 *
 *    Description:  Simple TOV test code to test tov.hh
 *
 *        Version:  1.0
 *        Created:  1/07/2021
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel David Tootle, tootle@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *   Notes: this is a minimal test for the tov.hh code
 *
 * =====================================================================================
 */


#define PWPOLY_SETUP
#define DEBUG
#include "cold_table.hh"
#include "cold_table_implementation.hh"
#include "cold_pwpoly.hh"
#include "cold_pwpoly_implementation.hh"
#include "setup_cold_table.cc"
#include "setup_polytrope.cc"
//#include "/home/user/Downloads/konrad/tov.hh"
#include "tov.hh"
#include "../FUKA_EOS_Wrapper.hh"
#include <memory>

int main() {
  using namespace Kadath::Margherita;
  using namespace Kadath::FUKA_EOS;
  // Margherita_setup_polytrope("gam2.polytrope");
  setup_Cold_Table("togashi.lorene",1000);
  // auto tov = std::make_unique<MargheritaTOV<Cold_PWPoly>>();
  auto tov = std::make_unique<MargheritaTOV<FUKA_EOS_Wrapper<margherita_eos_t, margherita_1d>>>();
//  tov->adaptive = false;
//  tov->rk45 = false;
  //tov->solve(1.37e-3);
  tov->solve_for_MADM(1.4003505615);
  auto& state = tov->state;
  std::cout << *tov << std::endl;
  return 0;
}
