#pragma once
#include "coord_fields.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"

/**
 * \addtogroup Syst_tools
 * \ingroup FUKA
 * Centralized tools to populate and manipulate System_of_eqs
 * @{*/

namespace Kadath {
namespace FUKA_Syst_tools {
  
/**
 * @brief Add various tensor contractions of fluid quantities
 * 
 * @param syst System of equations to modify
 */
inline
void syst_init_contraction_defs_hydro(System_of_eqs & syst) {
  #ifdef DEBUG
    std::cout << "Loading hydro contractions.\n";
  #endif
  // useful constraction definitions
  syst.add_def("drhodx = ex^i * D_i rho");
  syst.add_def("dHdx = ex^i * D_i H");
  syst.add_def("dHdx2 = ex^i * D_i dHdx");
}

/**
 * @brief Initialize fluid operators and definitions needed
 * for solving or accessing initial data
 * 
 * @tparam eos_t EOS type
 * @param syst System of equations to modify
 */
template<class eos_t>
void syst_init_defs_hydro(System_of_eqs& syst) {
  #ifdef DEBUG
    std::cout << "Loading hydro operators and definitions.\n";
  #endif
  // user defined OPEs to access EOS
  Param p;
  syst.add_ope("eps"  , &EOS<eos_t,EPSILON>::action, &p);
  syst.add_ope("press", &EOS<eos_t,PRESSURE>::action, &p);
  syst.add_ope("rho"  , &EOS<eos_t,DENSITY>::action, &p);
  syst.add_ope("dHdlnrho", &EOS<eos_t,DHDRHO>::action, &p);

  syst.add_def("h = exp(H)") ;
  syst.add_def("press = press(h)");
  syst.add_def("eps = eps(h)");
  syst.add_def("rho = rho(h)");
  syst.add_def("dHdlnrho = dHdlnrho(h)");
  syst.add_def("delta = h - eps - 1.");
}

/**
 * @brief Initialize fluid definitions in the case of corotation since these are
 * significantly simpler than the arbitrary rotation case
 * 
 * @param syst System of equations to modify
 * @param doms vector of domains to populate these definitions
 */
inline
void syst_init_eqdefs_hydro_corot(System_of_eqs& syst, std::vector<int> doms) {
  #ifdef DEBUG
    std::cout << "Setting hydro corotating constraints defintions.\n";
  #endif
  for(auto& d : doms) {
    
    // Only valid for corotation
    syst.add_def(d, "U^i = B^i / N");
    syst.add_def(d, "Usquare = P^4 * U_i * U^i") ;

    syst.add_def(d, "Wsquare = 1 / (1 - Usquare)");
    syst.add_def(d, "W = sqrt(Wsquare)");
    // end corot

    // Source term defintions
    syst.add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
    syst.add_def(d, "Stilde = 3 * press * delta \
                            + (Etilde + press * delta) * Usquare") ;
    syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

    // Constraint equation definitions inside the star
    syst.add_def(d, "eqP  = delta * D^i D_i P \
                          + A_ij * A^ij / P^7 / 8 * delta \
                          + 4piG / 2. * P^5 * Etilde") ;
    syst.add_def(d, "eqNP = delta * D^i D_i NP \
                          - 7. / 8. * NP / P^8 * delta * A_ij *A^ij \
                          - 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
    syst.add_def(d, "eqbet^i = delta * D_j D^j bet^i \
                              + delta * D^i D_j bet^j / 3. \
                              - 2. * delta * A^ij * D_j Ntilde \
                              - 4. * 4piG * N * P^4 * ptilde^i");
    syst.add_def(d, "firstint = H + log(N) - log(W)");
  }
}

/**
 * @brief Initialize fluid definitions in the case of arbitrary rotation
 * 
 * @param syst System of equations to modify
 * @param doms vector of domains to populate these definitions
 * @param spin_def string definition for the spin field
 */
inline
void syst_init_eqdefs_hydro(System_of_eqs& syst, std::vector<int> doms,
  std::string spin_def) {
  #ifdef DEBUG
    std::cout << "Setting hydro constraint defintions\n";
  #endif
  for(auto& d : doms) {
    syst.add_def(d, spin_def.c_str());
    syst.add_def(d, "eta_i  = D_i phi + P^4 * s_i");
    
    syst.add_def(d, "Wsquare = eta^i * eta_i / h^2 / P^4 + 1.");
    syst.add_def(d, "W = sqrt(Wsquare)");

    syst.add_def(d, "U^i = eta^i / P^4 / h / W");
    syst.add_def(d, "Usquare= P^4 * U_i * U^i") ;
    
    syst.add_def(d, "V^i = N * U^i - B^i");    

    // First Integral
    syst.add_def(d, "firstint = log(h * N / W + D_i phi * V^i)") ;

    // velocity potential equations
    syst.add_def(d, "eqphi  = P^6 * W * V^i * D_i H \
                            + dHdlnrho * D_i (P^6 * W * V^i)");

    // Source term defintions
    syst.add_def(d, "Etilde = press * h * Wsquare - press * delta") ;
    syst.add_def(d, "Stilde = 3 * press * delta \
                            + (Etilde + press * delta) * Usquare") ;
    syst.add_def(d, "ptilde^i = press * h * Wsquare * U^i") ;

    // Constraint equation definitions inside the star
    syst.add_def(d, "eqP  = delta * D^i D_i P \
                          + A_ij * A^ij / P^7 / 8 * delta \
                          + 4piG / 2. * P^5 * Etilde") ;
    syst.add_def(d, "eqNP = delta * D^i D_i NP \
                          - 7. / 8. * NP / P^8 * delta * A_ij *A^ij \
                          - 4piG / 2. * N * P^5 * (Etilde + 2. * Stilde)");
    syst.add_def(d, "eqbet^i = delta * D_j D^j bet^i \
                              + delta * D^i D_j bet^j / 3. \
                              - 2. * delta * A^ij * D_j Ntilde \
                              - 4. * 4piG * N * P^4 * ptilde^i");
  }
}

/**
 * @brief Initialize quasi-local definition of fluid diagnostics
 * 
 * @param syst System of equations to modify
 * @param doms vector of domains to populate these definitions
 */
inline
void syst_init_quasi_local_defs_hydro(System_of_eqs& syst, 
  std::vector<int> doms) {
  #ifdef DEBUG
    std::cout << "Setting QL hydro definitions.\n";
  #endif
  for(auto& d : doms){
    // Baryonic Mass
    syst.add_def(d, "intMb = P^6 * rho * W") ;
    // QLMADM
    syst.add_def(d, "intqlMadm  = - D_i D^i P * 2. / 4piG") ;
    // Integral of log specific enthalpy
    syst.add_def(d, "intH  = P^6 * H * W") ;
  }
}
/** @}*/
}}
