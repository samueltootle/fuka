#pragma once
#include "coord_fields.hpp"
#include "Configurator/config_bco.hpp"
#include "Configurator/config_binary.hpp"
#include "Solvers/sequences/ns_sequence.hpp"

/**
 * \addtogroup Syst_tools
 * \ingroup FUKA
 * Centralized tools to populate and manipulate System_of_eqs
 * @{*/

namespace Kadath {
namespace FUKA_Syst_tools {
using namespace ::Kadath::FUKA_Config;

/**
 * @brief Initialize constants that are used for all initial data
 * 
 * @tparam cfields_t container type holding coordinate fields
 * @param syst System of equations to modify
 * @param coord_vectors container of coordinate fields
 */
template<class cfields_t>
void syst_init_csts(System_of_eqs & syst, cfields_t& coord_vectors) {
  #ifdef DEBUG
    std::cout << "Loading constants and constant fields.\n";
  #endif
  syst.add_cst("PI", M_PI);
  
  // coordinate fields
  syst.add_cst("mg", *coord_vectors[GLOBAL_ROT]);

  syst.add_cst("ex"  , *coord_vectors[EX])  ;
  syst.add_cst("ey"  , *coord_vectors[EY])  ;
  syst.add_cst("ez"  , *coord_vectors[EZ])  ;
  syst.add_cst("einf", *coord_vectors[S_INF]) ;
}

/**
 * @brief Initialize system definitions relevant for all initial data
 * 
 * @param syst System of equations to modify
 */
inline
void syst_init_defs(System_of_eqs & syst) {
  int const ndom = syst.get_space().get_nbr_domains();
  #ifdef DEBUG
    std::cout << "Loading global definitions.\n";
  #endif
  //Useful definitions
  syst.add_def("NP = P*N");
  syst.add_def("Ntilde = N / P^6");
  // Conformal Extrinsic curvature.
  syst.add_def("A^ij  = (D^i bet^j + D^j bet^i \
                      - 2. / 3.* D_k bet^k * f^ij) / 2. / Ntilde");

  // ADM Linear Momentum
  syst.add_def(ndom - 1, "intPx = A_i^j * ex_j * einf^i / 8 / PI") ;
  syst.add_def(ndom - 1, "intPy = A_i^j * ey_j * einf^i / 8 / PI") ;
  syst.add_def(ndom - 1, "intPz = A_i^j * ez_j * einf^i / 8 / PI") ;

  syst.add_def(ndom - 1, "intMadm = -dr(P) / 2 / PI");
  syst.add_def(ndom - 1, "intMk =  dr(N) / 4 / PI");

  // Irreducible mass integrand
  syst.add_def("intArea = P^4 / 4piG") ;
  syst.add_def("intMirrsq  = intArea / 4") ;
}

/**
 * @brief Initialize various tensor contractions for spacetime quantities
 * 
 * @param syst System of equations to modify
 */
inline
void syst_init_contraction_defs_vac(System_of_eqs & syst) {
  #ifdef DEBUG
    std::cout << "Loading vacuum contraction definitions.\n";
  #endif
  int const ndom = syst.get_space().get_nbr_domains();

  // Contractions of definitions and fields for analysis
  syst.add_def("Axx     = A^ij * ex_i * ex_j ");
  syst.add_def("Ayy     = A^ij * ey_i * ey_j");
  syst.add_def("Azz     = A^ij * ez_i * ez_j");
  syst.add_def("Axy     = A^ij * ex_i * ey_j ");
  syst.add_def("Ayz     = A^ij * ey_i * ez_j");
  syst.add_def("Axz     = A^ij * ex_i * ez_j");
  syst.add_def("TraceA  = A^ij * f_ij");

  // Contractions of definitions and fields for analysis
  syst.add_def("Bx = B^i  * ex_i ");
  syst.add_def("By = B^i  * ey_i");
  syst.add_def("Bz = B^i  * ez_i");

  syst.add_def("betx = bet^i  * ex_i ");
  syst.add_def("bety = bet^i  * ey_i");
  syst.add_def("betz = bet^i  * ez_i");
}

/**
 * @brief Add definitions specific to inspiral binaries
 * 
 * @tparam config_t Config type
 * @param syst System of equations to modify
 * @param bconfig Config file
 * @param CART coordinate field relative to the coordinate distance to the center of mass
 */
template<class config_t>
void syst_init_inspiral(System_of_eqs & syst, config_t& bconfig, Vector& CART) {
  int const ndom = syst.get_space().get_nbr_domains();
  std::string eccstr{};
  if(!std::isnan(bconfig.set(ADOT))) {
    syst.add_cst("adot", bconfig(ADOT));
    syst.add_cst("r", CART);
    syst.add_def("comr^i = r^i - xaxis * ex^i + yaxis * ey^i");
    eccstr +=" + adot * comr^i";
  } 

  // Define COM shifted orbital rotation field
  syst.add_def("Morb^i = mg^i + xaxis * ey^i + yaxis * ex^i");
  
  // Define global shift (inertial + corotating)
  std::string Bstr {"B^i= bet^i + ome * Morb^i" + eccstr};
  syst.add_def(Bstr.c_str());

  // ADM Angular Momentum Def
  syst.add_def(ndom - 1, "intJ = multr(A_ij * Morb^j * einf^i) / 2 / 4piG");
}

/**
 * @brief Initialize definitions important for two-body problems
 * 
 * @tparam cfields_t Coordinate fields container type
 * @tparam config_t Config type
 * @param syst System of equations to modify
 * @param coord_vectors Container of coordinate fields
 * @param bconfig Config file
 */
template<class cfields_t, class config_t>
void syst_init_binary(System_of_eqs & syst, 
  cfields_t& coord_vectors, config_t& bconfig) {
  #ifdef DEBUG
    std::cout << "Initializing standard binary fields, constants, and definitions.\n";
  #endif
  int const ndom = syst.get_space().get_nbr_domains();
  syst.add_cst("4piG", bconfig(QPIG)) ;
  syst_init_csts(syst, coord_vectors);
  syst_init_defs(syst);

  // remaining coordinate fields
  syst.add_cst("mm", *coord_vectors[BCO1_ROT]) ;
  syst.add_cst("mp", *coord_vectors[BCO2_ROT]) ;

  syst.add_cst("sm", *coord_vectors[S_BCO1])  ;
  syst.add_cst("sp", *coord_vectors[S_BCO2])  ;
}

/**
 * @brief Initialize quasi-local definitions based purely on spacetime
 * 
 * @param syst System of equations to modify
 * @param doms vector of domains to populate these definitions
 * @param rotdef String related to the rotation field definition
 * @param surfdef String related to the surface element definition
 */
inline
void syst_init_quasi_local_defs(System_of_eqs& syst,
  std::vector<int> doms, std::string rotdef, std::string surfdef) {
  #ifdef DEBUG
    std::cout << "Initializing QL spin definition.\n";
  #endif
  std::string qlspin_def = "intS = A_ij * " + rotdef + "^i * " + surfdef + "^j";
  qlspin_def += " / 2 / 4piG";
  for(auto& d : doms){
    syst.add_def(d, qlspin_def.c_str()) ;
  }
}

/**
 * @brief Initialize definitions relevent for isolated objects
 * 
 * @tparam cfields_t Coordinate fields container type
 * @tparam config_t Config type
 * @param syst System of equations to modify
 * @param coord_vectors Container of coordinate fields
 * @param bconfig Config file
 */
template<class cfields_t, class config_t>
void syst_init_co(System_of_eqs & syst, 
  cfields_t& coord_vectors, config_t& bconfig) {
  // FIXME boosts not considered
  #ifdef DEBUG
    std::cout << "Initializing standard CO fields, constants, and definitions.\n";
  #endif
  int const ndom = syst.get_space().get_nbr_domains();
  syst.add_cst("4piG", bconfig(BCO_QPIG));
  syst_init_csts(syst, coord_vectors);

  syst.add_cst("sm", *coord_vectors[S_BCO1])  ;

  syst.add_def("B^i = bet^i + ome * mg^i");
  syst_init_defs(syst);

  // ADM Angular Momentum Def
  syst.add_def(ndom - 1, "intJ = multr(A_ij * mg^j * einf^i) / 2. / 4piG");

  // Local spin definition
  syst.add_def("intS = A_ij * mg^i * sm^j / 2 / 4piG") ;
}

/**
 * @brief Initialize vacuum equation definitions
 * 
 * @param syst System of equations to modify
 * @param doms vector of domains to populate these definitions
 */
inline
void syst_init_eqdefs_vac(System_of_eqs& syst, std::vector<int> doms) {
  #ifdef DEBUG
    std::cout << "Initializing vacuum constraint equation definitions.\n";
  #endif
  for(auto& d : doms) {
    syst.add_def(d,"eqP     = D^i D_i P + A_ij * A^ij / P^7 / 8") ;
    syst.add_def(d,"eqNP    = D^i D_i NP - 7. / 8. * NP / P^8 * A_ij * A^ij");
    syst.add_def(d,"eqbet^i = D_j D^j bet^i \
                          + D^i D_j bet^j / 3. - 2. * A^ij * D_j Ntilde");
  }
}

template<class config_t>
inline std::string set_ns_mass_fixing(System_of_eqs& syst, config_t& bconfig, std::unique_ptr<Kadath::FUKA_Solvers::ns_sequence const>& seq) {
    std::string central_fixing_definition{"h - hc"};
    auto idx{seq->mass_idx()};
    switch(idx) {
      case BCO_PARAMS::HC:
        syst.add_cst("hc"  , bconfig(BCO_PARAMS::HC));
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        break;
      case BCO_PARAMS::NC:
        syst.add_cst("Nc", bconfig(BCO_PARAMS::NC));
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        central_fixing_definition = "rho - Nc";
        break;
      case BCO_PARAMS::MADM:
        syst.add_var("hc", bconfig(BCO_PARAMS::HC));
        syst.add_var("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_cst("Madm", bconfig(BCO_PARAMS::MADM));
        break;
      case BCO_PARAMS::MB:
        syst.add_var("hc", bconfig(BCO_PARAMS::HC));
        syst.add_cst("Mb"  , bconfig(BCO_PARAMS::MB));
        syst.add_var("Madm", bconfig(BCO_PARAMS::MADM));
        break;
      default:
        std::string msg{"Sequence initialized, but not implemented for Mass index = " + std::to_string(int(idx))};
        throw std::runtime_error(msg.c_str());
        break;
    }
    return central_fixing_definition;
}

template<class config_t>
inline std::string set_ns_spin_fixing(System_of_eqs& syst, config_t& bconfig, std::unique_ptr<Kadath::FUKA_Solvers::ns_sequence const>& seq) {
    std::string spin_fixing_definition{"integ(intJ) - chi * Madm * Madm = 0"};
    auto idx{seq->spin_idx()};
    switch(idx) {
      case BCO_PARAMS::OMEGA:
        syst.add_var("chi" , bconfig(BCO_PARAMS::CHI));
        syst.add_cst("ome" , bconfig(BCO_PARAMS::OMEGA));
        break;
      case BCO_PARAMS::JADM:
        spin_fixing_definition = "integ(intJ) - Jadm = 0";
        syst.add_cst("Jadm", bconfig(BCO_PARAMS::JADM));
        syst.add_var("ome" , bconfig(BCO_PARAMS::OMEGA));
        bconfig.set(BCO_PARAMS::CHI) = bconfig(BCO_PARAMS::JADM) / bconfig(BCO_PARAMS::MADM) / bconfig(BCO_PARAMS::MADM);
        break;
      case BCO_PARAMS::CHI:
        syst.add_cst("chi" , bconfig(BCO_PARAMS::CHI));
        syst.add_var("ome" , bconfig(BCO_PARAMS::OMEGA));
        break;
      default:
        std::string msg{"Sequence initialized, but not implemented for Spin index = " + std::to_string(int(idx))};
        throw std::runtime_error(msg.c_str());
        break;
    }
    return spin_fixing_definition;
}
/** @}*/
}}
#include "fuka_syst_setup_hydro.hpp"