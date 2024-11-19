/*
 * =====================================================================================
 *
 *       Filename:  tov.hh
 *
 *    Description:  Margherita TOV: Simple TOV solving routines
 *
 *        Version:  2.0
 *        Created:  1/07/2021
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel David Tootle, tootle@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *   Notes: this is a rewrite of the original TOV solver by ERM based on the needs
 *   of FUKA.
 *
 * =====================================================================================
 */
#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <functional>

namespace Kadath {
namespace Margherita {
template <typename EOS> class MargheritaTOV {
  public:
  enum {RADIUS=0, RHOB, PRESS, MASSR, MASSB, PHI, RISO, CONF, ENTHALPY, TIDALH, TIDALBETA, NUMVAR};
  typedef EOS eos_t;

  private:
  using ary_t = std::array<double,NUMVAR>;
  using state_t = std::vector<ary_t>;
  
  static constexpr size_t max_iter = 40000;
  static constexpr double dr0 = 10. / 10000.;
  const double loop_eps = 1e-6;
  const double stop_at_fraction_press_c = 1e-13;

  public:
  bool adaptive{true};
  bool rk45{true};
  bool compact_output{false};
  double radius{0};
  double arealr{0};
  double mass{0};
  double rhoc{0};
  double baryon_mass{0};
  double tidal_love_k2{0};
  double press_c{0};
  state_t state{};

  void reset(const double& rho_init) {
    state.clear();
    state.reserve(max_iter);
    radius = 0;
    arealr = 0;
    mass = 0;
    rhoc = rho_init;
    baryon_mass = 0;
    tidal_love_k2 = 0;
    press_c = 0;
  }
  
  /**
   * gen_initial
   *
   * here the initial conditions are set.  Update this function
   * in the event more variables are added.
   */
  inline ary_t gen_initial() {
    using namespace Kadath::Margherita;
    double eps;
    typename EOS::error_t err;
    press_c = EOS::press_cold_eps_cold__rho(eps, rhoc, err);
    ary_t initial_conditions {};

    initial_conditions[RHOB] = rhoc;
    initial_conditions[PRESS] = press_c;
    initial_conditions[PHI] = 1.;
    return initial_conditions;
  }

  /** 
   * evolve
   *
   * update function for dy's - returns k * h
   * if more variables are added, one needs to add the relevant
   * equations here and update the class enum.
   *
   * @param f your input array of dy's
   * @param h step size
   */
  inline ary_t evolve(const ary_t& f, const double& h) const {
    typename EOS::error_t err;

    ary_t res{};
    res.fill(0);
    auto m = f[MASSR];
    double r = f[RADIUS];
    double r2 = r * r;
    double r3 = r2* r;
    double dedp, press, rho, rhoE;
    press = f[PRESS];
    rho = EOS::rho_energy_dedp__press_cold(rhoE, dedp, press, err);
    double sch_fac, mass_taylor;

    auto H = f[TIDALH];
    auto bet = f[TIDALBETA];
    if(state.size() == 1) {  
      H = r2;
      bet = 2. * r;
    }
    // eq(178) Tichy. https://arxiv.org/abs/1610.03805
    res[MASSR] = 4. * M_PI * r2 * rhoE * h;

    // dividing into two regimes based on: https://github.com/zachetienne/nrpytutorial/blob/master/Tutorial-ADM_Initial_Data-TOV.ipynb
    if(r<1e-4 || m<=0){
    
      res[PHI] = (4.*M_PI/3. * rhoE * r + 4. * M_PI * r * press) / (1. - 8. * M_PI * rhoE *r2) * h;
      res[RISO] = 1./std::sqrt(1. - 8. * M_PI * rhoE * r2) * h;
      sch_fac = 1./(1. - 8. * M_PI/3. * r2 * rhoE);
      
      double mass_taylor_by_r = 4. * M_PI/3. * r2 * rhoE;
      double mass_taylor_by_r2 = 4. * M_PI/3. * r * rhoE;

      auto tmpvar = mass_taylor_by_r2 + 4. * M_PI * r * press;
      res[TIDALBETA] = 2. * sch_fac * (
        -2. * M_PI * (5. * rhoE + 9. * press + dedp * (rhoE + press)) * r2
        + 3. + 2. * sch_fac * tmpvar * tmpvar)
        + 2. * (2.0) * sch_fac * (-1. + mass_taylor_by_r+ 2. * M_PI * r2 * (rhoE - press));

    }
    else{
      // eq(180) Tichy. https://arxiv.org/abs/1610.03805
      res[PHI] = (m + 4. * M_PI * r3 * press) / (r * r * (1. - 2 * m/r)) * h;
      res[RISO] = std::sqrt(1./(1. - 2.*m/r)) * f[RISO] / r  * h;
      sch_fac = r / (r - 2. * m) ;
      mass_taylor = m;
      auto tmpvar = mass_taylor / r2 + 4. * M_PI * r * press;
      res[TIDALBETA] = 2. * sch_fac * H * (
        -2. * M_PI * (5. * rhoE + 9. * press + dedp * (rhoE + press))
        + 3. / r2 + 2. * sch_fac * tmpvar * tmpvar)
        + 2. * bet / r * sch_fac * (-1. + mass_taylor / r + 2. * M_PI * r2 * (rhoE - press));
    }

    // eq(179) Tichy. https://arxiv.org/abs/1610.03805
    res[PRESS] = -res[PHI] * (rhoE + press);
    res[MASSB] = std::sqrt(sch_fac) * 4. * M_PI * r2 * rho * h;

    // eq (11,12) https://arxiv.org/pdf/0911.3535.pdf
    res[TIDALH] = bet * h;
    res[TIDALBETA] *= h;
    
    return res;
  }
  
  // RK45 with option for adaptive step sizes - enabled by default
  inline ary_t rk45_step(const ary_t& input, double& dr) const { 
    using namespace Kadath::Margherita;
    typename EOS::error_t err;
    auto eps = 5.e-8 * input[PRESS];

    auto k1 = evolve(input, dr);

    auto s = input;
    
    s[RADIUS] += (2./9.) * dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (2./9.) * k1[i];
    }
    
    auto k2 = evolve(s, dr);
    
    s = input;
    s[RADIUS] += dr / 3.;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += 1./12. * k1[i] + 0.25 * k2[i] ;
    }

    auto k3 = evolve(s, dr);
    
    s = input;
    s[RADIUS] += 0.75 * dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (69./128.) * k1[i] + (-243./128.) * k2[i] + (135./64.) * k3[i];
    } 

    auto k4 = evolve(s, dr);
    
    s = input;
    s[RADIUS] += dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (-17./12.) * k1[i] + (27./4.) * k2[i] + (-27./5.) * k3[i] + (16./15.) * k4[i];
    } 
    auto k5 = evolve(s, dr);

    s = input;
    s[RADIUS] += (5./6.) * dr ;
    for(auto i = 1; i < NUMVAR; ++i) {
      s[i] += (65./432.) * k1[i] + (-5./16.) * k2[i] + (13./16.) * k3[i] + (4./27.) * k4[i] + (5./144.) * k5[i];
    } 
    auto k6 = evolve(s, dr);
    
    auto res = input;
    for(int i = 1; i < NUMVAR; ++i)
      res[i] += (1./450.) * (47. * k1[i] + 216. * k3[i] + 64. * k4[i] + 15. * k5[i] + 108. * k6[i]);
    
    // manual update of tracked variables
    res[RADIUS] += dr;
    res[RHOB] = EOS::rho__press_cold(res[PRESS],err);
    
    ary_t error;
    for (int nn = 1; nn < NUMVAR; ++nn) {
      error[nn] = 1. / 300. *
                  (-2. * k1[nn] + 9. * k3[nn] - 64. * k4[nn] - 15. * k5[nn] +
                   72. * k6[nn]);
    };

    auto adaptive_step = [&]() {
      // Is the error too large?
      if (std::abs(error[PRESS]) > eps && state.size() > 2) {
        if (dr > 2. * dr0) {
          dr = dr * 0.5;
        } else {
          dr *= 0.9 * std::pow(eps / std::abs(error[PRESS]), 1./5.);
        }
        res = rk45_step(input, dr);
      } else if (std::abs(error[PRESS]) < 1. / 20. * eps) {
        dr *= 2.;
      }
    };
    if(adaptive)
      adaptive_step();
       
    return res;
  }
  
  // basic RK4 - mainly for testing
  inline ary_t rk_step(const ary_t& input, double& dr) const { 
    using namespace Kadath::Margherita;
    typename EOS::error_t err;
    // auto eps = 5.e-8 * input[PRESS];

    auto k1 = evolve(input, dr);

    auto s1 = input;
    
    s1[RADIUS] += dr / 2.;
    for(auto i = 1; i < NUMVAR; ++i)
      s1[i] += k1[i] / 2.;
    
    auto k2 = evolve(s1, dr);
    
    auto s2 = input;
    s2[RADIUS] += dr / 2.;
    for(auto i = 1; i < NUMVAR; ++i)
      s2[i] += k2[i] / 2.;

    auto k3 = evolve(s2, dr);
    
    auto s3 = input;
    s3[RADIUS] += dr;
    for(auto i = 1; i < NUMVAR; ++i)
      s3[i] += k3[i];
    
    auto k4 = evolve(s3, dr);
    
    auto res = input;
    for(int i = 1; i < NUMVAR; ++i)
      res[i] += 1./6. * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    
    // manual update of tracked variables
    res[RADIUS] += dr;
    res[RHOB] = EOS::rho__press_cold(res[PRESS],err);
   
    return res;
  }

  // correct PHI based on enforcing Schwarzschild BC
  void correct_phi() {
    // eq. (5) Baumgarte Notes
    auto get_correction = [&](double m, double r, double phi) {
      return 0.5 * std::log(1. - 2 * m / r) - phi;
    };
    
    const double correction = 
      get_correction(state.back()[MASSR], state.back()[RADIUS], state.back()[PHI]);

    for(auto& el : state) 
      el[PHI] += correction;
  }

  /**
   * correct_isotropic_r
   * calculate the final isotropic radius and correct it based
   * on enforcing the Schwarzschild BC
   */
  inline void correct_isotropic_r() {

     const double R_Schw=state.back()[RADIUS];
     const double M_Schw=state.back()[MASSR];
     const double RISO_Schw=state.back()[RISO];

     // based on the notebook of Etienne:
    for(auto& el : state)  
      el[RISO] *= (1./2.) * (std::sqrt(R_Schw * ( R_Schw - 2.0 * M_Schw)) + R_Schw - M_Schw)/ RISO_Schw;
  }

  /**
   * calculate_conformal_factor
   *
   * calculate conformal factor CONF such that
   * \gamma_ij = CONF^4 \eta_ij
   * eq(22) from Baumgarte notes
   */
  inline void calculate_conformal_factor() {
    for(auto& el : state) {
      const auto& r = el[RADIUS];
      const auto& rhat = el[RISO];
      el[CONF] = std::exp(0.5 * std::log(r / rhat));
    }
  }

  /**
   * calculate dimensionless love number k2
   *
   * eq(14) https://arxiv.org/pdf/0911.3535.pdf
   */
  inline void calculate_k2() {
    const auto& sol = state.back();
    const double y = sol[RADIUS] * sol[TIDALBETA] / sol[TIDALH];
    const double C = sol[MASSR] / sol[RADIUS];
    const double C2 = C*C;
    const double C3 = C2 * C;
    const double C5 = C3 * C2;
    const double tmp = (1. - 2. * C);
    const double tmp2 = tmp * tmp;
    tidal_love_k2 = 8. * C5 / 5. * tmp2 * (
        2. + 2. * C * (y - 1.) - y)
      / (2. * C * (6. - 3. * y + 3. * C * (5. * y - 8.))
      + 4. * C3 * (13. - 11. * y + C * (3. * y - 2.) + 2. * C2 * (1. + y))
      + 3. * tmp2 * (2. - y + 2. * C * (y - 1.)) * std::log(tmp));
  }


  /**
   * calculate enthalpy (e.g. to determine the actual areal radius by log(h)=0)
   *
   */
  void fillout_enthalpy(){
    for (auto &el: state ){
      typename EOS::error_t err;
      double rhoE, dedp;
      double press = el[PRESS];
      double rho = EOS::rho_energy_dedp__press_cold(rhoE, dedp, press, err);
      double eps;
      typename EOS::error_t err2,press_c;
      double rhoc=rho;
      press_c = EOS::press_cold_eps_cold__rho(eps, rhoc, err);
      el[ENTHALPY]=1+eps + press/rho ;
    }
  }

  /**
   * determine areal radius based on the h=1)
   *
   */
  void determine_ArealR_from_enthalpy(){
    #ifdef TOV1D_DEBUG
    double ArealR;
    for (auto &el: state ){
        if(el[ENTHALPY]<1){
          ArealR=el[RADIUS];
          break;
        }
    }
    std::cout << ArealR << std::endl;
    #endif
  }

  void determine_ArealR_from_integral(){
    #ifdef TOV1D_DEBUG
    double ArealR, IsotropicR, ConfAtR,Conf2AtR;//,Psi4AtR;
    for (auto &el: state ){

        if(el[ENTHALPY]<1){
            IsotropicR=el[RISO];
            ConfAtR=el[CONF];
            Conf2AtR=ConfAtR*ConfAtR;
            ArealR=IsotropicR*Conf2AtR;
            break;
        }
    }
    std::cout << ArealR << std::endl;
    #endif
  }

  /**
   * solve
   *
   * find a solution to the TOV equations for the given central and EOS
   *
   * @param rho_init central density input
   * @param eps precision desired for finding ADM mass
   */
  inline void solve(const double rho_init, const bool fin = true) { 
    using std::placeholders::_1;
    using std::placeholders::_2;
    // determines which RK function to use instead of doing the
    // logic check multiple times in the loop.
    std::function<ary_t(ary_t&, double&)> rk_iteration = (rk45) ? 
      std::bind(&MargheritaTOV::rk45_step, this, _1, _2) : 
      std::bind(&MargheritaTOV::rk_step, this, _1, _2);
    
    // reset all member variables to defaults
    reset(rho_init);

    // initialize our starting dr in case adaptive RK is used
    double dr = dr0;

    // add initial conditions
    state.push_back(gen_initial());
    
    size_t count = 0;
    while(true && count < max_iter) {
      auto res = rk_iteration(state.back(), dr);
      
      typename EOS::error_t err;
      double dedp, press, rho, rhoE;
      press = state.back()[PRESS];
      rho = EOS::rho_energy_dedp__press_cold(rhoE, dedp, press, err);
      
      // similar condition to the Etienne's TOV solver 
      // (https://github.com/zachetienne/nrpytutorial/blob/master/Tutorial-ADM_Initial_Data-TOV.ipynb)
      if(press < press_c*stop_at_fraction_press_c) 
        break; 
      state.push_back(res);
      count++;
    }
    if(count >= max_iter) {
      std::stringstream ss;
      ss << "Maximum iterations (" << max_iter << ") reached.\n"
        "The central density provided, " << rho_init << ", may not be able to produce a stable\n"
        "static TOV solution\n\n";
      throw std::runtime_error(ss.str().c_str());
    }

    // erase initial conditions if a stable solution was found at all
    // if this is neglected, it causes issues in calculating RISO since
    // you get division by 0
    if(count > 0)
      state.erase(state.begin());
    else
      return;
   
    // save some computation if we just need ArealR, M, and Mb
    if(fin) {
      correct_phi();
      correct_isotropic_r();
      calculate_conformal_factor();
      radius = state.back()[RISO];
      calculate_k2();
      fillout_enthalpy();
      determine_ArealR_from_enthalpy();
      determine_ArealR_from_integral();
    }
    arealr = state.back()[RADIUS];
    mass = state.back()[MASSR];
    baryon_mass = state.back()[MASSB];
  }
  
  /**
   * solve_for_MADM
   *
   * Root finder that takes an Madm value and desired precision and attempts
   * to find a solution to the TOV equations for the given Madm and EOS by
   * varying the central density
   *
   * @param M_fin desired final ADM mass
   * @param rho_max0 initial maximum density to search for MADM
   * @param rho_min0 initial minimum density to search for MADM
   * @param eps precision desired for finding ADM mass
   */
  inline bool solve_for_MADM(double M_fin, const double rho_max0=1e-2, 
    const double rho_min0=5e-4, const double eps = 1e-3) {
    // if M_fin < Mmax let the calling code know in case
    // the calling code intends to solve a rotating NS with M > Mtov
    bool use_Mmax = false; 

    // rho to evaluate
    double rho_eval=0.;

    // initial bracketing range
    double rho_min = rho_min0;
    double rho_max = rho_max0;

    // deltarho to find bracket around M_fin
    double delrho = (rho_max - rho_min) / 200.;

    double maxM = 0.;
    double maxMrho = 0.;
  
    // find bracketing range based on the central density to 
    // to make sure a root exists - i.e. Madm can be found
    auto rho_bracketing = [&]() {
      maxM = 0;
      maxMrho = 0;
      rho_min = 0;
      rho_max = 0;
      for(rho_eval = rho_min0; rho_eval < rho_max0; rho_eval += delrho) {
        solve(rho_eval, false);
        if(mass > maxM) {
          maxM = mass;
          maxMrho = rho_eval;
        }

        if(mass < M_fin)
          rho_min = rho_eval;
        else if(mass > M_fin) {
          rho_max = rho_eval;
          break;
        }
      }
      // std::cout << rho_min << ", " << rho_max << std::endl;
    };
    rho_bracketing();
   
    if(rho_eval >= rho_max0) {
      std::cerr << "Maximum density (" << rho_eval << ") reached.\n"
        << "Mass mass obtained: " << maxM << ", with density: "<< maxMrho << ".\n"
        "The density bracketing range may not be sufficiently constrained for the given EOS\n";
      if(maxM > 1.) {
        use_Mmax = true;
        solve(0.99*maxMrho, true);
        return use_Mmax;
        //M_fin = maxM * 0.95;
        //rho_bracketing();
      }
      else
        std::_Exit(EXIT_FAILURE);
    }

    // from the discovered bracketing range, we now find the ADM Mass using
    // a basic bisection of the central density
    size_t count = 0;
    double diff = 1;
    while( std::fabs(diff) > eps && count < max_iter) {
      rho_eval = (rho_max + rho_min) / 2.;      
      solve(rho_eval, false);
      diff = (mass-M_fin);
      
      #ifdef margerita_check_all
      std::cout << "M: " << mass
                << "\nMb: " << baryon_mass
                << "\nArealR: " << arealr
                << "\nNc: " << rho_eval
                << "\nDiff: " << diff 
                << "\nIter: " << count << "\n\n";
      #endif
      if( diff > 0) 
        rho_max = rho_eval;
      else
        rho_min = rho_eval;
      count++;
    }
    if(count >= max_iter) {
      std::cerr << "Maximum iterations (" << max_iter << ") reached .\n"
        "The mass may be too high or the range of rho may not be sufficiently wide\n";
      std::_Exit(EXIT_FAILURE);
    }
      
    correct_phi();
    correct_isotropic_r();
    calculate_conformal_factor();
    calculate_k2();
    fillout_enthalpy();
    radius = state.back()[RISO];
    
    return use_Mmax;
  }

  friend std::ostream& operator<< (std::ostream& stream, const MargheritaTOV& tov) {

    if(tov.compact_output){
    stream << std::setprecision(12);
    stream << tov.press_c << "\t"
           << tov.mass << "\t"
           << tov.baryon_mass << "\t"
           << tov.radius << "\t"
           << tov.tidal_love_k2 <<std::endl;
    }else{
    stream << std::setprecision(12);
    stream << " Central pressure: " << tov.press_c << std::endl;
    stream << " Central density: " << tov.rhoc << std::endl;
    stream << " Mass: " << tov.mass << std::endl;
    stream << " Baryon mass: " << tov.baryon_mass << std::endl;
    stream << " Radius: " << tov.arealr << std::endl;
    stream << " Isotropic Radius: " << tov.radius << std::endl;
    stream << " Love number k2: " << tov.tidal_love_k2 << std::endl;
    };

    return stream;
  }
};
}}
