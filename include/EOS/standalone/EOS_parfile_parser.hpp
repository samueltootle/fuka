#pragma once
#include <array>
#include <cmath>
#include <string>
#include <fstream>
#include "cold_pwpoly.hh"
#include "margherita.hh"
#ifdef WITH_GRHAYL_EOS
#include <grhayl/ghl.h>
#define FUKA_MAX_EOS_PARAMS \
  MAX(MAX_EOS_PARAMS, Kadath::Margherita::Cold_PWPoly::max_num_pieces)
#else
#define FUKA_MAX_EOS_PARAMS Kadath::Margherita::Cold_PWPoly::max_num_pieces
#endif

namespace Kadath {
namespace FUKA_EOS {
class parse_polytrope_file {
 public:
  double rhomin;
  double rhomax;
  double Pmin;
  double K0;
  int num_pieces;
  std::array<double, FUKA_MAX_EOS_PARAMS> gamma_tab;
  std::array<double, FUKA_MAX_EOS_PARAMS> rho_tab;
  std::string Units;

 private:
  inline void rescale_input() {
    using namespace Margherita_constants;
    const double gam0m1 = gamma_tab[0] - 1.0;
    double rho_unit = 1.;
    double K_unit = 1.;
    if (Units != "geometrised") {
      if (Units == "cgs") {
        rho_unit = 1.0 * RHOGF;
        K_unit = pow(INVRHOGF, gam0m1) / c2_cgs;
      } else {
        if (Units == "cgs_cgs_over_c2") {
          rho_unit = 1.0 * RHOGF;
          K_unit = pow(INVRHOGF, gam0m1);
        } else {
          std::cerr << "Unit system, " << Units << ", not recognised!\n";
          std::_Exit(EXIT_FAILURE);
        }
      }
    }

    K0 *= K_unit;
    for (int i = 0; i < num_pieces; ++i)
      rho_tab[i] = rho_unit * rho_tab[i];
  }

  inline void parse_file(std::string fname) {
    std::ifstream f(fname);

    if (!f.is_open()) {
      std::string msg = "failed to open " + fname;
      throw std::runtime_error(msg.c_str());
    }

    //string descriptor to ignore
    std::string descr;

    //lambda to ignore leading comments and blank lines
    //up to the next value to extract
    auto skip_comments = [&]() {
      auto peek_c = f.peek();
      while (peek_c == '#' || peek_c == '\n') {
        std::getline(f, descr, '\n');
        peek_c = f.peek();
      }
    };

    auto skip_and_grab = [&](auto& val) {
      skip_comments();
      f >> descr >> val;
    };

    skip_and_grab(num_pieces);
    assert(num_pieces <= FUKA_MAX_EOS_PARAMS);

    skip_and_grab(rhomin);
    skip_and_grab(rhomax);
    skip_and_grab(K0);
    skip_and_grab(Pmin);

    auto read_tab = [&](auto& ary) {
      skip_comments();
      f >> descr;
      //Can't use foreach since array is static length
      for (int i = 0; i < num_pieces; ++i) {
        if (!(f >> ary[i])) {
          std::cerr << "Not enough vars in " << descr << " for " << num_pieces
                    << "pieces.\n";
          std::_Exit(EXIT_FAILURE);
        }
      }
    };

    read_tab(gamma_tab);
    read_tab(rho_tab);
    skip_and_grab(Units);
  }

 public:
  void operator()(std::string polytrope_file) {
    parse_file(polytrope_file);
    rescale_input();
  }
};

class parse_3d_EOS_par_file {
 public:
  double rho_atm;
  double rho_min;
  double rho_max;
  double Ye_atm;
  double Ye_min;
  double Ye_max;
  double T_atm;
  double T_min;
  double T_max;
  double T_beta;
  std::string table_abspath;
  std::string table_format;

 private:
  inline void parse_file(std::string fname) {
    std::ifstream f(fname);

    if (!f.is_open()) {
      std::string msg = "failed to open " + fname;
      throw std::runtime_error(msg.c_str());
    }

    //string descriptor to ignore
    std::string descr;

    //lambda to ignore leading comments and blank lines
    //up to the next value to extract
    auto skip_comments = [&]() {
      auto peek_c = f.peek();
      while (peek_c == '#' || peek_c == '\n') {
        std::getline(f, descr, '\n');
        peek_c = f.peek();
      }
    };

    auto skip_and_grab = [&](auto& val) {
      skip_comments();
      f >> descr >> val;
    };

    skip_and_grab(table_abspath);
    skip_and_grab(table_format);

    // Grab atmospheric values
    skip_and_grab(rho_atm);
    skip_and_grab(Ye_atm);
    skip_and_grab(T_atm);

    // Grab min and max values
    skip_and_grab(rho_min);
    skip_and_grab(rho_max);
    skip_and_grab(Ye_min);
    skip_and_grab(Ye_max);
    skip_and_grab(T_min);
    skip_and_grab(T_max);
    skip_and_grab(T_beta);
  }

 public:
  void operator()(std::string table_par_file) { parse_file(table_par_file); }
};
}  // namespace FUKA_EOS
}  // namespace Kadath