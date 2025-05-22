/*
    Copyright 2019 Samuel Tootle

    This file is part of Kadath.

    Kadath is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Kadath is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Kadath.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <Configurator/config_bco.hpp>
#include <Configurator/config_binary.hpp>
namespace Kadath {
namespace FUKA_Config {
std::ostream & operator<<(std::ostream &out, const BCO_INFO &BCO) {
  std::string s = BCO.get_type() + " info";
  int n = ((42 - s.size()) > 0) ? 42 - s.size() : s.size() - 42;
  n /= 2;
  std::string title = std::string(n, '*') + s + std::string(n, '*');
  out << title << std::endl;
  Kadath::FUKA_Config_Utils::print_params(BCO.bco_map, BCO.bco_params, out);
  if(auto ns_ptr = dynamic_cast<const BCO_NS_INFO*>(&BCO)) {
      print_params(ns_ptr->get_eos_map(), ns_ptr->return_eos_params());
      if(!is_storage_all_nan(ns_ptr->return_diffrot_params())) {
        s = "differential rotation parameters";
        n = ((42 - s.size()) > 0) ? 42 - s.size() : s.size() - 42;
        n /= 2;
        title = std::string(n, '*') + s + std::string(n, '*');
        out << title << std::endl;
        print_params(ns_ptr->get_diffrot_map(), ns_ptr->return_diffrot_params());
      }
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, const BIN_INFO &BIN) {
  std::string s = BIN.get_type() + " info";
  int n = ((42 - s.size()) > 0) ? 42 - s.size() : s.size() - 42;
  n /= 2;
  std::string title = std::string(n, '*') + s + std::string(n, '*');
  out << title << std::endl;
  print_params(BIN.bin_map, BIN.bin_params);
  out << *BIN.BCOS[0];
  out << *BIN.BCOS[1];
  return out;
}
}}