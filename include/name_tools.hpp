/*
    Copyright 2017 Philippe Grandclement

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

#ifndef __NAME_TOOLS_HPP_
#define __NAME_TOOLS_HPP_

#include "headcpp.hpp"
#include "array.hpp"
#include <string>
#define LMAX 1000

namespace Kadath {
void trim_spaces (char* dest, const char* name) ;
void get_util (char*, char*) ;
int nbr_char (const char*, char) ;
void get_term (char*,char*,char) ;
void get_parts (const char*,char*,char*, char, int place = 0) ;
bool is_tensor (const char*, const char*, int&, char*&, Array<int>*&) ;
std::string extract_path(std::string fullvar);
std::string extract_filename(std::string fullvar);
std::string str_tolower(std::string s);
inline std::string stdio_header(std::string s, const char header_char = '*') {
    int n = ((42 - s.size()) > 0) ? 42 - s.size() : s.size() - 42;
    n /= 2;
    return std::string(n, header_char) + s + std::string(n, header_char);
}
}
#endif
