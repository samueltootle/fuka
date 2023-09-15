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

#include "headcpp.hpp"
#include "name_tools.hpp"
#include "tensor.hpp"
#include <string>
#include <algorithm>
namespace Kadath {
// Removes the spaces in excess from a string
void trim_spaces (char* dest, const char* name) {
	int in_output = 0 ;
	int in_input = 0 ;
	bool inspace = (name[0] == ' ') ? true : false ;
	bool isspace ;

	// Get rid of the spaces
	char auxi[LMAX] ;
	while (name[in_input]!='\0') {
		isspace = (name[in_input]==' ') ? true : false ;
		if ((isspace==false) || ((isspace==true) && (inspace==false))) {
			auxi[in_output] = name[in_input] ;
			in_output ++ ;
		}
		inspace = isspace ;
		in_input++ ;
	}

	if (auxi[in_output-1]!= ' ') {
		auxi[in_output] = ' ' ;
		in_output ++ ;
		}
	auxi[in_output] = '\0' ;

	//Get rid of the parenthesis not needed :	  
	int first = 0 ;
	int last = in_output-2 ;
	bool doloop = true ; 

	while (doloop)  {

	  bool needtodo = (auxi[first]== '(') ? true : false ;
	  if (auxi[last]!=')')
		needtodo = false ;
	  if (needtodo) {
	    int ninside = 0 ;
	    for (int i=first ; i<last ; i++) {
		if (auxi[i]== '(')
		  ninside ++ ;
		if (auxi[i]== ')')
		    ninside-- ;
		if ((auxi[i]== ')') && (ninside==0))
		  needtodo = false ;
	    }

	    if ((needtodo) && (auxi[last] ==')')) {
		first ++ ;
		last -- ;
	    }
	  }
	    doloop = needtodo ;
	}

	int pos = 0 ;
	for (int i=first ; i<in_output-first-1 ; i++) {
		dest[pos] = auxi[i] ;
		pos ++ ;
	}
 
	dest[pos] = ' ' ;
	pos ++ ;
	dest[pos] = '\0' ;
}

void get_util (char* res, char* input) {

	char auxi[LMAX] ;
	trim_spaces(auxi, input) ;

	int first = 0 ;
	int conte = 0 ;
	while ((first==0) && (auxi[conte]!='\0')) {
		if (input[conte]=='_')
			first = 1 ;
		if (input[conte]=='^')
			first = 2 ;
		conte ++ ;
	}
	
	if (first==0) {
		strcpy(res, auxi) ;
	}
	else if (first==1) {
		get_term (auxi, res, '_') ;
		}
			else get_term(auxi, res, '^') ;
}

// Our many given char in a string ?
int nbr_char (const char* name, char xx) {

	const char* copie = name ;
	int conte = 0 ;
	while (*copie!='\0') {
		if (*copie==xx)
			conte ++ ;
		copie ++ ;
	}
	return conte ;
}

// Get some part of a string (up to the character sep)
void get_term (char* input, char* output, char sep) {

	char* found = strchr(input, sep) ;

	char auxi[LMAX] ;
	for (int i=0 ; i<int(strlen(input)-strlen(found)) ; i++)
		auxi[i] = input[i] ;
	auxi[int(strlen(input)-strlen(found))]='\0' ;

	trim_spaces(output, auxi) ;


	found++ ;
	int i=0 ;
	while (*found != '\0') {
		input[i] = *found ;
		i++ ;
		found++ ;
	}
	input[i]='\0' ;
}

void get_parts (const char* input, char* first, char* second, char sep, int place) {
	// Put ourselves in the right place :
	int length = static_cast<int>(strlen(input)) ;
	int nfound = -1 ;
	bool finloop = false ;
	int pos = length - 1 ;
	while (!finloop) {
		if (input[pos]==sep)
			nfound++ ;
		pos -- ;
		if (nfound==place)
			finloop = true ;
		if (pos==-1)
			finloop = true ;
	}
	if (pos==-1) {
	      // Not good !
	      first[0] = '\0' ;
	      second[0] = '0' ;
	}
	else {
	  pos ++ ;
	  char auxi1[LMAX] ;
	  for (int i=0 ; i<pos ; i++)
		auxi1[i] = input[i] ;
	  auxi1[pos]='\0' ;
	  trim_spaces(first, auxi1) ;

	  char auxi2[LMAX] ;
	  for (int i=pos+1 ; i<length ;  i++)
	      auxi2[i-pos-1] = input[i] ;
	  auxi2[length-1-pos]='\0' ;
	  trim_spaces(second, auxi2) ;
    }
}


bool is_tensor (const char* input, const char* name_tensor, int& valence, char*& name_ind, Array<int>*& type_ind) {
	
	bool res = (nbr_char(input, ' ')==1) ? true : false ;
	if (!res)
		return false ;
	else {
		int first = 0 ;
		int conte = 0 ;
		while ((first==0) && (input[conte]!='\0')) {
			if (input[conte]=='_')
				first = 1 ;
			if (input[conte]=='^')
				first = 2 ;
			conte ++ ;
		}

		char name[LMAX] ;
		char ind[LMAX] ;
	
		if (first==0) {
			// No indices 
			valence = 0 ;
			int same = strcmp(input, name_tensor);
			if (same!=0)
				return false;
			else
				return true ;
		}

		// get the name :
		int kant = 0 ;
		while ((input[kant]!='_') && (input[kant]!='^')) {
			name[kant] = input[kant] ;
			kant++ ;
		}
		name[kant] = ' ' ;
		name[kant+1] = '\0' ;

		int kkant = 0 ;
		while (input[kant]!='\0') {
			ind[kkant] = input[kant] ;
			kkant++ ;
			kant ++ ;
		}
		ind[kkant] = ' ' ;
		ind[kkant+1] = '\0' ;


		int same = strcmp(name, name_tensor);
		if (same!=0)
			return false ;
		else {
			valence = 0 ;
			int cc = 0 ;
			while (ind[cc]!='\0') {
				if ((ind[cc]!=' ') && (ind[cc]!='_') && (ind[cc]!='^'))
					valence ++ ;
				cc ++ ;
			}
			
			if (type_ind!=0x0)
			  delete type_ind ;
			type_ind = new Array<int> (valence) ;
			if (name_ind != 0x0)
			  delete [] name_ind ;
			name_ind = new char [valence] ;
			
			// On remplit les arrays
			cc = 0 ;
			int nind = 0 ;
			while (ind[cc]!='\0') {
				if ((ind[cc]!=' ') && (ind[cc]!='_') && (ind[cc]!='^')) {
					name_ind[nind] = ind[cc] ;
					if (first == 1)
						type_ind->set(nind) = COV ;
					else	
						type_ind->set(nind) = CON ;
					nind ++ ; 
				}
				if (ind[cc]=='_')
					first = 1 ;
				if (ind[cc]=='^')
					first = 2 ;
				cc ++ ;
			}
			return true ;
		}
	}
}
std::string extract_path(std::string fullvar) {
  return (fullvar.rfind("/") == std::string::npos) ? \
          "./" : fullvar.substr(0,fullvar.rfind("/")+1);
}

std::string extract_filename(std::string fullvar) {
  return fullvar.substr(fullvar.rfind("/")+1, fullvar.size());
}

// This is taken from cppreference
// https://en.cppreference.com/w/cpp/string/byte/tolower
std::string str_tolower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), 
                   [](unsigned char c){ return std::tolower(c); }
                  );
    return s;
}

}
