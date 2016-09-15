/*Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

zentas is a k-medoids library written in C++ and Python. This file is part of zentas.
zentas is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License version 3 as published by the Free Software Foundation.
zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of
the GNU General Public License along with zentas. If not, see
<http://www.gnu.org/licenses/>.
*/
#ifndef ZENTAS_FASTA_HPP
#define ZENTAS_FASTA_HPP


#include <string>
#include <map>


namespace txt{

void append_txt(
std::string filename, 
std::vector<char> & data, 
std::vector<size_t> & sizes);
}

namespace fasta{

void 
append_fasta(
std::string filename, 
std::vector<char> & data, 
std::vector<std::string> & names, 
std::vector<size_t> & sizes );


bool
get_isfasta(std::string filename);

}

#endif
