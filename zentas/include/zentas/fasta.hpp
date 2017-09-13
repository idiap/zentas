// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_FASTA_HPP
#define ZENTAS_FASTA_HPP

#include <map>
#include <string>
#include <vector>

namespace txt
{

void append_txt(std::string filename, std::vector<char>& data, std::vector<size_t>& sizes);
}

namespace fasta
{

void append_fasta(std::string               filename,
                  std::vector<char>&        data,
                  std::vector<std::string>& names,
                  std::vector<size_t>&      sizes);

bool get_isfasta(std::string filename);
}

#endif
