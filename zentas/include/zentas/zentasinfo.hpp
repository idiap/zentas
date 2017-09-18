// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTASINFO_HPP
#define ZENTASINFO_HPP

#include <sstream>
#include <map>
namespace nszen
{


std::map<std::string, std::string> init_rf_dict();
const std::map<std::string, std::string> & get_rf_dict();

std::map<std::string, std::string> init_seq_dict();
const std::map<std::string, std::string> & get_seq_dict();

std::string get_equals_line(size_t);
std::string get_char_line(size_t n, char c);
std::string get_python_init_string();
std::string get_python_outdict_string();
std::string get_output_verbose_string();
std::string get_python_txt_seq_string();
std::string get_python_seq_string();
std::string get_python_spa_string();
std::string get_python_den_string();
}

#endif
