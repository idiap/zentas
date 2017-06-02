#ifndef ZENTASINFO_HPP
#define ZENTASINFO_HPP

#include <sstream>
namespace nszen{

std::string get_output_info_string();
std::string get_python_paramater_string();
std::string get_python_function_decl();
std::string get_equals_line(size_t n);

}

#endif
