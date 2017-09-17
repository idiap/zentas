// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTASERROR_HPP
#define ZENTASERROR_HPP

#include <stdexcept>

namespace zentas
{

class zentas_error : public std::runtime_error
{
  public:
  zentas_error(const std::string& what_arg);
  //~zentas_error() = default;
};

void zentas_warning(const std::string& warning);
}

#endif
