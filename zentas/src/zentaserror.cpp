// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <zentas/stringutilbase.hpp>
#include <zentas/zentaserror.hpp>

#include <iostream>
#include <sstream>

namespace zentas
{

std::string tgformat(const std::string& what_arg, std::string prefix, std::string suffix)
{
  std::stringstream fms;
  fms << "\n\n";
  unsigned l_terminal = 95;

  auto frags = stringutil::split(what_arg, "\n");

  std::string space;
  space.resize(2 * l_terminal, ' ');
  std::string interspace("   ");

  unsigned line_current = 1;
  for (auto& x : frags)
  {
    if (x.size() == 0)
    {
      x = " ";
    }

    unsigned p_current = 0;
    while (p_current < x.size())
    {
      auto l_current = std::min<unsigned>(l_terminal, x.size() - p_current);
      fms << prefix << interspace << x.substr(p_current, l_current)
          << space.substr(0, l_terminal - l_current) << interspace << suffix << " ("
          << line_current << ")\n";
      p_current += l_current;
      ++line_current;
    }
  }

  fms << "\n";

  return fms.str();
}

zentas_error::zentas_error(const std::string& what_arg)
  : std::runtime_error(tgformat(what_arg, "zentas", "ERROR"))
{
}

void zentas_warning(const std::string& warning)
{
  std::cerr << tgformat(warning, "zentas", "WARNING") << std::flush;
}
}
