// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <iostream>
#include <string>

#include <zentas/outputwriter.hpp>
#include <zentas/zentaserror.hpp>

namespace zentas
{
namespace outputwriting
{

OutputWriter::OutputWriter(bool to_terminal, bool to_file, std::string filename)
  : to_terminal(to_terminal), to_file(to_file), filename(filename)
{
  if (to_file == true)
  {
    if (filename.compare("") == 0)
    {
      throw zentas::zentas_error(
        "empty filename passed to OutputWrite, with to_file flag true. This is not possible.");
    }

    file.open(filename, std::ios::out);

    // std::ofstream file(filename, std::ios::out);
    if (file.good() == false)
    {
      std::string errm =
        "bad filename in constructor of OutputWriter object. The filename provided is `";
      errm += filename;
      errm += "'. The directory of the file must exist, OutputWriters do not create directories. "
              "Either create all directories in the path, or change the provided path.  ";
      throw zentas::zentas_error(errm);
    }
  }
}

// make default destructor.

OutputWriter::~OutputWriter() { file.close(); }

template <>
OutputWriter& OutputWriter::operator<<(Flusher f)
{
  f.increment();

  if (to_terminal)
  {

#ifndef COMPILE_FOR_R
    std::cout << std::flush;
#else
    Rcpp::Rout << std::flush;
#endif
  }

  if (to_file)
  {
    file.flush();
  }
  return *this;
}

template <>
OutputWriter& OutputWriter::operator<<(Endline e)
{
  e.increment();
  if (to_terminal)
  {
#ifndef COMPILE_FOR_R
    std::cout << std::endl;
#else
    Rcpp::Rout << std::endl;
#endif
  }

  if (to_file)
  {
    file << "\n";
    file.flush();
  }
  return *this;
}
}

outputwriting::Endline Endl;
outputwriting::Flusher Flush;
}
