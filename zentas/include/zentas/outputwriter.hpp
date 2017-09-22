// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include <fstream>
#include <iostream>
#include <string>

/* here is where you should define COMPILE_FOR_R.
 * TODO : the Rcpp trick does not work with templates.
 * mowri should jst be preprocessor defined.
 * ....
 *  */

namespace zentas
{
namespace outputwriting
{

class Flusher
{
  public:
  void increment(){}
};

class Endline
{
  public:
  void increment(){}
};

class OutputWriter
{

  public:
  bool          to_terminal;
  bool          to_file;
  std::ofstream file;
  std::string   filename;

  OutputWriter();
  ~OutputWriter();
  OutputWriter(bool to_terminal, bool to_file, std::string filename = "");
  void operator()(std::string);

  template <typename T>
  OutputWriter& operator<<(T t)
  {

    if (to_terminal)
    {

#ifndef COMPILE_FOR_R
      std::cout << t;
#else
      Rcpp::Rout << t;
#endif
    }

    if (to_file)
    {
      file << t;
    }

    return *this;
  }
};

template <>
OutputWriter& OutputWriter::operator<<(Flusher f);

template <>
OutputWriter& OutputWriter::operator<<(Endline e);
}

extern outputwriting::Flusher Flush;
extern outputwriting::Endline Endl;
}

#endif
