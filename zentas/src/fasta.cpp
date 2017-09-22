// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>
#include <zentas/zentaserror.hpp>
#include <zentas/fasta.hpp>

namespace txt
{
void append_txt(std::string filename, std::vector<char>& data, std::vector<size_t>& sizes)
{

  std::ifstream input(filename);
  if (!input.good())
  {
    throw zentas::zentas_error("Error opening '" + filename + "'. ");
  }

  std::string line, content;

  while (std::getline(input, line).good())
  {

    // remove trailing white space
    line = std::regex_replace(line, std::regex(" +$"), "");

    if (line.empty() == true)
    {
    }

    else if (line[0] == '#')
    {  // A comment to be ignored
    }

    else if (line[0] == '>')
    {
      throw zentas::zentas_error("'>' is not a valid first character in a regular text file. This "
                                 "has been been done so as to prevent confusion between FASTA and "
                                 "regular txt files. If this is problematic, I can easily change "
                                 "it, please contact me : james . ne wling @ idi ap.ch");
    }

    else
    {
      data.insert(data.end(), line.begin(), line.end());
      sizes.push_back(line.size());
    }
  }
}
}

namespace fasta
{

enum State
{
  just_had_name,
  just_had_content,
  just_had_space,
  at_start
};

bool get_isfasta(std::string filename)
{

  bool          isfasta = false;
  std::ifstream input(filename);
  if (!input.good())
  {
    throw zentas::zentas_error("Error opening '" + filename + "'. ");
  }

  std::string line, name, content;
  while (std::getline(input, line).good())
  {
    if (line.empty() == false and line[0] == '>')
    {
      isfasta = true;
      break;
    }
  }

  return isfasta;
}

void append_fasta(std::string               filename,
                  std::vector<char>&        data,
                  std::vector<std::string>& names,
                  std::vector<size_t>&      sizes)
{

  std::ifstream input(filename);
  if (!input.good())
  {
    throw zentas::zentas_error("Error opening '" + filename + "'. ");
  }

  std::string line, name, content;

  State state = at_start;

  while (std::getline(input, line).good())
  {

    if (line.empty() == true)
    {
      if (state == just_had_name)
      {
        throw zentas::zentas_error("Empty lines may not follow name lines");
      }
      state = just_had_space;
    }

    else if (line[0] == '#')
    {  // A comment to be ignored
    }

    else if (line[0] == '>')
    {  // Identifier marker

      if (state == just_had_name)
      {
        throw zentas::zentas_error("Multiline names are not allowed");
      }

      if (line.size() == 1)
      {
        throw zentas::zentas_error("Empty names (>'') are not allowed");
      }

      if (content.size() != 0)
      {
        data.insert(data.end(), content.begin(), content.end());
        sizes.push_back(content.size());
        content.resize(0);
      }

      names.push_back(line);
      state = just_had_name;
    }

    else
    {

      if (state == just_had_space)
      {
        throw zentas::zentas_error("Blank lines are not allowed within content, or between name "
                                   "and content (only allowed between content and name). The "
                                   "current line is  ->[" +
                                   line + "]<-");
      }

      if (line.find(' ') != std::string::npos)
      {  // Invalid sequence--no spaces allowed
        throw zentas::zentas_error("A space detected in line : ->[" + line +
                                   "]<- . This line appears to be content, and thus should not "
                                   "contain any spaces (even at the end).");
      }

      content += line;
      state = just_had_content;
    }
  }

  if (content.size() == 0)
  {
    throw zentas::zentas_error("The file has been completed, expecting final content to append");
  }
  std::transform(content.begin(), content.end(), content.begin(), ::toupper);
  data.insert(data.end(), content.begin(), content.end());
  sizes.push_back(content.size());
}

}  // namespace
