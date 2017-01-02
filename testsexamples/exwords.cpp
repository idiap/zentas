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
#include <vector>
#include <iostream>
#include <map>

#include "fasta.hpp"
#include "zentas.hpp"


/* Test case : word clustering. */
int cluster_words(){
  
  std::vector<std::string> filenames; 
  
  std::string root = PATH_DATADIR;
  
  //filenames.push_back(root + "nucleic.txt");
  filenames.push_back(root + "words1.txt");
  filenames.push_back(root + "words2.txt");
  std::string outfilename = root + "output.txt";  
  std::string  costfilename = root + "costs.txt";
  size_t K = 3;
  std::string algorithm("clarans");
  size_t level = 3;
  size_t max_proposals = 10000;
  size_t capture_output = false;
  std::string text;
  size_t seed = 1012;
  double maxtime = 10000.;
  std::string metric("normalised levenshtein");
  size_t nthreads = 1;
  size_t maxrounds = 100;
  bool patient = false;
  std::string energy("identity");
  bool rooted = true;
  double critical_radius = 0;
  double exponent_coeff = 0;
  nszen::textfilezentas(filenames, outfilename, costfilename, K, algorithm, level, max_proposals, capture_output, text, seed, maxtime, metric, nthreads, maxrounds, patient, energy, rooted, critical_radius, exponent_coeff);
 
  return 0; 
}

int main(){

  return cluster_words();
}
