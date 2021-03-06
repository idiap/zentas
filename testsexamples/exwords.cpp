// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <iostream>
#include <map>
#include <vector>
#include <zentas/fasta.hpp>
#include <zentas/zentas.hpp>

/* Test case : word clustering. */
/* many of the parameters used here are the same as in exdense.cpp, if in doubt look there*/
int cluster_words()
{

  std::vector<std::string> filenames;

  std::string root = PATH_DATADIR;

  // filenames.push_back(root + "nucleic.txt");
  filenames.push_back(root + "words1.txt");
  filenames.push_back(root + "words2.txt");
  std::string outfilename  = root + "output.txt";
  std::string costfilename = root + "costs.txt";
  size_t      K            = 3;
  std::string algorithm("clarans");
  size_t      level          = 3;
  size_t      max_proposals  = 10000;
  size_t      capture_output = false;
  std::string text;
  size_t      seed     = 1012;
  double      max_time = 10000.;
  double      min_mE   = 0.;
  double      max_itok = 100.0;
  std::string metric("normalised levenshtein");
  size_t      nthreads   = 1;
  size_t      max_rounds = 100;
  bool        patient    = false;
  std::string energy("identity");
  bool        rooted          = true;
  double      critical_radius = 0;
  double      exponent_coeff  = 0;
  std::string initialisation_method("kmeans++-10");
  bool        with_tests        = false;
  bool        do_balance_labels = false;
  nszen::textfilezentas(filenames,
                        outfilename,
                        costfilename,
                        K,
                        algorithm,
                        level,
                        max_proposals,
                        capture_output,
                        text,
                        seed,
                        max_time,
                        min_mE,
                        max_itok,
                        metric,
                        nthreads,
                        max_rounds,
                        patient,
                        energy,
                        with_tests,
                        rooted,
                        critical_radius,
                        exponent_coeff,
                        initialisation_method,
                        do_balance_labels);

  return 0;
}

int main() { return cluster_words(); }
