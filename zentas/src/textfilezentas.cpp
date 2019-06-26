// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>
#include <zentas/costs.hpp>
#include <zentas/fasta.hpp>
#include <zentas/zentas.hpp>
#include <zentas/zentaserror.hpp>

namespace nszen
{

void textfilezentas(std::vector<std::string> filenames,
                    std::string              outfilename,
                    std::string              costfilename,
                    size_t                   K,
                    std::string              algorithm,
                    size_t                   level,
                    size_t                   max_proposals,
                    bool                     capture_output,
                    std::string&             text,
                    size_t                   seed,
                    double                   max_time,
                    double                   min_mE,
                    double                   max_itok,
                    std::string              metric,
                    size_t                   nthreads,
                    size_t                   max_rounds,
                    bool                     patient,
                    std::string              energy,
                    bool                     with_tests,
                    bool                     rooted,
                    double                   critical_radius,
                    double                   exponent_coeff,
                    std::string              initialisation_method,
                    bool                     do_balance_labels)
{

  /* Input : filenames, outfilename,  costfilename
   *
   *
   * To populate :
   *
   * size_t ndata,
   * const size_t * const sizes,
   * const int * const ptr_datain,
   * const size_t * const indices_init
   * size_t * const indices_final,
   * size_t * const labels
   * bool with_cost_matrices,
   * size_t dict_size,
   * double c_indel,
   * double c_switch,
   * const double * const c_indel_arr,
   * const double * const c_switch_arr
   *
   * */

  if (initialisation_method == "from_indices_init")
  {
    throw zentas::zentas_error("cannonot initialise from indices (from_indices_init) in "
                               "textfilezentas (indices not clearly defined until data in array)");
  }

  std::mt19937_64                         gen(seed);
  std::uniform_int_distribution<unsigned> dis;

  /* (0) check that all specified files are valid */
  std::vector<std::string> all_filenames = filenames;
  all_filenames.push_back(costfilename);

  for (auto& fn : all_filenames)
  {
    std::ifstream input(fn);
    if (!input.good())
    {
      throw zentas::zentas_error("Error opening '" + fn + "'. ");
    }
  }

  std::ofstream ofs;
  ofs.open(outfilename, std::ofstream::out);
  ofs.close();

  /* (1) get ndata, sizes, ptr_datain from filenames */
  std::vector<char>        v_datain;
  std::vector<std::string> v_names;
  std::vector<size_t>      v_sizes;

  /* (1.01) determine if isfasta */
  bool isfasta = false;
  for (auto& fn : filenames)
  {
    isfasta = (isfasta || fasta::get_isfasta(fn));
  }

  if (isfasta == true)
  {
    for (auto& fn : filenames)
    {
      fasta::append_fasta(fn, v_datain, v_names, v_sizes);
    }
  }

  else
  {
    for (auto& fn : filenames)
    {
      txt::append_txt(fn, v_datain, v_sizes);
    }
  }

  /* (1.1) determine which letters are used (nucleotides or amino acids), convert to uppercase,
   * catch bad characters */
  std::set<char> chars_used;
  for (auto& c : v_datain)
  {
    // c = toupper(c);
    chars_used.insert(c);
  }
  const char* const   ptr_datain = v_datain.data();
  const size_t* const sizes      = v_sizes.data();
  size_t              ndata      = v_sizes.size();

  /* (2) get indices_init from seed, ndata, K */
  // std::vector<size_t> v_indices_init(K, 0);
  // if (K >= ndata/2){
  // throw zentas::zentas_error("Request for " + std::to_string(K) + " clusters with " +
  // std::to_string(ndata) + " samples rejected. Select a smaller K ( less than " +
  // std::to_string(ndata/2) + " ).");
  //}

  /* (3) initialise indices_final, labels */
  std::unique_ptr<size_t[]> uptr_indices_final(new size_t[K]);
  auto                      indices_final = uptr_indices_final.get();
  std::unique_ptr<size_t[]> uptr_labels(new size_t[ndata]);
  auto                      labels = uptr_labels.get();

  /* (4) get with_cost_matrices, dict_size, c_indel, c_switch, c_indel_arr, c_switch_arr from and
   * costfilename.
   * if with_cost_matrices, characters are mapped to the contiguous character set : {char(0), ...
   * char(n_distinct_input_chars)},
   * this allows the matrix c_switch_arr to be compact, resulting in (hopefully) quicker lookup.
   * variables map_into_contiguous and map_into_input manage the char <-> char mapping. I do not
   * completely understand chars, parts of this code may be stupid */
  std::map<char, char> map_into_contiguous;
  std::vector<char> map_into_input;

  /* default values, to be correctly set in a moment */
  double  c_indel            = 0.;
  double  c_switch           = 0.;
  size_t  dict_size          = 0;
  double* c_indel_arr        = nullptr;
  double* c_switch_arr       = nullptr;
  bool    with_cost_matrices = false;

  /* if with_cost_matrices is true, these vectors will be used and c_indel_arr and c_switch_arr
   * will point to them */
  std::vector<double> v_c_indel;
  std::vector<double> v_c_switch;

  /* set the cost variables, reading from costfilename */
  std::map<std::pair<char, char>, double> substitution_costs;
  std::map<char, double> indel_costs;
  costs::set_costs(costfilename, substitution_costs, indel_costs);

  /* case 1 : constant indel and switch costs */
  if (indel_costs.count('*') == 1)
  {
    c_indel            = indel_costs['*'];
    c_switch           = substitution_costs[std::pair<char, char>{'*', '*'}];
    with_cost_matrices = false;
  }

  /* case 2 : variable indel and switch costs */
  else
  {
    /* confirm that we have all the costs we need */
    for (auto& x : chars_used)
    {
      if (indel_costs.count(x) == 0)
      {
        throw zentas::zentas_error("It appears that the indel cost for char " + std::to_string(x) +
                                   " is not present in indel_costs. Conclusion : there is a "
                                   "missing line [char] [value] in the file " +
                                   costfilename);
      }
      for (auto& y : chars_used)
      {
        if (x != y)
        {
          if (substitution_costs.count(std::pair<char, char>{x, y}) == 0)
          {
            std::string error_message = "It appears that the substitution cost for  ";
            error_message = error_message + x + " " + y + " is not present in substitution_costs. "
                                                          "Conclusion : there is a missing line "
                                                          "[char] [char] [value] in the file " +
                            costfilename;
            throw zentas::zentas_error(error_message);
          }
        }
      }
    }

    dict_size = chars_used.size();

    /* set the char <-> char maps */
    int counter = 0;
    for (auto& x : chars_used)
    {
      map_into_contiguous[x] = static_cast<char>(counter);
      map_into_input.push_back(x);
      ++counter;
    }

    /* map data to contiguous character set */
    for (auto& c : v_datain)
    {
      c = map_into_contiguous[c];
    }

    /* set v_c_indel and v_c_switch */
    v_c_indel.resize(dict_size, 0);
    v_c_switch.resize(dict_size * dict_size, 0);
    for (size_t i = 0; i < dict_size; ++i)
    {
      v_c_indel[i] = indel_costs[map_into_input[i]];
      for (size_t j = 0; j < dict_size; ++j)
      {
        if (j != i)
        {
          v_c_switch[i * dict_size + j] =
            substitution_costs[std::pair<char, char>{map_into_input[i], map_into_input[j]}];
        }
      }
    }
    c_indel_arr        = v_c_indel.data();
    c_switch_arr       = v_c_switch.data();
    with_cost_matrices = true;
  }

  /* (5) call szentas */
  const size_t* const indices_init = nullptr;
  szentas(ndata,
          sizes,
          ptr_datain,
          K,
          indices_init,
          initialisation_method,
          algorithm,
          level,
          max_proposals,
          capture_output,
          text,
          seed,
          max_time,
          min_mE,
          max_itok,
          indices_final,
          labels,
          metric,
          nthreads,
          max_rounds,
          patient,
          energy,
          with_tests,
          rooted,
          with_cost_matrices,
          dict_size,
          c_indel,
          c_switch,
          c_indel_arr,
          c_switch_arr,
          critical_radius,
          exponent_coeff,
          do_balance_labels);

  /* (6) write results to outfilename */
  if (with_cost_matrices == true)
  {
    /* (6.1) convert from contiguous characters back to to input characters */
    for (auto& c : v_datain)
    {
      c = map_into_input[c];
    }
  }

  /* (6.2) write output to file */
  ofs.open(outfilename, std::ofstream::out);
  size_t i_start = 0;

  if (isfasta == true)
  {
    for (size_t i = 0; i < ndata; ++i)
    {
      ofs << v_names[i] << "  (" << labels[i] << ")\n";
      for (size_t d = 0; d < sizes[i]; ++d)
      {
        ofs << v_datain[i_start + d];
      }
      ofs << "\n";
      i_start += sizes[i];
    }
  }

  else
  {
    for (size_t i = 0; i < ndata; ++i)
    {
      for (size_t d = 0; d < sizes[i]; ++d)
      {
        ofs << v_datain[i_start + d];
      }
      ofs << " \t : " << labels[i] << "\n";
      i_start += sizes[i];
    }
  }
  ofs.close();
}
}
