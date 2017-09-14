// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <atomic>
#include <cmath>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <vector>
#include <zentas/levenshtein.hpp>
#include <zentas/zentaserror.hpp>

namespace nszen
{

LevenshteinInitializer::LevenshteinInitializer(const size_t        dict_size_,
                                               const double        c_indel_,
                                               const double        c_switch_,
                                               const double* const c_indel_arr_,
                                               const double* const c_switch_arr_,
                                               bool                normalised_)
  : dict_size(dict_size_),
    c_indel(c_indel_),
    c_switch(c_switch_),
    c_indel_arr(c_indel_arr_),
    c_switch_arr(c_switch_arr_),
    normalised(normalised_)
{

  bool do_tests = true;

  if (do_tests && c_switch_arr != nullptr)
  {

    std::stringstream ss;
    bool              is_fine = true;

    /* normalised Levenshtein has special constrains on indel */
    if (do_tests && normalised == true)
    {
      for (unsigned i = 0; i < dict_size; ++i)
      {
        if (c_indel_arr[i] != c_indel_arr[0])
        {
          ss << "The normalised Levenshtein `metric' is only a true metric when the indel costs "
                "are all the same (see abstract of Yujian et al. 2007). ";
          ss << "However, cost_indel[" << i << "] is not the same as cost_indel[0]";
          ss << " (" << c_indel_arr[i] << ") vs (" << c_indel_arr[0] << "). ";
          ss << "As an example of why this constraint is required, consider c_switch = 10, "
                "c_indel(1) = 9, c_indel(2) = 10, ";
          ss << "and the sequences S1 = 12  S2 = 121 S3 = 21. The d(S1, S2) + d(S2, S3) = 36/59 "
                "and d(S1, S3) = 36/58 : not a metric!\n\n";
          is_fine = false;
        }
      }
    }

    /* confirm symmetry */
    for (unsigned i = 0; i < dict_size; ++i)
    {
      for (unsigned j = 0; j < dict_size; ++j)
      {
        if (c_switch_arr[i * dict_size + j] != c_switch_arr[j * dict_size + i])
        {
          ss << "cost_switch is not symmetric, it should be ";
          ss << "(" << c_switch_arr[i * dict_size + j] << " != " << c_switch_arr[j * dict_size + i]
             << ")\n\n";
          is_fine = false;
        }
        if (c_switch_arr[i * dict_size + j] < 0)
        {
          ss << "cost_switch should contain no negative values, it does. \n\n";
          is_fine = false;
        }
      }
    }

    if (is_fine == false)
    {
      throw zentas::zentas_error(ss.str());
    }

    /* confirm triangle inequality holds. may be a quicker way to do this. */
    /* i to j via k */

    for (unsigned i = 0; i < dict_size; ++i)
    {
      for (unsigned j = 0; j < dict_size; ++j)
      {
        double d_ij = c_switch_arr[dict_size * i + j];
        /* (1) switch matrix (checking switch(a,b) + switch(b,c) >= switch(a,c))*/
        for (unsigned k = 0; k < dict_size; ++k)
        {
          double d_ik = c_switch_arr[dict_size * i + k];
          double d_jk = c_switch_arr[dict_size * j + k];
          if (d_ij > d_ik + d_jk)
          {
            ss << "d[" << i << "][" << j << "]"
               << " > "
               << "d[" << i << "][" << k << "]"
               << " + "
               << "d[" << j << "][" << k << "], ";
            ss << d_ij << " > " << d_ik << " + " << d_jk << ".  \n\n";
            is_fine = false;
          }
        }
        /* (2) i to j with 2 indels must be more expensive than switch from i to j */
        if (d_ij > c_indel_arr[i] + c_indel_arr[j])
        {
          ss << "switch[" << i << "][" << j << "]"
             << " > "
             << "indel[" << i << "]"
             << " + "
             << "indel[" << j << "], ";
          ss << d_ij << " > " << c_indel_arr[i] << " + " << c_indel_arr[j] << ". \n\n";
          is_fine = false;
        }
        /* (3) finally, switch a to b then indel b must be more expensive than indel a */
        if (d_ij + c_indel_arr[i] < c_indel_arr[j])
        {
          ss << "switch[" << i << "][" << j << "] + indel[" << i << "]"
             << " < "
             << "indel[" << j << "], ";
          ss << d_ij << " + " << c_indel_arr[i] << " < " << c_indel_arr[j] << " .\n\n";
          is_fine = false;
        }
      }
    }

    if (is_fine == false)
    {
      std::string indel_broken =
        "the cost_switch/cost_indel matrices do not obey the triangle inequality.  ";
      indel_broken +=
        "Consider `reducing' the `inflated' costs. See Wagner and Fischer, 1974 for details. ";
      throw zentas::zentas_error(indel_broken + ss.str());
    }
  }
}

LevenshteinInitializer::LevenshteinInitializer(const double c_indel_,
                                               const double c_switch_,
                                               bool         normalised_)
  : LevenshteinInitializer(0, c_indel_, c_switch_, nullptr, nullptr, normalised_)
{
}

LevenshteinInitializer::LevenshteinInitializer(const size_t        dict_size_,
                                               const double* const c_indel_arr_,
                                               const double* const c_switch_arr_,
                                               bool                normalised_)
  : LevenshteinInitializer(dict_size_, 0., 0., c_indel_arr_, c_switch_arr_, normalised_)
{
}

LevenshteinInitializer::LevenshteinInitializer()
  : LevenshteinInitializer(0, 0., 0., nullptr, nullptr, false)
{
}

class MultiIndel
{
  private:
  const double* values;

  public:
  MultiIndel(const double* const vs) : values(vs) {}
  MultiIndel() {}

  double operator()(size_t i) const { return values[i]; }
};

/* dict_size x dict_size switch costs for Levenshtein */
class MultiSwitch
{
  private:
  const double* values;
  const size_t  dict_size;

  public:
  MultiSwitch(const double* const vals, size_t ds) : values(vals), dict_size(ds)
  {
  }

  double operator()(size_t i, size_t j) const { return values[i * dict_size + j]; }
};

class ConstIndel
{
  private:
  double value;

  public:
  ConstIndel(double val) : value(val) {}
  ConstIndel() {}

  double operator()(size_t) const { return value; }
};

class ConstSwitch
{
  private:
  double value;

  public:
  ConstSwitch(double v) : value(v) {}

  double operator()(int i, int j) const { return value * (i != j); }
};

/* Random comment. The use of min_c_indel might not be the optimal strategy in Levenshtein */
template <class TSample, class TIndelCost, class TSwitchCost>
void set_levenshtein_distance(const TSample&    v_vertical,
                              const TSample&    v_horizontal,
                              double            threshold,
                              size_t            dict_size,
                              const TIndelCost  f_c_indel,
                              double            min_c_indel,
                              double            max_c_indel,
                              const TSwitchCost f_c_switch,
                              double*           A_prev,
                              double*           A_acti,
                              int               nrows,
                              int               ncols,
                              size_t&           n_cells_visited_local,
                              double&           distance)
{

  (void)dict_size;

  int length_difference = ncols - nrows;
  if (length_difference < 0)
  {
    throw zentas::zentas_error("ncols must not be smaller than nrows in "
                               "set_levenshtein_distnance, ie shorter sequence first");
  }

  threshold = std::min(threshold, max_c_indel * (nrows + ncols));
  if (min_c_indel * length_difference >= threshold)
  {
    distance = threshold;
  }

  else
  {
    int j_start;
    int j_end;
    int half_width = std::min(ncols - 1, static_cast<int>(std::floor(threshold / min_c_indel)));
    int width      = 3 + 2 * half_width;
    std::fill(A_prev, A_prev + width, threshold);
    std::fill(A_acti, A_acti + width, threshold);

    A_acti[1 + half_width] = 0;
    for (int j = 2 + half_width; j < width; ++j)
    {
      A_acti[j] = A_acti[j - 1] + f_c_indel(v_horizontal.values[j - 2 - half_width]);
    }
    int    row          = 0;
    int    column       = 0;
    double min_distance = 0;
    while (row < nrows && min_distance < threshold)
    {

      std::swap(A_acti, A_prev);
      ++row;
      min_distance = threshold;
      j_start      = std::max(1, 1 + half_width - row);
      j_end        = std::min(ncols - row + 2 + half_width, width - 1);
      for (int j = j_start; j < j_end; ++j)
      {
        column = row + j - (1 + half_width);

        A_acti[j] = std::min(
          std::min(f_c_indel(v_horizontal.values[column - 1]) + A_acti[j - 1],
                   f_c_indel(v_vertical.values[row - 1]) + A_prev[j + 1]),
          A_prev[j] + f_c_switch(v_vertical.values[row - 1], v_horizontal.values[column - 1]));

        min_distance = std::min(min_distance, A_acti[j]);
      }
      n_cells_visited_local += j_end - j_start;
    }

    if (std::abs(nrows - ncols) <= half_width + 1)
    {
      distance = std::min(threshold, A_acti[1 + half_width + ncols - nrows]);
    }
    else
    {
      distance = threshold;
    }
  }
}

template <typename TSDataIn>
/* See levenshtein_prototyping for equivalent python version */
class LevenshteinMetric_X
{

  private:
  std::vector<size_t> v_ncalcs;
  const size_t        dict_size;
  bool                normalised;

  ConstIndel  cv_c_indel;
  ConstSwitch cv_c_switch;
  MultiIndel  av_c_indel;
  MultiSwitch av_c_switch;

  std::vector<size_t>                    v_n_cells_visited;
  std::vector<size_t>                    v_n_cells_visitable;
  int                                    max_size;
  int                                    memory_size;
  size_t                                 nthreads;
  std::vector<std::mutex>                v_mutex0;
  std::mutex                             counter_mutex;
  std::vector<std::unique_ptr<double[]>> v_up_A_acti;
  std::vector<double*>                   v_A_acti;
  std::vector<std::unique_ptr<double[]>> v_up_A_prev;
  std::vector<double*>                   v_A_prev;

  double min_c_indel;
  double max_c_indel;

  bool test_where_possible;

  protected:
  std::mt19937_64                         gen;
  std::uniform_int_distribution<unsigned> dis;

  public:
  typedef typename TSDataIn::Sample Sample;
  typedef LevenshteinInitializer    Initializer;

  LevenshteinMetric_X(const TSDataIn& datain, size_t nthreads_, const LevenshteinInitializer& li)
    : v_ncalcs(nthreads_, 0),

      dict_size(li.dict_size),
      normalised(li.normalised),

      cv_c_indel(li.c_indel),
      cv_c_switch(li.c_switch),
      av_c_indel(li.c_indel_arr),
      av_c_switch(li.c_switch_arr, dict_size),

      v_n_cells_visited(nthreads_, 0),
      v_n_cells_visitable(nthreads_, 0),
      max_size(datain.get_max_size()),

      /* below : making 10*size + 10 makes no difference to performance (speed) */
      memory_size(4 * datain.get_max_size() + 10),
      nthreads(nthreads_),
      v_mutex0(nthreads_)

  {
    for (size_t ti = 0; ti < nthreads; ++ti)
    {
      v_up_A_acti.emplace_back(new double[memory_size]);
      v_A_acti.push_back(v_up_A_acti[ti].get());

      v_up_A_prev.emplace_back(new double[memory_size]);
      v_A_prev.push_back(v_up_A_prev[ti].get());
    }

    if (dict_size > 0)
    {
      min_c_indel = std::numeric_limits<double>::max();
      max_c_indel = 0.;
      for (size_t w = 0; w < dict_size; ++w)
      {
        min_c_indel = std::min(min_c_indel, av_c_indel(w));
        max_c_indel = std::max(max_c_indel, av_c_indel(w));
      }
    }

    else
    {
      min_c_indel = cv_c_indel(0);
      max_c_indel = cv_c_indel(0);
    }

    test_where_possible = false;
    if (test_where_possible == true)
    {
      std::cerr << "\n\nLEVENSHTEIN WITH (LIMITED) TESTS ENABLED : WILL BE SLOWER" << std::endl;
    }
  }

  void set_distance_simple_test(const Sample& v_vertical,
                                const Sample& v_horizontal,
                                double        threshold,
                                double&       distance)
  {
    (void)threshold;
    distance = std::abs(int(v_vertical.size) - int(v_horizontal.size));
  }

  void set_distance(const Sample& v_vertical,
                    const Sample& v_horizontal,
                    double        threshold,
                    double&       distance)
  {

    /* make sure the shorter vector comes first */
    if (v_vertical.size < v_horizontal.size)
    {
      set_distance_tiffany(v_vertical, v_horizontal, threshold, distance);
    }
    else
    {
      set_distance_tiffany(v_horizontal, v_vertical, threshold, distance);
    }

    if (test_where_possible && v_vertical.size == v_horizontal.size)
    {
      double d1;
      double d2;
      set_distance_tiffany(v_vertical, v_horizontal, threshold, d1);
      set_distance_tiffany(v_horizontal, v_vertical, threshold, d2);

      if (d1 != d2)
      {
        std::stringstream ss;
        ss << "Haha! Levenshtein distances not the same when order reversed  \n";
        ss << v_vertical.str() << " ->  " << v_horizontal.str() << " : " << d1 << "\n";
        ss << v_horizontal.str() << " ->  " << v_vertical.str() << " : " << d2 << "\n";
        throw zentas::zentas_error(ss.str());
      }
    }
  }

  void set_distance_tiffany(const Sample& v_vertical,
                            const Sample& v_horizontal,
                            double        threshold,
                            double&       distance)
  {

    /* numerical issues */
    threshold *= 1.0000230507000110130001701900023;

    /* n_d = 2d / ( alpha (L1 + L2) + d )  where alpha = max indel cost
     * => d = alpha n_d (L1 + L2) / (2 - n_d)
     * */
    if (normalised == true)
    {
      // threshold in normalised space

      threshold = std::min(1., threshold);

      // threshold in non-normalised space, where the distance calculation is to take place.
      threshold = max_c_indel * threshold *
                  static_cast<double>(v_vertical.size + v_horizontal.size) / (2. - threshold);
    }

    int    nrows                 = static_cast<int>(v_vertical.size);
    int    ncols                 = static_cast<int>(v_horizontal.size);
    size_t n_cells_visited_local = 0;

    // get a mutex and keep it
    size_t mutex_i = dis(gen) % nthreads;

    while (v_mutex0[mutex_i].try_lock() == false)
    {
      ++mutex_i;
      mutex_i %= nthreads;
    }

    std::lock_guard<std::mutex> lock(v_mutex0[mutex_i], std::adopt_lock);
    double*                     A_acti = v_A_acti[mutex_i];
    double*                     A_prev = v_A_prev[mutex_i];

    /* do some fast tracking tests */
    if (nrows == 0 || ncols == 0)
    {
      throw zentas::zentas_error(
        "empty string, I need to confirm that this is not a special case. remind me to do this !");
    }

    else
    {
      /* constant indel and switch cost */
      if (dict_size == 0)
      {
        set_levenshtein_distance(v_vertical,
                                 v_horizontal,
                                 threshold,
                                 dict_size,
                                 cv_c_indel,
                                 min_c_indel,
                                 max_c_indel,
                                 cv_c_switch,
                                 A_prev,
                                 A_acti,
                                 nrows,
                                 ncols,
                                 n_cells_visited_local,
                                 distance);
      }

      /* matrices of costs */
      else
      {
        set_levenshtein_distance(v_vertical,
                                 v_horizontal,
                                 threshold,
                                 dict_size,
                                 av_c_indel,
                                 min_c_indel,
                                 max_c_indel,
                                 av_c_switch,
                                 A_prev,
                                 A_acti,
                                 nrows,
                                 ncols,
                                 n_cells_visited_local,
                                 distance);
      }
    }

    if (normalised == true)
    {
      /* return to the normalised space */
      distance =
        2. * distance /
        (max_c_indel * static_cast<double>(v_horizontal.size + v_vertical.size) + distance);
    }

    ++v_ncalcs[mutex_i];
    v_n_cells_visitable[mutex_i] += nrows * ncols;
    v_n_cells_visited[mutex_i] += n_cells_visited_local;

    distance = static_cast<double>(static_cast<float>(distance));
  }

  size_t get_ncalcs() const
  {
    size_t ncalcs = 0;
    for (auto& x : v_ncalcs)
    {
      ncalcs += x;
    }
    return ncalcs;
  }

  double get_rel_calccosts() const
  {

    /* How well have we done as compared to the O(row*column) algorithm ?*/
    size_t n_cells_visited = 0;
    for (auto& x : v_n_cells_visited)
    {
      n_cells_visited += x;
    }

    size_t n_cells_visitable = 0;
    for (auto& x : v_n_cells_visitable)
    {
      n_cells_visitable += x;
    }

    return static_cast<double>(n_cells_visited) / static_cast<double>(n_cells_visitable);
  }
};

template <typename TSDataIn>
LevenshteinMetric<TSDataIn>::LevenshteinMetric(const TSDataIn&               datain,
                                               size_t                        nthreads,
                                               const LevenshteinInitializer& li)
{
  lvsm.reset(new LevenshteinMetric_X<TSDataIn>(datain, nthreads, li));
}

template <typename TSDataIn>
void LevenshteinMetric<TSDataIn>::set_distance(const Sample& v_vertical,
                                               const Sample& v_horizontal,
                                               double        threshold,
                                               double&       distance)
{
  lvsm->set_distance(v_vertical, v_horizontal, threshold, distance);
}

template <typename TSDataIn>
size_t LevenshteinMetric<TSDataIn>::get_ncalcs()
{
  return lvsm->get_ncalcs();
}

template <typename TSDataIn>
double LevenshteinMetric<TSDataIn>::get_rel_calccosts()
{
  return lvsm->get_rel_calccosts();
}

template <typename TSDataIn>
LevenshteinMetric<TSDataIn>::~LevenshteinMetric() = default;

template class LevenshteinMetric<nszen::StringDataUnrootedIn<char>>;
template class LevenshteinMetric<nszen::StringDataUnrootedIn<int>>;
template class LevenshteinMetric<nszen::StringDataRootedIn<char>>;
template class LevenshteinMetric<nszen::StringDataRootedIn<int>>;
}
