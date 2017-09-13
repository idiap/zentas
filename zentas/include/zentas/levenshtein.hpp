// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_LEVENSTHTEIN_HPP
#define ZENTAS_LEVENSTHTEIN_HPP

#include <zentas/tdatain.hpp>

namespace nszen
{

class LevenshteinInitializer
{

  public:
  size_t dict_size;
  double c_indel;
  double c_switch;
  /*  will use either the above (if dict_size == 0) or the below (otherwise). */
  const double* c_indel_arr;
  const double* c_switch_arr;
  bool          normalised;

  LevenshteinInitializer(const size_t        dict_size,
                         const double        c_indel,
                         const double        c_switch,
                         const double* const c_indel_arr,
                         const double* const c_switch_arr,
                         bool                normalised);
  LevenshteinInitializer(const double c_indel, const double c_switch, bool normalised);
  LevenshteinInitializer(const size_t        dict_size,
                         const double* const c_indel_arr,
                         const double* const c_switch_arr,
                         bool                normalised);
  LevenshteinInitializer();
};

/* using pimpl for Levenshtein. This is the working class, */
template <typename TSDataIn>
class LevenshteinMetric_X;

/* using pimpl for Levenshtein. This is the wrapping class, */
template <typename TSDataIn>
class LevenshteinMetric
{

  private:
  std::unique_ptr<LevenshteinMetric_X<TSDataIn>> lvsm;

  public:
  typedef typename TSDataIn::Sample Sample;
  typedef LevenshteinInitializer    Initializer;
  LevenshteinMetric(const TSDataIn& datain, size_t nthreads, const LevenshteinInitializer& li);
  void set_distance(const Sample& v_vertical,
                    const Sample& v_horizontal,
                    double        threshold,
                    double&       distance);
  // void set_distance(const Sample & a, const Sample & b, double & distance);
  size_t get_ncalcs();
  double get_rel_calccosts();
  /* destructor declared elsewhere as unique_ptr to undefined class, as described at
   * https://stackoverflow.com/questions/27336779/unique-ptr-and-forward-declaration*/
  ~LevenshteinMetric();
};

}  // namespace nszen

#endif
