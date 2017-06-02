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
#ifndef ZENTAS_TMETRIC_H
#define ZENTAS_TMETRIC_H

/* 
 * 
 * define classes satisfying TMetric template

template <typename TDataIn> 
class TMetric{
  public: 
    TMetric(const TVDataIn & datain);
    inline void set_distance(const TVDataIn::Sample &, const TVDataIn::Sample &, double &) ; 
    get_ncalcs(); 
    get_rel_calccosts();
};

*/

#include <algorithm>
#include <mutex>
#include <atomic>
#include <vector>
#include <iostream>

#include "tdatain.hpp" 
/* the above is included for this guy:
 * template <typename TNumber>
struct SparseVectorSample; */

namespace nszen{
  


    

class LpMetricInitializer{

  
  public:
  
    char p;
    
    LpMetricInitializer(char p): p(p){
      
    }
    
    LpMetricInitializer(): p('!') {}
    
    void reset(std::string metric){
      if (metric.compare("li") == 0){
        p = 'i';
      }
        
      else if (metric.compare("l2") == 0){
        p = '2';
      }
      
      else if (metric.compare("l1") == 0){
        p = '1';
      }
    
      else if (metric.compare("l0") == 0){
        p = '0';
      }
      
      else{
        std::stringstream ss;
        ss << "Currently, only li (inf) & l0 & l1 & l2 metrics are implemented for vector data, not `" << metric << "'."; 
        throw zentas::zentas_error(ss.str());
      }
    }
};




template <typename TNumber>
class LpDistance{
  
  public:
  
    LpDistance(size_t dimension):dimension(dimension) {}
  
    inline size_t get_ncalcs(){
      return ncalcs;
    }
    
    inline size_t get_calccosts(){
      return calccosts;
    }
  
    size_t dimension;
    size_t ncalcs = 0;
    size_t calccosts = 0;

    /* the inherited classes have very similar set_distance functions, but with differences scattered within them.
     * code could be compactified by having set_distance be non-virtual, with a few calls to virtual functions within.
     * but virtual calls are slower than non-virtual calls (I have measured this) and for now sticking to current approach  */

    virtual inline void set_distance(const TNumber * const &, const TNumber * const &, double, double & ) {}
    
    virtual inline void set_distance(const SparseVectorSample<TNumber> &, const SparseVectorSample<TNumber> &, double, double &) {}
    
        
};

template <typename TNumber>
class L2Distance : public LpDistance<TNumber>{
      
  public:
    
    L2Distance(size_t dimension):LpDistance<TNumber>(dimension){} 

    
    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    
    virtual inline void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final{
      
      /* Experiment with blas : 
       * with d = 1000, the speed-up was only 10s -> 7s. 
       * Not interesting enough to warrant the additional compilation hassle. 
       * it looked like this:
      >> ++ncalcs;
      >> wblas::copy(dimension, a, 1, worker, 1);
      >> wblas::axpy(dimension, static_cast<TNumber>(-1.), b, 1, worker, 1);
      >> distance = wblas::dot(dimension, worker, 1, worker, 1);
      >> calccosts += dimension;
      */
              
      ++ncalcs;
      distance = 0;       
      double diff;
      size_t d;

      /* Comment on squaring the threshold : We square the threshold 
       * as we'll be working with the squared distance. This might look
       *  dodgey: what if the threshold squared oveflows? This is not
       *  serious as it overflows to std::infinity, which is larger than
       *  any double anyway. Besides, it would be silly to complain 
       * here, as we're assuming that the distance squared is a valid 
       * double (which is a necessary assumption to make the computation).
       *  */
      threshold *= threshold;
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        distance += diff*diff;
        if (distance > threshold){
          break;
        }
      }
      
      calccosts += (d == dimension ? d : 1 + d);
      distance = std::sqrt(distance);

    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;      
      double diff;
      size_t a_pos = 0;
      size_t b_pos = 0;
      bool indices_remain = (a.size > 0 && b.size > 0);

      threshold *= threshold;
      
      while (indices_remain && distance <= threshold){
        if (a.indices[a_pos] == b.indices[b_pos]){
          diff = a.values[a_pos] - b.values[b_pos];
          ++a_pos;
          ++b_pos;
          indices_remain = (a_pos < a.size && b_pos < b.size);
        }
        
        else if (a.indices[a_pos] < b.indices[b_pos]){
          diff = a.values[a_pos];
          ++a_pos;
          indices_remain = a_pos < a.size;
        }
        
        else{
          diff = b.values[b_pos];
          ++b_pos;
          indices_remain = b_pos < b.size;
        }
        
        distance += diff*diff;
        calccosts += 1;
      }
      
      while (a_pos != a.size and distance < threshold){
        diff = a.values[a_pos];
        distance += diff*diff;
        calccosts += 1;
        ++a_pos;
      }
      
      while (b_pos != b.size and distance < threshold){
        diff = b.values[b_pos];
        distance += diff*diff;
        calccosts += 1;
        ++b_pos;        
      }
      
      distance = std::sqrt(distance);
    }
  
};





template <typename TNumber>
class L1Distance : public LpDistance<TNumber>{
      
  public:
    L1Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;

    virtual inline void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final{
      
      ++ncalcs;
      distance = 0;       
      double diff;
      size_t d;
      
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        distance += std::abs(diff);
        if (distance > threshold){
          break;
        }
      }
      calccosts += (d == dimension ? d : 1 + d);
    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;      
      double diff;
      size_t a_pos = 0;
      size_t b_pos = 0;
      
      /* if both a and b have unchecked indices remaining */
      bool indices_remain = (a.size > 0 && b.size > 0); 

      while (indices_remain && distance <= threshold){
        if (a.indices[a_pos] == b.indices[b_pos]){
          diff = a.values[a_pos] - b.values[b_pos];
          ++a_pos;
          ++b_pos;
          indices_remain = (a_pos < a.size && b_pos < b.size);
        }
        
        else if (a.indices[a_pos] < b.indices[b_pos]){
          diff = a.values[a_pos];
          ++a_pos;
          indices_remain = a_pos < a.size;
        }
        
        else{
          diff = b.values[b_pos];
          ++b_pos;
          indices_remain = b_pos < b.size;
        }
        
        distance += std::abs(diff);
        calccosts += 1;
      }
      
      
      /* either or a or b or both have exhausted their indices, noe make both exhausted */
      while (a_pos != a.size && distance <= threshold){
        distance += std::abs(a.values[a_pos]);
        calccosts += 1;
        ++a_pos;
      }
      
      while (b_pos != b.size && distance <= threshold){
        distance += std::abs(b.values[b_pos]);
        calccosts += 1;
        ++b_pos;        
      }      
    }
};



template <typename TNumber>
class L0Distance : public LpDistance<TNumber>{
      
  public:
    L0Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;

    virtual inline void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final{
      ++ncalcs;
      distance = 0;       
      size_t d;
      for (d = 0; d < dimension; ++d){
        distance += (*(a + d) != *(b + d));
        if (distance > threshold){
          break;
        }
      }
      calccosts += (d == dimension ? d : 1 + d);
    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;
      size_t a_pos = 0;
      size_t b_pos = 0;
      bool indices_remain = (a.size > 0 && b.size > 0);
      
      //threshold += 0; //why :| :| :| ?
      while (indices_remain && distance <= threshold){
        if (a.indices[a_pos] == b.indices[b_pos]){
          distance += (a.values[a_pos] != b.values[b_pos]);
          ++a_pos;
          ++b_pos;
          indices_remain = (a_pos < a.size && b_pos < b.size);
        }
        
        else if (a.indices[a_pos] < b.indices[b_pos]){
          distance += 1;
          ++a_pos;
          indices_remain = a_pos < a.size;
        }
        
        else{
          distance += 1;
          ++b_pos;
          indices_remain = b_pos < b.size;
        }
        calccosts += 1;
      }
      
      if (distance <= threshold){
        distance += a.size - a_pos;
        distance += b.size - b_pos;
        calccosts += 1;
      }
      
    
    }
};

template <typename TNumber>
class L_oo_Distance : public LpDistance<TNumber>{
      
  public:
    L_oo_Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;

    virtual inline void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final{
      ++ncalcs;
      distance = 0;       
      size_t d;
      double diff;
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        distance = std::max(distance, std::abs(diff));
        if (distance > threshold){
          break;
        }
      }
      calccosts += (d == dimension ? d : 1 + d);
    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;
      size_t a_pos = 0;
      size_t b_pos = 0;
      bool indices_remain = (a.size > 0 && b.size > 0);
      double diff;
      while (indices_remain && distance <= threshold){
        if (a.indices[a_pos] == b.indices[b_pos]){
          diff = a.values[a_pos] - b.values[b_pos];
          ++a_pos;
          ++b_pos;
          indices_remain = (a_pos < a.size && b_pos < b.size);
        }
        
        else if (a.indices[a_pos] < b.indices[b_pos]){
          diff = a.values[a_pos];
          ++a_pos;
          indices_remain = a_pos < a.size;
        }
        
        else{
          diff = a.values[b_pos];
          ++b_pos;
          indices_remain = b_pos < b.size;
        }
        distance = std::max(distance, std::abs(diff));
        calccosts += 1;
      }


      
      while (a_pos != a.size && distance < threshold){
        diff = a.values[b_pos];
        distance = std::max(distance, std::abs(diff));
        calccosts += 1;
        ++a_pos;
      }
      
      while (b_pos != b.size and distance < threshold){
        diff = b.values[b_pos];
        distance = std::max(distance, std::abs(diff));
        calccosts += 1;
        ++b_pos;        
      }            
      
    }
};


  
/* where TVDataIn has public size_t dimension and typename TVDataIn::Sample */
template<typename TVDataIn>
class LpMetric{
  
  public:
    typedef typename TVDataIn::Sample Sample;
    typedef typename TVDataIn::AtomicType AtomicType;
    typedef LpMetricInitializer Initializer;    
  
  private:
    std::unique_ptr<LpDistance<AtomicType>>  uptr_lpdistance;

  public:
    LpMetric(const TVDataIn & datain, size_t nthreads, const LpMetricInitializer & l2mi):
    p(l2mi.p)  {
      
      //quelch warning.
      nthreads = 0;
      nthreads += 1;
      
      if (p == '2'){
        uptr_lpdistance.reset(new L2Distance<AtomicType>(datain.dimension));
      }
      
      else if (p == '1'){
        uptr_lpdistance.reset(new L1Distance<AtomicType>(datain.dimension));        
      }

      else if (p == '0'){
        uptr_lpdistance.reset(new L0Distance<AtomicType>(datain.dimension));        
      }

      else if (p == 'i'){
        uptr_lpdistance.reset(new L_oo_Distance<AtomicType>(datain.dimension));        
      }      
      
      else{
        std::string err = "No implementation for p = ";
        err = err + p;
        throw zentas::zentas_error(err);
      }
    }
    
    template <typename T>
    inline void set_distance(const T & a, const T & b, double threshold, double & distance){
      /* This virtual function call costs 
       * about 2% when dimension = 2. 2% 
       * slowdown is negligible, I'm going
       * for the clean code.  */
      uptr_lpdistance->set_distance(a, b, threshold, distance);
    }
    
    
    template <typename T>
    inline void set_distance(const T & a, const T & b, double & distance){
      set_distance(a, b, std::numeric_limits<double>::max(), distance); 
    }
    
    inline size_t get_ncalcs() const{
      return uptr_lpdistance->ncalcs;
    }
    
    inline double get_rel_calccosts() const{
      return static_cast<double> (uptr_lpdistance->calccosts)  / static_cast<double> (uptr_lpdistance->ncalcs);
    }
  
  private:
    const char p;

};







class MultiIndel{
  private:
    const double * values;
    
  public:
    MultiIndel(const double * const values):values(values) {}
    MultiIndel() {}
    
    inline double operator() (size_t i) const {
      return values[i];
    }
};

/* dict_size x dict_size switch costs for Levenshtein */
class MultiSwitch{
  private:
    const double * values;
    const size_t dict_size;
     
  public:
    MultiSwitch(const double * const values, size_t dict_size):values(values), dict_size(dict_size) {}
    
    inline double operator() (size_t i, size_t j) const {
      return values[i*dict_size + j];
    }
};

class ConstIndel{
  private:
    double value;
  
  public:
    ConstIndel(double value):value(value) {}
    ConstIndel() {}
    
    inline double operator() (size_t) const {
      return value;
    }
};

class ConstSwitch{
  private:
    double value;
    
    
  public:
    ConstSwitch(double value):value(value) {
    }
    
    inline double operator() (int i, int j) const {
      return value*(i!=j);
    }
};





/* Random comment. The use of min_c_indel might not be the optimal strategy in Levenshtein */
template <class TSample, class TIndelCost, class TSwitchCost>
inline void set_levenshtein_distance(const TSample & v_vertical, const TSample & v_horizontal, double threshold,
                                     size_t dict_size, const TIndelCost f_c_indel, double min_c_indel, 
                                     double max_c_indel,  const TSwitchCost f_c_switch, double * A_prev, 
                                     double *  A_acti, int nrows, int ncols, size_t & n_cells_visited_local, 
                                     double & distance){

  (void)dict_size;
  
  int length_difference = ncols - nrows; 
  if (length_difference < 0){
    throw zentas::zentas_error("ncols must not be smaller than nrows in set_levenshtein_distnance, ie shorter sequence first");
  }
  
  threshold = std::min(threshold, max_c_indel*(nrows + ncols));
  if (min_c_indel*length_difference >= threshold){
    distance = threshold;
  }
    
  else{
    int j_start; 
    int j_end;    
    int half_width = std::min(ncols - 1, static_cast<int>(std::floor(threshold / min_c_indel)));
    int width = 3 + 2*half_width;
    std::fill(A_prev, A_prev + width, threshold);
    std::fill(A_acti, A_acti + width, threshold);

    A_acti[1 + half_width] = 0;
    for (int j = 2 + half_width; j < width; ++j){
      A_acti[j] = A_acti[j-1] + f_c_indel(v_horizontal.values[j - 2 - half_width]);
    }
    int row = 0; 
    int column = 0;    
    double min_distance = 0;
    while (row < nrows && min_distance < threshold){
      
      std::swap(A_acti, A_prev);
      ++row;
      min_distance = threshold;
      j_start = std::max(1, 1 + half_width - row);
      j_end = std::min(ncols - row + 2 + half_width, width - 1);
      for (int j = j_start; j < j_end; ++j){
        column = row + j - (1 + half_width);

        A_acti[j] = std::min(
          std::min(
          f_c_indel(v_horizontal.values[column - 1]) + A_acti[j - 1], 
          f_c_indel(v_vertical.values[row - 1]) + A_prev[j + 1]),
          A_prev[j] + 
          f_c_switch(v_vertical.values[row - 1], v_horizontal.values[column - 1])
          );

        min_distance = std::min(min_distance, A_acti[j]);
      }
      n_cells_visited_local += j_end - j_start;
    }
    
    
    if (std::abs(nrows - ncols) <= half_width + 1){
      distance = std::min(threshold, A_acti[1 + half_width + ncols - nrows]);
    }
    else{
      distance = threshold;
    }
  }
}










class LevenshteinInitializer{
  
  public:
    size_t dict_size;
    double c_indel;
    double c_switch;
    /*  will use either the above (if dict_size == 0) or the below (otherwise). */
    const double *  c_indel_arr;
    const double *  c_switch_arr;
    bool normalised;

    LevenshteinInitializer(const size_t dict_size, const double c_indel, const double c_switch, const double * const c_indel_arr, const double * const c_switch_arr, bool normalised): dict_size(dict_size), c_indel(c_indel), c_switch(c_switch), c_indel_arr(c_indel_arr), c_switch_arr(c_switch_arr), normalised(normalised) {
      
      
      if (c_switch_arr != nullptr){
        /* confirm symmetry */
        for (unsigned i = 0; i < dict_size; ++i){
          for (unsigned j = 0; j < dict_size; ++j){  
            if (c_switch_arr[i*dict_size + j] != c_switch_arr[j*dict_size + i]){
              std::stringstream ss;
              ss << "cost_switch is not symmetric, it should be ";
              ss << "(" << c_switch_arr[i*dict_size + j] << " != " << c_switch_arr[j*dict_size + i] << ")";
              throw zentas::zentas_error(ss.str());
            }
            if (c_switch_arr[i*dict_size + j] < 0){
              throw zentas::zentas_error("cost_switch should contain no negative values, it does ");
            }
          }
        }
        
        /* confirm triangle inequality holds. may be a quicker way to do this. */
        for (unsigned i = 0; i < dict_size; ++i){
          for (unsigned j = 0; j < dict_size; ++j){
            if (i != j){
              double d_ij = c_switch_arr[dict_size*i + j];
              for (unsigned k = 0; k < dict_size; ++k){
                double d_ik = c_switch_arr[dict_size*i + k];
                double d_jk = c_switch_arr[dict_size*j + k];
                if (d_ij > d_ik + d_jk){
                  std::stringstream ss;
                  ss << "the cost_switch matrix does not obey the triangle inequality, \n";
                  ss << "d[" << i << "][" << j << "]" <<  " > "  << "d[" << i << "][" << k <<  "]" << " + " << "d[" << j << "][" << k << "], \n";
                  ss << d_ij << " > " << d_ik  << " + " << d_jk << ".\n";
                  throw zentas::zentas_error(ss.str());
                }
              }
            }
          }
        }
        
        
        if (normalised == true){
          for (unsigned i = 0; i < dict_size; ++i){
            if (c_indel_arr[i] != c_indel_arr[0]){
              std::stringstream ss;
              ss << "The normalised Levenshtein metric is only a true metric when the indel costs are the same, but cost_indel[" << i << "] is not the same as cost_indel[0]";
              ss << " (" << c_indel_arr[i] << ") vs (" << c_indel_arr[0] << ").";
              throw zentas::zentas_error(ss.str());
            }
          }
        }
      }
    }
    
    LevenshteinInitializer(const double c_indel, const double c_switch, bool normalised): LevenshteinInitializer(0, c_indel, c_switch, nullptr, nullptr, normalised) {}
    
    LevenshteinInitializer(const size_t dict_size, const double * const c_indel_arr, const double * const c_switch_arr, bool normalised): LevenshteinInitializer(dict_size, 0., 0., c_indel_arr, c_switch_arr, normalised){}
    
    LevenshteinInitializer():LevenshteinInitializer(0,0.,0.,nullptr, nullptr, false) {}
    
};

template<typename TSDataIn>
/* See levenshtein_prototyping for equivalent python version */
class LevenshteinMetric{
    
  private:
    std::vector<size_t> v_ncalcs;   
    const size_t dict_size; 
    bool normalised;
    
    ConstIndel cv_c_indel;
    ConstSwitch cv_c_switch;
    MultiIndel av_c_indel;
    MultiSwitch av_c_switch;
    
    std::vector<size_t> v_n_cells_visited;
    std::vector<size_t> v_n_cells_visitable;
    int max_size;
    int memory_size;
    size_t nthreads;
    std::vector<std::mutex> v_mutex0;
    std::mutex counter_mutex;    
    std::vector<std::unique_ptr<double []> > v_up_A_acti;
    std::vector<double * > v_A_acti;
    std::vector< std::unique_ptr<double []> > v_up_A_prev;
    std::vector< double * > v_A_prev;
    
    double min_c_indel;
    double max_c_indel;

    bool test_where_possible;
    
  protected:  
    std::mt19937_64 gen;
    std::uniform_int_distribution<unsigned> dis;
  
  public:
    typedef typename TSDataIn::Sample Sample;
    typedef LevenshteinInitializer Initializer;
    
    LevenshteinMetric(const TSDataIn & datain, size_t nthreads, const LevenshteinInitializer & li): 
    v_ncalcs(nthreads, 0),

    dict_size(li.dict_size),
    normalised(li.normalised),

    cv_c_indel(li.c_indel), 
    cv_c_switch(li.c_switch), 
    av_c_indel(li.c_indel_arr), 
    av_c_switch(li.c_switch_arr, dict_size),    
    
     
    v_n_cells_visited(nthreads, 0), 
    v_n_cells_visitable(nthreads, 0), 
    max_size(datain.get_max_size()), 
    
    /* below : making 10*size + 10 makes no difference to performance (speed) */
    memory_size(4*datain.get_max_size() + 10), 
    nthreads(nthreads),
    v_mutex0(nthreads)
    
    {
      for (size_t ti = 0; ti < nthreads; ++ti){
        v_up_A_acti.emplace_back(new double [memory_size]);
        v_A_acti.push_back(v_up_A_acti[ti].get());

        v_up_A_prev.emplace_back(new double [memory_size]);
        v_A_prev.push_back(v_up_A_prev[ti].get());
      }
      
      if (dict_size > 0){
        min_c_indel = std::numeric_limits<double>::max();
        max_c_indel = 0.;
        for (size_t w = 0; w < dict_size; ++w){
          min_c_indel = std::min(min_c_indel, av_c_indel(w));
          max_c_indel = std::max(max_c_indel, av_c_indel(w));
        }
      }
      
      else{
        min_c_indel = cv_c_indel(0);
        max_c_indel = cv_c_indel(0);
      }
      
      
      test_where_possible = false;
      if (test_where_possible == true){
        std::cerr << "\n\nLEVENSHTEIN WITH (LIMITED) TESTS ENABLED : WILL BE SLOWER" << std::endl;        
      }
    }
    
    inline void set_distance_simple_test(const Sample & v_vertical, const Sample & v_horizontal, double threshold, double & distance){
      (void)threshold;      
      distance = std::abs(int(v_vertical.size) - int(v_horizontal.size));
    }
    

    inline void set_distance(const Sample & v_vertical, const Sample & v_horizontal, double threshold, double & distance){      
      
      /* make sure the shorter vector comes first */
      if (v_vertical.size < v_horizontal.size) {
        set_distance_tiffany(v_vertical, v_horizontal, threshold, distance);
      }
      else{
        set_distance_tiffany(v_horizontal, v_vertical, threshold, distance);
      }
      
      if (test_where_possible && v_vertical.size == v_horizontal.size){
        double d1;
        double d2;
        set_distance_tiffany(v_vertical, v_horizontal, threshold, d1);
        set_distance_tiffany(v_horizontal, v_vertical, threshold, d2);
        
        if (d1 != d2){
          std::stringstream ss;
          ss << "Haha! Levenshtein distances not the same when order reversed  \n";
          ss << v_vertical.str() << " ->  " << v_horizontal.str() <<  " : " << d1 << "\n";
          ss << v_horizontal.str() << " ->  " << v_vertical.str() <<  " : " << d2 << "\n";
          throw zentas::zentas_error(ss.str());
        }        
      }
    }

    inline void set_distance_tiffany(const Sample & v_vertical, const Sample & v_horizontal, double threshold, double & distance) {
      
 
      /* numerical issues */
      threshold *= 1.0000230507000110130001701900023;

      /* n_d = 2d / ( alpha (L1 + L2) + d )  where alpha = max indel cost
       * => d = alpha n_d (L1 + L2) / (2 - n_d)
       * */
      if (normalised == true){
        // threshold in normalised space
        threshold = std::min(1., threshold);         
        //threshold in non-normalised space, where the distance calculation is to take place.
        threshold = max_c_indel * threshold * static_cast<double>(v_vertical.size + v_horizontal.size) / (2. - threshold); 
      }
           
     
      int nrows = static_cast<int>(v_vertical.size);
      int ncols = static_cast<int>(v_horizontal.size);
      size_t n_cells_visited_local = 0;

      // get a mutex and keep it
      size_t mutex_i = dis(gen)%nthreads;
      
      while (v_mutex0[mutex_i].try_lock() == false){
        ++mutex_i;
        mutex_i %= nthreads;
      }
        
      std::lock_guard<std::mutex> lock (v_mutex0[mutex_i], std::adopt_lock);
      double * A_acti = v_A_acti[mutex_i];
      double * A_prev = v_A_prev[mutex_i];

      /* do some fast tracking tests */
      if (nrows == 0 || ncols == 0){
        throw zentas::zentas_error("empty string, I need to confirm that this is not a special case. remind me to do this !");
      }
      
      else{
        /* constant indel and switch cost */
        if (dict_size == 0){ 
         set_levenshtein_distance(v_vertical, v_horizontal, threshold, dict_size, cv_c_indel, min_c_indel, max_c_indel, cv_c_switch, A_prev, A_acti, nrows, ncols,  n_cells_visited_local, distance);
        }
        
        /* matrices of costs */
        else{
          set_levenshtein_distance(v_vertical, v_horizontal, threshold, dict_size, av_c_indel, min_c_indel, max_c_indel, av_c_switch, A_prev, A_acti, nrows, ncols, n_cells_visited_local, distance);
        }
      }


      if (normalised == true){
        /* return to the normalised space */
        distance = 2.*distance / (max_c_indel * static_cast<double> ( v_horizontal.size + v_vertical.size ) + distance ); 
      }
      
      ++v_ncalcs[mutex_i];
      v_n_cells_visitable[mutex_i] += nrows*ncols;
      v_n_cells_visited[mutex_i] += n_cells_visited_local;
      
      distance = static_cast<double>(static_cast<float> (distance));
    }
          
    inline void set_distance(const Sample & a, const Sample & b, double & distance) {
      set_distance(a, b, std::numeric_limits<double>::max(), distance);
    }
    
    
    
    
    inline size_t get_ncalcs() const{
      size_t ncalcs = 0;
      for (auto & x : v_ncalcs){
        ncalcs += x;
      }
      return ncalcs;
    }
    
    inline double get_rel_calccosts() const{
      
      /* How well have we done as compared to the O(row*column) algorithm ?*/
      size_t n_cells_visited = 0;
      for (auto & x : v_n_cells_visited){
        n_cells_visited += x;
      }
 
      size_t n_cells_visitable = 0;
      for (auto & x : v_n_cells_visitable){
        n_cells_visitable += x;
      }
      
      return static_cast<double> (n_cells_visited) / static_cast<double> (n_cells_visitable);
    }
    
};

} //namespace nszen

#endif
