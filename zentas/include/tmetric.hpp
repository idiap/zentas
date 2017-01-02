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

//#include "blastemplates.h"

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
        throw std::runtime_error("Currently, only li (inf) & l0 & l1 & l2 metrics are implemented for vector data");
      }
    }
};




template <typename TNumber>
class LpDistance{
  
  public:
  
    LpDistance(size_t dimension):dimension(dimension) {}
    //LpDistance(){}
  
    inline size_t get_ncalcs(){
      return ncalcs;
    }
    
    inline size_t get_calccosts(){
      return calccosts;
    }
  
    size_t dimension;
    size_t ncalcs = 0;
    size_t calccosts = 0;

    //virtual inline void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) = 0;
    virtual inline void set_distance(const TNumber * const &, const TNumber * const &, double, double & ) {}
    
    //virtual inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) = 0;
    virtual inline void set_distance(const SparseVectorSample<TNumber> &, const SparseVectorSample<TNumber> &, double, double &) {}
    
};

template <typename TNumber>
class L2Distance : public LpDistance<TNumber>{
/* Comment on squaring the threshold : We square the threshold as we'll be working with the squared distance. This might look dodgey: what if the threshold squared oveflows? This is not so serious as it overflows to std::infi nity, which is larger than any double anyway. Besides, it would be hypocritical to complain here, as we're assuming that the distance squared is a valid double (which is a necessary assumption to make the computation). */
      
  public:
  
    //std::unique_ptr<TNumber []> uptr_worker;
    //TNumber * const worker;
    
    L2Distance(size_t dimension):LpDistance<TNumber>(dimension){} //, uptr_worker (new TNumber [dimension]), worker(uptr_worker.get()) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;

    virtual inline void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final{
      
      
      /* Experiment with blas : with d = 1000, the speed-up was only 10s -> 7s. Not interesting enough to warrant the additional compilation hassle. 
      ++ncalcs;
      wblas::copy(dimension, a, 1, worker, 1);
      wblas::axpy(dimension, static_cast<TNumber>(-1.), b, 1, worker, 1);
      distance = wblas::dot(dimension, worker, 1, worker, 1);
      calccosts += dimension;
      */



      ++ncalcs;
      distance = 0;       
      threshold *= threshold;
      double diff;
      size_t d;      
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        distance += diff*diff;
        //checking every time is too much, TODO think about changing this. TODO blasify this.
        if (distance > threshold){
          break;
        }
      }
      calccosts += d;
      distance = std::sqrt(distance);



    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;      
      threshold *= threshold;
      double diff;
      size_t a_pos = 0;
      size_t b_pos = 0;
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
        //checking every time is too much, TODO think about changing this
        if (distance > threshold){
          break;
        }
      }
      calccosts += d;
    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;      
      double diff;
      size_t a_pos = 0;
      size_t b_pos = 0;
      bool indices_remain = (a.size > 0 && b.size > 0); // both a and b have uncheck indices.

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
      
      
      //either or a or b or both have exhausted their indices, make both exhausted.
      
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
        //checking every time is too much, TODO think about changing this
        if (distance > threshold){
          break;
        }
      }
      calccosts += d;
    }

    inline void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
      ++ncalcs;
      distance = 0;
      size_t a_pos = 0;
      size_t b_pos = 0;
      bool indices_remain = (a.size > 0 && b.size > 0); // neither have exhausted their indices
      
      threshold += 100;
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
        //either or a or b or both have exhausted their indices      
        distance += a.size - a_pos;
        distance += b.size - b_pos;
        calccosts += 1;
        //calccosts += a.size - a_pos;
        //calccosts += b.size - b_pos;
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
        //checking every time is too much, TODO think about changing this
        if (distance > threshold){
          break;
        }
      }
      calccosts += d;
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
        throw std::runtime_error(err);
      }
    }
    //TODO : mutex lock for calccosts increment NO! same approach as Levenshtein : vector of calccosts, added when get_ncalcs called. 
    
    template <typename T>
    inline void set_distance(const T & a, const T & b, double threshold, double & distance){
      //This virtual function call costs about 2% when dimension = 2. 2% slowdown is negligible, I'm going for the clean code.  
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
//    const size_t dimension;
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

class MultiSwitch{
  private:
    const double * values;
    const size_t dict_size;
     
  public:
    MultiSwitch(const double * const values, size_t dict_size):values(values), dict_size(dict_size) {}
//    MultiSwitch(): values(nullptr), dict_size(0) {}
    
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




template <class TSample, class TIndelCost, class TSwitchCost>
inline void set_levenshtein_distance(const TSample & v_vertical, const TSample & v_horizontal, double threshold,

size_t dict_size, const TIndelCost f_c_indel, double min_c_indel, double max_c_indel,  const TSwitchCost f_c_switch, 
// where (above) c_indel has dict_size values and c_switch has dict_size * dict_size values

double * A_prev, double *  A_acti, int nrows, int ncols, size_t & n_cells_visited_local, double & distance){
  /* Comment to put somewhere propicious :  The use of min_c_indel is perhaps not the optimal strategy */

  //quelch warning
  dict_size += 0;
  
  int length_difference = std::abs(nrows - ncols);
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
      A_acti[j] = A_acti[j-1] + f_c_indel(v_horizontal.values[j - 2 - half_width]); // [v_horizontal[j - 1 - half_width]]*
    }
    int row = 0; 
    int column = 0;    
    double min_distance = 0;
    while (row < nrows && min_distance < threshold){
      
      //std::fill(scores.begin(), scores.end(), threshold);
      
      std::swap(A_acti, A_prev);
      ++row;
      min_distance = threshold;
      j_start = std::max(1, 1 + half_width - row);
      j_end = std::min(ncols - row + 2 + half_width, width - 1);
      for (int j = j_start; j < j_end; ++j){
        column = row + j - (1 + half_width);
        //A_acti[j] = std::min(
        //std::min(c_indel + A_acti[j - 1], c_indel + A_prev[j + 1]),
          //A_prev[j] + c_switch*(v_vertical.values[row - 1] != v_horizontal.values[column - 1])
        //);
        A_acti[j] = std::min(
          std::min(
          f_c_indel(v_horizontal.values[column - 1]) + A_acti[j - 1], 
          f_c_indel(v_vertical.values[row - 1]) + A_prev[j + 1]),
          A_prev[j] + 
          f_c_switch(v_vertical.values[row - 1], v_horizontal.values[column - 1])
          //f_c_switch(v_vertical.values[row - 1]*dict_size + v_horizontal.values[column - 1])
          //*(v_vertical.values[row - 1] != v_horizontal.values[column - 1]) //TODO : make faster and cleaner (as before).
          
          );
                
        min_distance = std::min(min_distance, A_acti[j]);
        
        //scores[column] = std::min(A_acti[j], threshold);
 
      }
      
      
      n_cells_visited_local += j_end - j_start;
    }
    
    
    if (std::abs(nrows - ncols) <= half_width + 1){
      distance = std::min(threshold, A_acti[1 + half_width + ncols - nrows]);
    }
    else{
      //throw std::runtime_error("length_difference <= half_width +1 should not be possible, there is a flaw in algorithm's logic");
      distance = threshold;
    }
  }
}










class LevenshteinInitializer{
  
  public:
    size_t dict_size;
    double c_indel;
    double c_switch;
    // will use either the above (if dict_size == 0) or the below (otherwise).
    const double *  c_indel_arr;
    const double *  c_switch_arr;
    bool normalised;

    LevenshteinInitializer(const size_t dict_size, const double c_indel, const double c_switch, const double * const c_indel_arr, const double * const c_switch_arr, bool normalised): dict_size(dict_size), c_indel(c_indel), c_switch(c_switch), c_indel_arr(c_indel_arr), c_switch_arr(c_switch_arr), normalised(normalised) {}
    
    LevenshteinInitializer(const double c_indel, const double c_switch, bool normalised): LevenshteinInitializer(0, c_indel, c_switch, nullptr, nullptr, normalised) {}
    
    LevenshteinInitializer(const size_t dict_size, const double * const c_indel_arr, const double * const c_switch_arr, bool normalised): LevenshteinInitializer(dict_size, 0., 0., c_indel_arr, c_switch_arr, normalised){}
    
    LevenshteinInitializer():LevenshteinInitializer(0,0.,0.,nullptr, nullptr, false) {}
  
//    LevenshteinInitializer& operator= ( const LevenshteinInitializer & ) = default;	
  
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
    //// will use either the above (if dict_size > 0) or the below (if dict_size == 0).
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
    memory_size(4*datain.get_max_size() + 10), //making 10*size + 10 makes no difference to performance (speed)
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
    }
    
    inline void set_distance_simple_test(const Sample & v_vertical, const Sample & v_horizontal, double threshold, double & distance){
      //quelch warning
      threshold += 0;
      
      distance = std::abs(int(v_vertical.size) - int(v_horizontal.size));
    }
    

    inline void set_distance(const Sample & v_vertical, const Sample & v_horizontal, double threshold, double & distance) {
      
 
      /* numerical issues */
      threshold *= 1.00001;

      /* n_d = 2d / ( alpha (L1 + L2) + d )  where alpha = max indel cost
       * => d = alpha n_d (L1 + L2) / (2 - n_d)
       * 
       * */
      if (normalised == true){
        
        // threshold in normalised space
        threshold = std::min(1., threshold); 
        
        //threshold in non-normalised space, where the distance calculation is to take place.
        threshold = max_c_indel * threshold * (v_vertical.size + v_horizontal.size) / (2. - threshold); 
      }
           
     
      int nrows = static_cast<int>(v_vertical.size);
      int ncols = static_cast<int>(v_horizontal.size);
      size_t n_cells_visited_local = 0;

      
      //get a mutex and keep
      size_t mutex_i = dis(gen)%nthreads;
      
      while (v_mutex0[mutex_i].try_lock() == false){
        ++mutex_i;
        mutex_i %= nthreads;
      }
        
      std::lock_guard<std::mutex> lock (v_mutex0[mutex_i], std::adopt_lock);
      double * A_acti = v_A_acti[mutex_i];
      double * A_prev = v_A_prev[mutex_i];


      /* do some fast tracking tests */
      double length_difference  = std::abs(nrows - ncols);
      if (nrows == 0 || ncols == 0){
        throw std::runtime_error("empty string, I need to confirm that this is not a special case. remind me to do this !");
        distance = 3.1415*length_difference;
      }
      
      else{
        if (dict_size == 0){ // constant indel and switch cost (switch is either 0 or c_switch).
         set_levenshtein_distance(v_vertical, v_horizontal, threshold, dict_size, cv_c_indel, min_c_indel, max_c_indel, cv_c_switch, A_prev, A_acti, nrows, ncols,  n_cells_visited_local, distance);

        }
        
        else{ // variable indel and switch costs. 
          set_levenshtein_distance(v_vertical, v_horizontal, threshold, dict_size, av_c_indel, min_c_indel, max_c_indel, av_c_switch, A_prev, A_acti, nrows, ncols, n_cells_visited_local, distance);
        }
      }


      if (normalised == true){
        //return to the normalised space
        distance = 2.*distance / (max_c_indel * ( v_horizontal.size + v_vertical.size ) + distance ); 
      }
      
      ++v_ncalcs[mutex_i];
      v_n_cells_visitable[mutex_i] += nrows*ncols;
      v_n_cells_visited[mutex_i] += n_cells_visited_local;
      
  
      
    }
          
    inline void set_distance(const Sample & a, const Sample & b, double & distance) {
      set_distance(a, b, std::numeric_limits<double>::max(), distance);
      //TODO : If there is no threshold, one might do better with a simpler algorithm. 
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
