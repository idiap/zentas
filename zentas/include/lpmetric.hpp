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




#include "tdatain.hpp" 
/* the above is included for this guy:
 * template <typename TNumber>
struct SparseVectorSample; */

namespace nszen{
  
class LpMetricInitializer{

  public:
    char p;
    LpMetricInitializer(char p);
    LpMetricInitializer();
    void reset(std::string metric);
};



    /* Experiment with blas : 
     * with d = 1000, the speed-up was only 10s -> 7s. 
     * Not interesting enough to warrant the additional compilation hassle. 
     * it looked like this:
    >> ++ncalcs;
    >> wblas::copy(dimension, a, 1, worker, 1);
    >> wblas::axpy(dimension, static_cast<TNumber>(-1.), b, 1, worker, 1);
    >> distance = wblas::dot(dimension, worker, 1, worker, 1);
    >> calccosts += dimension; */    



template <typename TNumber>
class LpDistance{

  public:

    size_t dimension;
    size_t ncalcs = 0;
    size_t calccosts = 0;

    LpDistance(size_t dimension):dimension(dimension) {}
  
     size_t get_ncalcs(){
      return ncalcs;
    }
    
     size_t get_calccosts(){
      return calccosts;
    }
    
    virtual void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & a_distance) = 0;
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) = 0;
    
};


template <typename TNumber, class CorrectThreshold, class CorrectDistance, class UpdateDistance>
class TLpDistance : public LpDistance<TNumber>{

  public:
  
    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;

    TLpDistance(size_t dimension): LpDistance<TNumber>(dimension) {}

  private:
    CorrectThreshold correct_threshold;
    CorrectDistance correct_distance;
    UpdateDistance update_distance;
  
    
    void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & a_distance){
      ++ncalcs;
      a_distance = 0;       
      double diff;
      size_t d;
      correct_threshold(threshold);
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        update_distance(a_distance, diff);
        if (a_distance > threshold){
          break;
        }
      }
      calccosts += (d == dimension ? d : 1 + d);
      correct_distance(a_distance);
    }
    

    void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) {
              
      ++ncalcs;
      distance = 0;
      size_t a_pos = 0;
      size_t b_pos = 0;
      /* if both a and b have unchecked indices remaining */
      bool indices_remain = (a.size > 0 && b.size > 0);
      double diff;
    
      correct_threshold(threshold);
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
        
        update_distance(distance, diff);
        calccosts += 1;
      }
    
    
    
      while (a_pos != a.size && distance < threshold){
        diff = a.values[b_pos];
        update_distance(distance, diff);
        calccosts += 1;
        ++a_pos;
      }
      
      while (b_pos != b.size and distance < threshold){
        diff = b.values[b_pos];
        update_distance(distance, diff);
        calccosts += 1;
        ++b_pos;        
      }            
      
      correct_distance(threshold);
      
    }
};

class NullOpDouble{
  public:
    inline void operator() (double & x){ (void)x;}
};

class SquareDouble{
  public:
    inline void operator() (double & x){ x *= x;}
  
};

class RootDouble{
  public:
    inline void operator() (double & x){ x = std::sqrt(x);}  
};


class NonZero{
  public:
    inline void operator() (double & dist, const double & diff){ dist += (diff != 0);}  
};

class AbsDiff{
  public:
    inline void operator() (double & dist, const double & diff){ dist += std::abs(diff);}  
};

class SquareDiff{
  public:
    inline void operator() (double & dist, const double & diff){ dist += diff*diff;}  
};

class MaxDiff{
  public:
    inline void operator() (double & dist, const double & diff){ dist = std::max(dist, std::abs(diff));}  
};





template <typename TNumber>
using L0Distance = TLpDistance<TNumber, NullOpDouble, NullOpDouble, NonZero>;

template <typename TNumber>
using L1Distance = TLpDistance<TNumber, NullOpDouble, NullOpDouble, AbsDiff>;

template <typename TNumber>
using L2Distance = TLpDistance<TNumber, SquareDouble, RootDouble, SquareDiff>;

template <typename TNumber>
using L_oo_Distance = TLpDistance<TNumber, NullOpDouble, NullOpDouble, MaxDiff>;





  
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
    LpMetric(const TVDataIn & datain, size_t nthreads, const LpMetricInitializer & l2mi):p(l2mi.p) 
      {
      
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
    void set_distance(const T & a, const T & b, double threshold, double & distance){
      /* This virtual function call costs 
       * about 2% when dimension = 2. 2% 
       * slowdown is negligible, I'm going
       * for this clean code.  */
      uptr_lpdistance->set_distance(a, b, threshold, distance);
    }
        
     size_t get_ncalcs() const{
      return uptr_lpdistance->ncalcs;
    }
    
     double get_rel_calccosts() const{
      return static_cast<double> (uptr_lpdistance->calccosts)  / static_cast<double> (uptr_lpdistance->ncalcs);
    }
  
  private:
    const char p;

};

} //namespace nszen

#endif

