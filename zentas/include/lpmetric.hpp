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
  
  
  



    
    virtual void correct_threshold(double & threshold){
      (void)threshold;
    }
    

    /* Comment on squaring the threshold : We square the threshold 
     * as we'll be working with the squared distance. This might look
     *  dodgey: what if the threshold squared oveflows? This is not
     *  serious as it overflows to std::infinity, which is larger than
     *  any double anyway. Besides, it would be silly to complain 
     * here, as we're assuming that the distance squared is a valid 
     * double (which is a necessary assumption to make the computation).
     *  */    
    virtual void correct_distance(double & a_distance){
      (void)a_distance;
    }

    virtual void update_distance_partial(double & distance, double & diff) = 0;    

    virtual void set_distance(const SparseVectorSample<TNumber> &, const SparseVectorSample<TNumber> &, double, double &) = 0;


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

    void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & a_distance){
    //void set_distance(const TNumber * const &, const TNumber * const &, double, double & ) {
      ++ncalcs;
      a_distance = 0;       
      double diff;
      size_t d;
      correct_threshold(threshold);
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        update_distance_partial(a_distance, diff);
        if (a_distance > threshold){
          break;
        }
      }
      calccosts += (d == dimension ? d : 1 + d);
      correct_distance(a_distance);
    }    
        
};

template <typename TNumber>
class L2Distance : public LpDistance<TNumber>{
      
  public:
    L2Distance(size_t dimension):LpDistance<TNumber>(dimension){} 

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
     
    virtual void correct_threshold(double & threshold) override final{
      threshold *= threshold;
    }
    virtual void correct_distance(double & distance) override final{
      distance = std::sqrt(distance);
    }
    virtual void update_distance_partial(double & a_distance, double & diff) override final{
      a_distance += diff*diff;
    }
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;  
};





template <typename TNumber>
class L1Distance : public LpDistance<TNumber>{
      
  public:
    L1Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    
    
    virtual void update_distance_partial(double & a_distance, double & diff) override final{
      a_distance += std::abs(diff);
    }    

    
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;
};



template <typename TNumber>
class L0Distance : public LpDistance<TNumber>{
      
  public:
    L0Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    
    
    virtual void update_distance_partial(double & a_distance, double & diff) override final{
      a_distance += (diff != 0);
    }    

    
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;
};

template <typename TNumber>
class L_oo_Distance : public LpDistance<TNumber>{
      
  public:
    L_oo_Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    
    
    virtual void update_distance_partial(double & a_distance, double & diff) override final{
      a_distance = std::max(a_distance, std::abs(diff));
    }
    
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;
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
