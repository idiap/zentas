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
    LpMetricInitializer(char p);
    LpMetricInitializer();
    void reset(std::string metric);
};




template <typename TNumber>
class LpDistance{
  
  public:
  
    LpDistance(size_t dimension):dimension(dimension) {}
  
     size_t get_ncalcs(){
      return ncalcs;
    }
    
     size_t get_calccosts(){
      return calccosts;
    }
  
    size_t dimension;
    size_t ncalcs = 0;
    size_t calccosts = 0;

    /* the inherited classes have very similar set_distance functions, but with differences scattered within them.
     * code could be compactified by having set_distance be non-virtual, with a few calls to virtual functions within.
     * but virtual calls are slower than non-virtual calls (I have measured this) and for now sticking to current approach  */
    virtual void set_distance(const TNumber * const &, const TNumber * const &, double, double & ) = 0;
    virtual void set_distance(const SparseVectorSample<TNumber> &, const SparseVectorSample<TNumber> &, double, double &) = 0;
    
        
};

template <typename TNumber>
class L2Distance : public LpDistance<TNumber>{
      
  public:
    L2Distance(size_t dimension):LpDistance<TNumber>(dimension){} 

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    virtual void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final;
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;  
};





template <typename TNumber>
class L1Distance : public LpDistance<TNumber>{
      
  public:
    L1Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    virtual void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final;
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;
};



template <typename TNumber>
class L0Distance : public LpDistance<TNumber>{
      
  public:
    L0Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    virtual void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final;
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final;
};

template <typename TNumber>
class L_oo_Distance : public LpDistance<TNumber>{
      
  public:
    L_oo_Distance(size_t dimension):LpDistance<TNumber>(dimension) {}

    using LpDistance<TNumber>::dimension;
    using LpDistance<TNumber>::ncalcs;
    using LpDistance<TNumber>::calccosts;
    virtual void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) override final;
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
    
    
    //template <typename T>
     //void set_distance(const T & a, const T & b, double & distance){
      //set_distance(a, b, std::numeric_limits<double>::max(), distance); 
    //}
    
     size_t get_ncalcs() const{
      return uptr_lpdistance->ncalcs;
    }
    
     double get_rel_calccosts() const{
      return static_cast<double> (uptr_lpdistance->calccosts)  / static_cast<double> (uptr_lpdistance->ncalcs);
    }
  
  private:
    const char p;

};




class LevenshteinInitializer{
  
  public:
    size_t dict_size;
    double c_indel;
    double c_switch;
    /*  will use either the above (if dict_size == 0) or the below (otherwise). */
    const double *  c_indel_arr;
    const double *  c_switch_arr;
    bool normalised;

    LevenshteinInitializer(const size_t dict_size, const double c_indel, const double c_switch, const double * const c_indel_arr, const double * const c_switch_arr, bool normalised);
    LevenshteinInitializer(const double c_indel, const double c_switch, bool normalised);
    LevenshteinInitializer(const size_t dict_size, const double * const c_indel_arr, const double * const c_switch_arr, bool normalised);
    LevenshteinInitializer();
    
};


/* using pimpl for Levenshtein. This is the working class, */
template<typename TSDataIn>
class LevenshteinMetric_X;


/* using pimpl for Levenshtein. This is the wrapping class, */
template<typename TSDataIn>
class LevenshteinMetric{

  private:
   std::unique_ptr<LevenshteinMetric_X<TSDataIn>> lvsm;
  
  public:  
    typedef typename TSDataIn::Sample Sample;
    typedef LevenshteinInitializer Initializer;  
    LevenshteinMetric(const TSDataIn & datain, size_t nthreads, const LevenshteinInitializer & li);
    void set_distance(const Sample & v_vertical, const Sample & v_horizontal, double threshold, double & distance);
    //void set_distance(const Sample & a, const Sample & b, double & distance);
    size_t get_ncalcs();
    double get_rel_calccosts();
    /* destructor declared elsewhere as unique_ptr to undefined class, as described at
     * https://stackoverflow.com/questions/27336779/unique-ptr-and-forward-declaration*/ 
    ~LevenshteinMetric();
    

};


} //namespace nszen

#endif
