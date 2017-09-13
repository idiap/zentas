// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_LPLPLPMETRIC_H
#define ZENTAS_LPLPLPMETRIC_H

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
#include "sparsevectorrfcenter.hpp"
/* the above is included for this guy:
 * template <typename TNumber>
struct SparseVectorSample; */

namespace nszen{
  
class LpMetricInitializer{

  public:
    char p;
    bool do_refinement;
    std::string rf_alg;
    size_t rf_max_rounds;
    double rf_max_time;


    LpMetricInitializer(char p, bool do_refinement, std::string rf_alg, size_t rf_max_rounds, double rf_max_time);
    LpMetricInitializer();
    void reset(std::string metric, bool do_refinement, std::string rf_alg, size_t rf_max_rounds, double rf_max_time);
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
class BaseLpDistance{

  public:

    size_t dimension;
    size_t ncalcs = 0;
    size_t calccosts = 0;

    BaseLpDistance(size_t dimension):dimension(dimension) {}
  
     size_t get_ncalcs(){
      return ncalcs;
    }
    
     size_t get_calccosts(){
      return calccosts;
    }
    
    virtual void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & a_distance) = 0;
    virtual void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) = 0;
    virtual void set_distance(const SparseVectorRfCenter<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) = 0;
    virtual void set_distance(const SparseVectorRfCenter<TNumber> & a, const SparseVectorRfCenter<TNumber> & b, double threshold, double & distance) = 0;
    

    virtual void set_center(const std::string & energy, 
    const std::function<const TNumber * const (size_t)> & f_sample, 
    size_t ndata, TNumber * const ptr_center){
      (void)energy; (void)f_sample; (void)ndata; (void)ptr_center;
      throw zentas::zentas_error("cannot set vector center : from BaseLpDistance base class dense");
    }
    
    virtual void set_center(const std::string & energy, 
    const std::function<const SparseVectorSample<TNumber> & (size_t)> & f_sample, 
    size_t ndata, SparseVectorRfCenter<TNumber> * const ptr_center){
      (void)energy; (void)f_sample; (void)ndata; (void)ptr_center;
      throw zentas::zentas_error("cannot set vector center : from BaseLpDistance base class sparse");
    }


    void update_sum(const std::string & energy, 
    const std::function<const TNumber * const (size_t)> & f_add, size_t n_add, 
    const std::function<const TNumber * const (size_t)> & f_subtract,  size_t n_subtract,
    TNumber * const ptr_sum){
      (void)energy; (void)f_add; (void)f_subtract; (void)n_add; (void)n_subtract; (void)ptr_sum;
      throw zentas::zentas_error("cannot set vector center : from BaseLpDistance base class dense");
    }
    
    void update_sum(const std::string & energy, 
    const std::function<const SparseVectorSample<TNumber> & (size_t)> & f_add, size_t n_add,
    const std::function<const SparseVectorSample<TNumber> & (size_t)> & f_subtract, size_t n_subtract,
    SparseVectorRfCenter<TNumber> * const ptr_sum){
      (void)energy; (void)f_add; (void)f_subtract; (void)n_add; (void)n_subtract; (void)ptr_sum;
      throw zentas::zentas_error("cannot set vector center : from BaseLpDistance base class sparse");
    }

    
};

template <typename TNumber>
class L2Minimiser{
  public:
  
  size_t dimension;
  L2Minimiser(size_t dimension_):dimension(dimension_) {};

    void set_center(const std::string & energy, 
    const std::function<const TNumber * const (size_t)> & f_sample, 
    size_t ndata, TNumber * const ptr_center){
      if (energy != "quadratic"){
        if (energy == "identity"){
          throw zentas::zentas_error(
          "request to find the center with L2 metric and identity energy not implemented."
          "I'm considering the algorithm in 'Geometric Median in Nearly Linear Time' of Cohen et al. "
          "Currently only quadratic energy is supported for L2 metric (so vanilla k-means). ");
        }
        throw zentas::zentas_error("cannot yet perform L2 metric minimisation unless the energy is quadratic");
      }
      
      for (size_t d = 0; d < dimension; ++d){
        ptr_center[d] = 0;
      }
      
      const TNumber * sample;
      for (size_t j = 0; j < ndata; ++j){
        sample = f_sample(j);
        for (size_t d = 0; d < dimension; ++d){
          ptr_center[d] += sample[d];
        }
      }
      for (size_t d = 0; d < dimension; ++d){
        ptr_center[d] /= static_cast<TNumber>(ndata);
      }
    }
    
    void set_center(const std::string & energy, 
    const std::function<const SparseVectorSample<TNumber> & (size_t)> & f_sample, 
    size_t ndata, SparseVectorRfCenter<TNumber> * const minimiser){
      (void)energy; (void)f_sample; (void)ndata; (void)minimiser;
      throw zentas::zentas_error("cannot set vector centers for sparse data yet : from L2Minimiser base class sparse. ");
    }  
};


template <typename TNumber>
class L1Minimiser{
  public:
  
  size_t dimension;
  L1Minimiser(size_t dimension_):dimension(dimension_) {};

    void set_center(const std::string & energy, 
    const std::function<const TNumber * const (size_t)> & f_sample, 
    size_t ndata, TNumber * const ptr_center){
      if (energy != "identity"){
        throw zentas::zentas_error("cannot yet perform L1 metric minimisation unless the energy is identity (ie loss is sum of l1 distances to nearest centers)");
      }
      
      std::vector<TNumber> data_by_dimension(dimension*ndata);
      const TNumber * sample;
      for (size_t j = 0; j < ndata; ++j){
        sample = f_sample(j);
        for (size_t d = 0; d < dimension; ++d){
          data_by_dimension[d*ndata + j] += sample[d];
        }
      }
      
      //get medians in each dimension.
      for (size_t d = 0; d < dimension; ++d){
        std::nth_element(data_by_dimension.begin() + d*ndata, data_by_dimension.begin() + d*ndata + ndata/2, data_by_dimension.begin() + (d+1)*ndata);
        ptr_center[d] = data_by_dimension[d*ndata + ndata/2];
      }
    }
    
    void set_center(const std::string & energy, 
    const std::function<const SparseVectorSample<TNumber> & (size_t)> & f_sample, 
    size_t ndata, SparseVectorRfCenter<TNumber> * const minimiser){
      (void)energy; (void)f_sample; (void)ndata; (void)minimiser;
      throw zentas::zentas_error("cannot set vector centers for sparse data yet : from L1Minimiser base class sparse. ");
    }  
};




class NoimplMinimiser{
  public:
  size_t dimension;
  NoimplMinimiser(size_t dimension_):dimension(dimension_) {};
 
  template <typename TPointerCenter, typename TFunction>
  void set_center(const std::string & energy, const TFunction & f_sample, size_t ndata, TPointerCenter ptr_center){
    (void)energy; (void)f_sample; (void)ndata; (void)ptr_center;
    throw zentas::zentas_error("unable to set vector space center. "
    "Recall that refinement consists of alternating between updating labels, and updating centers. "
    "The problem here is in updating the centers, zentas does not currently support all metric-energy combinations here. "
    "(Should I implement a generic convex solver for all convex energies?)");
  } 
};

template <typename TNumber, class CorrectThreshold, class CorrectDistance, class UpdateDistance, class FMinimiser>
class TLpDistance : public BaseLpDistance<TNumber> {

  public:
  
    using BaseLpDistance<TNumber>::dimension;
    using BaseLpDistance<TNumber>::ncalcs;
    using BaseLpDistance<TNumber>::calccosts;

  private:
    CorrectThreshold correct_threshold;
    CorrectDistance correct_distance;
    UpdateDistance update_distance;
    FMinimiser f_minimiser;

  public:
    TLpDistance(size_t dimension): BaseLpDistance<TNumber>(dimension), f_minimiser(dimension) {}

  
    
    void set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & a_distance) override final{
      ++ncalcs;
      a_distance = 0;       
      double diff;
      size_t d;
      correct_threshold(threshold);
      for (d = 0; d < dimension; ++d){
        diff = *(a + d) - *(b + d);
        update_distance(a_distance, diff);
        if (a_distance > threshold){
          a_distance = std::numeric_limits<double>::max();
          calccosts += (1 + d);
          return;
        }
      }
      calccosts += dimension;
      correct_distance(a_distance);
    }
    

    void set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) override final{

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
          diff = b.values[b_pos];
          ++b_pos;
          indices_remain = b_pos < b.size;
        }
        
        update_distance(distance, diff);
        calccosts += 1;
      }
    
    
    
      while (a_pos != a.size && distance < threshold){
        diff = a.values[a_pos];
        update_distance(distance, diff);
        calccosts += 1;
        ++a_pos;
      }
      
      while (b_pos != b.size && distance < threshold){
        diff = b.values[b_pos];
        update_distance(distance, diff);
        calccosts += 1;
        ++b_pos;        
      }
      
      correct_distance(distance);
    }
    
    void set_distance(const SparseVectorRfCenter<TNumber> & c_in, const SparseVectorSample<TNumber> & x, double threshold, double & distance) {

      ++ncalcs;
      distance = 0;
      correct_threshold(threshold);
      
      SparseVectorRfCenter<TNumber> c(c_in);
      
      for (size_t d = 0; d < x.size; ++d){
        
        if (c.count(x.indices[d]) == 1){
          update_distance(distance, c.at(x.indices[d]) - x.values[d]);
          c[x.indices[d]] = 0;
        }
        else{
          update_distance(distance, x.values[d]);
          calccosts += 1;
        }
      }
      
      if (distance < threshold){
        for(auto iter = c.cbegin(); iter != c.cend(); ++iter){
          update_distance(distance, iter->second);
          calccosts += 1;
        }
      }
      correct_distance(distance);
    }

    void set_distance(const SparseVectorRfCenter<TNumber> & c1, const SparseVectorRfCenter<TNumber> & c2, double threshold, double & distance) {
      
      ++ncalcs;
      distance = 0;
      correct_threshold(threshold);
      
      for(auto iter = c1.cbegin(); iter != c1.cend(); ++iter){
        calccosts += 1;
        if (c2.count(iter->first) == 0){
          update_distance(distance, iter->second);
        }
        else{
          update_distance(distance, iter->second - c2.at(iter->first));
        }
        if (distance > threshold){
          break;
        }
      }
      
      for(auto iter = c2.cbegin(); iter != c2.cend(); ++iter){
        if (c1.count(iter->first) == 0){
          update_distance(distance, iter->second);
          calccosts += 1;
        }
        if (distance > threshold){
          break;
        }
      }
      
      correct_distance(distance);
    }


    virtual void set_center(const std::string & energy, const std::function<const TNumber * const (size_t)> & f_sample, size_t ndata, TNumber * const ptr_center) override final{
      f_minimiser.set_center(energy, f_sample, ndata, ptr_center);
    }
    
    virtual void set_center(const std::string & energy, const std::function<const SparseVectorSample<TNumber> & (size_t)> & f_sample, size_t ndata, SparseVectorRfCenter<TNumber> * const ptr_center) override final{
      f_minimiser.set_center(energy, f_sample, ndata, ptr_center);
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
using L0Distance = TLpDistance<TNumber, NullOpDouble, NullOpDouble, NonZero, NoimplMinimiser>;

template <typename TNumber>
using L1Distance = TLpDistance<TNumber, NullOpDouble, NullOpDouble, AbsDiff, L1Minimiser<TNumber>>;

template <typename TNumber>
using L2Distance = TLpDistance<TNumber, SquareDouble, RootDouble, SquareDiff, L2Minimiser<TNumber>>;

template <typename TNumber>
using L_oo_Distance = TLpDistance<TNumber, NullOpDouble, NullOpDouble, MaxDiff, NoimplMinimiser>;





  
/* where TVDataIn has public size_t dimension and typename TVDataIn::Sample */
template<typename TVDataIn>
class LpMetric{
  
  public:
 //   typedef typename TVDataIn::Sample Sample;
    typedef typename TVDataIn::AtomicType AtomicType;
    typedef LpMetricInitializer Initializer;    
  
  private:
    std::unique_ptr<BaseLpDistance<AtomicType>>  uptr_lpdistance;

  public:
  
    char get_p(){
      return p;
    }
    
    LpMetric(const TVDataIn & datain, size_t nthreads, const LpMetricInitializer & l2mi):p(l2mi.p) 
      {
      
      //quelch warning.
      (void)nthreads;
      
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
    
    template <typename T, typename V>
    void set_distance(const T & a, const V & b, double threshold, double & distance){
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

    
    template <typename TPointerCenter, typename TFunction>
    void set_center(const std::string & energy, const TFunction & f_sample, size_t ndata, TPointerCenter ptr_center){
      uptr_lpdistance->set_center(energy, f_sample, ndata, ptr_center);
    }

    template <typename TPointerCenter, typename TFunction>
    void update_sum(const std::string & energy, const TFunction & f_add, size_t n_add, const TFunction & f_subtract, size_t n_subtract, TPointerCenter ptr_sum){
      uptr_lpdistance->set_center(energy, f_add, n_add, f_subtract, n_subtract, ptr_sum);
    }
  
  private:
    const char p;

};

} //namespace nszen

#endif

