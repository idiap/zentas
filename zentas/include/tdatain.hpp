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
#ifndef ZENTAS_TDATAIN_H
#define ZENTAS_TDATAIN_H

#include <iostream>
#include <memory>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "zentaserror.hpp"


/* 
 * 
 * define classes satisfying TDataIn template

template <typename T> 
class TDataIn{
  public: 
    typedef Sample;
    size_t ndata;
    const Sample & at (size_t) const;
};

*/



namespace nszen{

/* For vector data : array of size ndata x dimension */
template <typename TAtomic>
struct VariableLengthSample{ /* Used for strings */
  public:
    size_t size;
    const TAtomic * const values;
    VariableLengthSample(size_t size, const TAtomic * const values): size(size), values(values) {}
    std::string str() const{
      std::stringstream ss;
      for (size_t i = 0; i < size; ++i){
        ss << values[i];
      }
      return ss.str();   
    }
    
};


template <typename TAtomic>
struct SparseVectorSample{ /* Sparse vector sample */
  public:
    size_t size;
    const TAtomic * const values;
    const size_t * const indices;
    SparseVectorSample(size_t size, const TAtomic * const values, const size_t * const indices): size(size), values(values), indices(indices) {}
};


template <typename TAtomic>
struct VariableLengthInitBundle{
  public:
    size_t ndata;
    const size_t * const sizes;
    const TAtomic * const data;
    VariableLengthInitBundle(size_t ndata, const size_t * const sizes, const TAtomic * const data): ndata(ndata), sizes(sizes), data(data) {}
};




template <typename TAtomic>
struct BaseInitBundle{
  public:    
    size_t ndata;
    size_t dimension;
    const TAtomic * const data;
    BaseInitBundle(size_t ndata, size_t dimension, const TAtomic * const data): ndata(ndata), dimension(dimension), data(data) {}
};


template <typename TAtomic>
struct BaseDataIn{
  public:
    typedef TAtomic AtomicType;
    typedef const TAtomic * Sample;
    typedef BaseInitBundle<TAtomic> InitBundle;
    size_t ndata;
    size_t dimension;
    const Sample data;


  public:
    BaseDataIn(size_t ndata, size_t dimension, const TAtomic * const data): ndata(ndata), dimension(dimension), data(data) {}
    BaseDataIn(const InitBundle & ib): BaseDataIn(ib.ndata, ib.dimension, ib.data) {}
    
     size_t get_ndata() const {
      return ndata;
    }
    
    Sample at_for_metric(size_t i) const {
      return data + dimension*i;
    }
    

    std::string string_for_sample(size_t i) const{
      (void)i;
      return "currently no string function for BaseDataIn";          
    }
    
};


/* The basis for string and sparse vector input data : ndata samples, each of a distinct length */
template <typename TAtomic>
class BaseVarLengthDataIn {
  
  public:
    typedef VariableLengthInitBundle<TAtomic> InitBundle;
    typedef const VariableLengthSample<TAtomic> Sample;
    typedef TAtomic AtomicType;
    size_t ndata;
    /* a hack to get acceptance by metric (in l2 for sparse). TODO : redesign / clean-up */
    size_t dimension = std::numeric_limits<size_t>::max();
    
  protected:
    const size_t * const sizes;
    const TAtomic * const data;
    std::unique_ptr<size_t> up_c_sizes;
    size_t * c_sizes;
    size_t max_size;
    double mean_size;
   
  public:
    BaseVarLengthDataIn(size_t ndata, const size_t * const sizes, const TAtomic * const data):
    ndata(ndata), sizes(sizes), data(data), up_c_sizes(new size_t [ndata + 1]), max_size(0), mean_size(0){
      
      
      c_sizes = up_c_sizes.get();
      c_sizes[0] = 0;
      for (size_t i = 0; i < ndata; ++i){
        c_sizes[i+1] = c_sizes[i] + sizes[i];
        max_size = std::max(max_size, sizes[i]);
      }
      mean_size = c_sizes[ndata] / static_cast<double> (ndata);

      if (ndata == 0){
        throw zentas::zentas_error("ndata == 0 in BaseVarLengthDataIn, this is strange");
      }
    }
    
    BaseVarLengthDataIn(const InitBundle & ib): BaseVarLengthDataIn(ib.ndata, ib.sizes, ib.data) {}

     size_t get_ndata() const{
      return ndata;
    }
    

     size_t get_size(size_t i) const{
      return sizes[i];
    }
    
     const TAtomic * get_data(size_t i) const{
      return data + c_sizes[i];
    }
        
    size_t get_max_size() const{
      return max_size;
    }
    size_t get_mean_size() const{
      return mean_size;
    }

    std::string string_for_sample(size_t i) const{
      std::stringstream ss;
      ss << "{";
      for (size_t d = 0; d < sizes[i]; ++d){
        ss << data[c_sizes[i] + d];
      }
      ss << "}";
      return ss.str();          
    }  

};


template <typename TAtomic> 
struct DenseVectorDataUnrootedIn : public BaseDataIn<TAtomic> {
  
  public:
    using BaseDataIn<TAtomic>::data;
    using BaseDataIn<TAtomic>::dimension;
    typedef typename BaseDataIn<TAtomic>::InitBundle InitBundle;
    typedef typename BaseDataIn<TAtomic>::Sample Sample;

    
    DenseVectorDataUnrootedIn(const InitBundle & ib):BaseDataIn<TAtomic>(ib) {}
    
    const Sample at_for_move(size_t i) const {
      return data + dimension*i;
    }
};

template <typename TAtomic> 
struct DenseVectorDataRootedIn : public BaseDataIn<TAtomic> {
  
  public:
    typedef typename BaseDataIn<TAtomic>::InitBundle InitBundle;
    typedef typename BaseDataIn<TAtomic>::Sample Sample;

    
    DenseVectorDataRootedIn(const InitBundle & ib):BaseDataIn<TAtomic>(ib) {}
    
    size_t at_for_move(size_t i) const {
      return i;
    }  
};






template <typename TAtomic>
class BaseStringDataIn: public BaseVarLengthDataIn<TAtomic> {
  
  
  public:

    using BaseVarLengthDataIn<TAtomic>::sizes;
    using BaseVarLengthDataIn<TAtomic>::data;
    using BaseVarLengthDataIn<TAtomic>::c_sizes;
    
    typedef VariableLengthInitBundle<TAtomic> InitBundle;
    typedef const VariableLengthSample<TAtomic> Sample;
    typedef TAtomic AtomicType;

    
    BaseStringDataIn(const InitBundle & ib): BaseVarLengthDataIn<TAtomic>(ib) {}
        
  
    const Sample at_for_metric(size_t i) const{
      return Sample(sizes[i], data + c_sizes[i]);
    }
    

};

template <typename TAtomic>
struct SparseVectorDataUnrootedInInitBundle{
  
  public:
    size_t ndata;
    const size_t * const sizes;
    const TAtomic * const data;
    const size_t * const indices_s;
    
    SparseVectorDataUnrootedInInitBundle(size_t ndata, const size_t * const sizes, const TAtomic * const data, const size_t * const indices_s):ndata(ndata), sizes(sizes), data(data), indices_s(indices_s) {}
};



template <typename TAtomic>
class StringDataUnrootedIn : public BaseStringDataIn<TAtomic>{
  
  public:
    using BaseVarLengthDataIn<TAtomic>::max_size;
    using BaseVarLengthDataIn<TAtomic>::mean_size;
    using BaseVarLengthDataIn<TAtomic>::sizes;
    using BaseVarLengthDataIn<TAtomic>::data;
    using BaseVarLengthDataIn<TAtomic>::c_sizes;
    typedef typename BaseVarLengthDataIn<TAtomic>::Sample Sample;
    typedef typename BaseVarLengthDataIn<TAtomic>::InitBundle InitBundle;
    
  
    StringDataUnrootedIn(const InitBundle & ib): BaseStringDataIn<TAtomic>(ib) {
      if (static_cast<double>(max_size) > 5.*mean_size){
        
        std::stringstream ss;
        ss << "mean size : " << mean_size << "  max size : " << max_size << "\n" << "The max size is more than 5 times the mean size. This means a potentially huge waste of memory. Consider the rooted version (rooted = true).";
        throw zentas::zentas_error(ss.str());
      }
    }
    
    const Sample at_for_move(size_t i) const{
      return Sample(sizes[i], data + c_sizes[i]);
    }
};


template <typename TAtomic>
class StringDataRootedIn : public BaseStringDataIn<TAtomic>{
  
  public:
    typedef typename BaseVarLengthDataIn<TAtomic>::InitBundle InitBundle;
    StringDataRootedIn(const InitBundle & ib): BaseStringDataIn<TAtomic>(ib){
      
    }
    size_t at_for_move(size_t i) const{
      return i;
    }
};


/* For sparse vector data */
template <typename TAtomic>
class SparseVectorDataUnrootedInBase : public BaseVarLengthDataIn<TAtomic> {  
  
  public:
    typedef SparseVectorDataUnrootedInInitBundle<TAtomic> InitBundle;
    typedef const SparseVectorSample<TAtomic> Sample;
    typedef TAtomic AtomicType;  
    using BaseVarLengthDataIn<TAtomic>::sizes;
    using BaseVarLengthDataIn<TAtomic>::data;
    using BaseVarLengthDataIn<TAtomic>::c_sizes;
    
  protected:
    const size_t * const indices_s;

  public:
    SparseVectorDataUnrootedInBase(const InitBundle & ib): BaseVarLengthDataIn<TAtomic>(ib.ndata, ib.sizes, ib.data), indices_s(ib.indices_s) {}
        
     const size_t * get_indices_s(size_t i) const{
      return indices_s + c_sizes[i];
    }
    
    const Sample at_for_metric(size_t i) const{
      return Sample(sizes[i], data + c_sizes[i], indices_s + c_sizes[i]);
    }
};


template <typename TAtomic>
class SparseVectorDataUnrootedIn : public SparseVectorDataUnrootedInBase<TAtomic>{


  public:

    using SparseVectorDataUnrootedInBase<TAtomic>::max_size;
    using SparseVectorDataUnrootedInBase<TAtomic>::mean_size;
    using SparseVectorDataUnrootedInBase<TAtomic>::sizes;
    using SparseVectorDataUnrootedInBase<TAtomic>::data;
    using SparseVectorDataUnrootedInBase<TAtomic>::indices_s;
    using SparseVectorDataUnrootedInBase<TAtomic>::c_sizes;
    typedef typename SparseVectorDataUnrootedInBase<TAtomic>::Sample Sample;
    typedef typename SparseVectorDataUnrootedInBase<TAtomic>::InitBundle InitBundle;
  
    SparseVectorDataUnrootedIn(const InitBundle & ib): SparseVectorDataUnrootedInBase<TAtomic>(ib) {
      if (static_cast<double>(max_size) > 5.*mean_size){
        std::stringstream ss;
        ss << "mean size : " << mean_size << "  max size : " << max_size << "\n" << "The max size is more than 5 times the mean size. This means a potentially huge waste of memory. Consider the rooted version (rooted = true) if implemented.";
        throw zentas::zentas_error(ss.str());
      }
    }

    const Sample at_for_move(size_t i) const{
      return Sample(sizes[i], data + c_sizes[i], indices_s + c_sizes[i]);
    }
};



template <typename TAtomic>
class SparseVectorDataRootedIn : public SparseVectorDataUnrootedInBase<TAtomic>{
  
  public:
    typedef typename SparseVectorDataUnrootedInBase<TAtomic>::InitBundle InitBundle;
    
    SparseVectorDataRootedIn(const InitBundle & ib): SparseVectorDataUnrootedInBase<TAtomic>(ib){}
    
    size_t at_for_move(size_t i) const{
      return i;
    }
};



} //namespace nszen

#endif
