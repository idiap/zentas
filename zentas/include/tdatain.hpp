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

// TODO rename this file.
#include "sparsevectorrfcenter.hpp"


/* ******************************************** 
* define classes satisfying TDataIn template  *
***********************************************

template <typename T> 
class TDataIn{
  public: 
    typedef Sample;
    size_t ndata;
    const Sample & at (size_t) const;
};

***********************************************/



namespace nszen{

template <typename TAtomic>
struct VariableLengthSample{
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
struct VariableLengthInitBundle{
  public:
    size_t ndata;
    const size_t * const sizes;
    const TAtomic * const data;
    VariableLengthInitBundle(size_t ndata, const size_t * const sizes, const TAtomic * const data): ndata(ndata), sizes(sizes), data(data) {}
};

template <typename TAtomic>
struct ConstLengthInitBundle{
  public:    
    size_t ndata;
    size_t dimension;
    const TAtomic * const data;
    ConstLengthInitBundle(size_t ndata, size_t dimension, const TAtomic * const data): ndata(ndata), dimension(dimension), data(data) {}
};


template <typename TAtomic>
struct SparseVectorDataInitBundle : public VariableLengthInitBundle<TAtomic> {
  public:
    const size_t * const indices_s;
    SparseVectorDataInitBundle(size_t ndata, const size_t * const sizes, const TAtomic * const data, const size_t * const indices_s):
    VariableLengthInitBundle<TAtomic>(ndata, sizes, data), indices_s(indices_s) {}
};


template <typename TAtomic>
struct BaseDataIn{
  public:
    
    /* used downstream ? */
    using AtomicType = TAtomic;
    
    size_t ndata;
    /* dimension takes on a different meaning, 
     * for ConstLength (true dimension) and 
     * VarLength (maximum index)  */
    size_t dimension;
    const TAtomic * data;

    BaseDataIn(size_t ndata, size_t dimension, const TAtomic * const data): ndata(ndata), dimension(dimension), data(data) {}
    BaseDataIn(const ConstLengthInitBundle<TAtomic> & ib): BaseDataIn(ib.ndata, ib.dimension, ib.data) {}

    BaseDataIn() = default;
    
    size_t get_ndata() const {
      return ndata;
    }
    
    virtual std::string string_for_sample(size_t i) const{
      std::stringstream ss;
      ss << "currently no string function, this is from BaseDataIn. To display sample " << i << " requires that the virtual method be specialised";
      return ss.str();
    }
};


/* The basis for dense vector data */
template <typename TAtomic>
struct BaseConstLengthDataIn : public BaseDataIn<TAtomic>{
  public:
  
    using Sample = const TAtomic * const;

  public:
    BaseConstLengthDataIn(const ConstLengthInitBundle<TAtomic> & ib): BaseDataIn<TAtomic>(ib) {}
    const TAtomic * const at_for_metric(size_t i) const {
      return BaseConstLengthDataIn::data + BaseConstLengthDataIn::dimension*i;
    }
    
    BaseConstLengthDataIn() = default;
    
};

template <typename TAtomic> 
struct DenseVectorDataUnrootedIn : public BaseConstLengthDataIn<TAtomic> {
  public:
    DenseVectorDataUnrootedIn(const ConstLengthInitBundle<TAtomic> & ib):BaseConstLengthDataIn<TAtomic>(ib) {}    
    const TAtomic * const at_for_move(size_t i) const {
      return BaseConstLengthDataIn<TAtomic>::data + BaseConstLengthDataIn<TAtomic>::dimension*i;
    }
    
    DenseVectorDataUnrootedIn() = default;
};

template <typename TAtomic> 
struct DenseVectorDataRootedIn : public BaseConstLengthDataIn<TAtomic> {
  public:    
    DenseVectorDataRootedIn(const ConstLengthInitBundle<TAtomic> & ib):BaseConstLengthDataIn<TAtomic>(ib) {}
    size_t at_for_move(size_t i) const {
      return i;
    }
    
    DenseVectorDataRootedIn() = default; 
};




/* The basis for string and sparse vector input data  */
template <typename TAtomic>
class BaseVarLengthDataIn : public BaseDataIn<TAtomic>{
  
  public:
    BaseVarLengthDataIn() = default; 
  
  protected:
    const size_t * sizes;
    /* changed from up to vector, as easier default copying */
    //std::unique_ptr<size_t []> up_c_sizes;
    std::vector<size_t> v_c_sizes;
    size_t * c_sizes;
    size_t max_size;
    double mean_size;
   
  public:
    BaseVarLengthDataIn(size_t ndata, const size_t * const sizes, const TAtomic * const data): BaseDataIn<TAtomic>(ndata, std::numeric_limits<size_t>::max(), data), sizes(sizes), v_c_sizes(ndata + 1,0), max_size(0), mean_size(0) {
      c_sizes = v_c_sizes.data();
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
    
    BaseVarLengthDataIn(const VariableLengthInitBundle<TAtomic> & ib): BaseVarLengthDataIn(ib.ndata, ib.sizes, ib.data) {}

     size_t get_size(size_t i) const{
      return sizes[i];
    }
    
     const TAtomic * get_data(size_t i) const{
      return BaseDataIn<TAtomic>::data + c_sizes[i];
    }
        

    std::string string_for_sample(size_t i) const{
      std::stringstream ss;
      ss << "{";
      for (size_t d = 0; d < sizes[i]; ++d){
        ss << BaseDataIn<TAtomic>::data[c_sizes[i] + d];
      }
      ss << "}";
      return ss.str();          
    }
    
    void unrooted_viability_test(){
      if (static_cast<double>(max_size) > 5.*mean_size){
        std::stringstream ss;
        ss << "In unrooted_viability_test, with mean size=" << mean_size << " and max size=" << max_size << ". The max size is more than 5 times larger than mean size. This means a potentially large waste of memory, and `cold cache'. Consider the rooted version of this algorithm. (rooted = true).";
        throw zentas::zentas_error(ss.str());
      }
    }
  
    size_t get_max_size() const{
      return max_size;
    }
    size_t get_mean_size() const{
      return mean_size;
    }

};


template <typename TAtomic>
class BaseStringDataIn: public BaseVarLengthDataIn<TAtomic> {
  
  
  public:

    using Sample = VariableLengthSample<TAtomic>;
    
    BaseStringDataIn(const VariableLengthInitBundle<TAtomic> & ib): BaseVarLengthDataIn<TAtomic>(ib) {}
    const Sample at_for_metric(size_t i) const{
      return Sample(BaseVarLengthDataIn<TAtomic>::sizes[i], BaseDataIn<TAtomic>::data + BaseVarLengthDataIn<TAtomic>::c_sizes[i]);
    }
    
    BaseStringDataIn() = default;
};




template <typename TAtomic>
class StringDataUnrootedIn : public BaseStringDataIn<TAtomic>{
  
  public:
    StringDataUnrootedIn(const VariableLengthInitBundle<TAtomic> & ib): BaseStringDataIn<TAtomic>(ib) { BaseVarLengthDataIn<TAtomic>::unrooted_viability_test(); }
    const VariableLengthSample<TAtomic> at_for_move(size_t i) const{
      return VariableLengthSample<TAtomic>(BaseVarLengthDataIn<TAtomic>::sizes[i], BaseDataIn<TAtomic>::data + BaseVarLengthDataIn<TAtomic>::c_sizes[i]);
    }
    
    StringDataUnrootedIn() = default;
};


template <typename TAtomic>
class StringDataRootedIn : public BaseStringDataIn<TAtomic>{
  
  public:

    StringDataRootedIn(const VariableLengthInitBundle<TAtomic> & ib): BaseStringDataIn<TAtomic>(ib){
      
    }
    size_t at_for_move(size_t i) const{
      return i;
    }
    
    StringDataRootedIn() = default;
};


/* For sparse vector data */
template <typename TAtomic>
class SparseVectorDataInBase : public BaseVarLengthDataIn<TAtomic> {  
  
  public:
    using Sample = SparseVectorSample<TAtomic>; 
    
  protected:
    const size_t * indices_s;

  public:
    SparseVectorDataInBase(const SparseVectorDataInitBundle<TAtomic> & ib): BaseVarLengthDataIn<TAtomic>(ib.ndata, ib.sizes, ib.data), indices_s(ib.indices_s) {}
        
     const size_t * get_indices_s(size_t i) const{
      return indices_s + BaseVarLengthDataIn<TAtomic>::c_sizes[i];
    }
    
    const SparseVectorSample<TAtomic> at_for_metric(size_t i) const{
      return SparseVectorSample<TAtomic>(
      BaseVarLengthDataIn<TAtomic>::sizes[i], 
      BaseDataIn<TAtomic>::data + BaseVarLengthDataIn<TAtomic>::c_sizes[i], 
      indices_s + BaseVarLengthDataIn<TAtomic>::c_sizes[i]);
    }
    
    SparseVectorDataInBase() = default;
};


template <typename TAtomic>
class SparseVectorDataUnrootedIn : public SparseVectorDataInBase<TAtomic>{
  
  public:  
    SparseVectorDataUnrootedIn(const SparseVectorDataInitBundle<TAtomic> & ib): 
    SparseVectorDataInBase<TAtomic>(ib) { BaseVarLengthDataIn<TAtomic>::unrooted_viability_test(); }

    SparseVectorSample<TAtomic> at_for_move(size_t i) const{
      return SparseVectorSample<TAtomic>(
      BaseVarLengthDataIn<TAtomic>::sizes[i], 
      BaseDataIn<TAtomic>::data + BaseVarLengthDataIn<TAtomic>::c_sizes[i], 
      SparseVectorDataInBase<TAtomic>::indices_s + BaseVarLengthDataIn<TAtomic>::c_sizes[i]);
    }
    
    SparseVectorDataUnrootedIn() = default;
    
    
};



template <typename TAtomic>
class SparseVectorDataRootedIn : public SparseVectorDataInBase<TAtomic>{
  
  public:
    SparseVectorDataRootedIn(const SparseVectorDataInitBundle<TAtomic> & ib): SparseVectorDataInBase<TAtomic>(ib){}    
    size_t at_for_move(size_t i) const{
      return i;
    }
    
    SparseVectorDataRootedIn() = default;
};



} //namespace nszen

#endif
