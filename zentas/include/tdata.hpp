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
#ifndef ZENTAS_TDATA_H
#define ZENTAS_TDATA_H

/* 
 * 
 * define classes satisfying TData template

 TODO : the requirements on TData should be minimal. Should have another class wrapping the minimal class (with throwing etc).
 
template <typename TDataIn> 
class TData{
  public: 
    size_t ndata;
    TData (const TDataIn &, bool);
    void append(const TDataIn::Sample &);
    const TDataIn::Sample & at (size_t) const; //by reference optional...
    void remove_last ();
    void replace_with(size_t, const TDataIn::Sample &);
};

*/

//#include "outputwriter.hpp"

#include <sstream>
#include <vector>
#include "zentaserror.hpp"
#include "sparsevectorrfcenter.hpp"

namespace nszen{  



template <typename AtomicType>
class SparseRefinementCenterData{

 public:
  
  /* Rooted or unrooted sparse data-in */
  //typedef TSDataIn DataIn;
  //typedef typename DataIn::Sample Sample;
  
  //using AtomicType = typename DataIn::AtomicType;
  
  typedef SparseVectorRfCenter<AtomicType> RfCenter;
  
  SparseRefinementCenterData(size_t dimension){(void)dimension;} //maybe use dimension to initialise size of maps

  
  void add(size_t, const SparseVectorSample<AtomicType> & ){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY add");
  }

  void subtract(size_t, const SparseVectorSample<AtomicType> &){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY subtract");
  }
  
  void set_to_zero(size_t){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY set to zero");
  }

  void scale(size_t, const RfCenter &, double){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY scale");
  }

  const RfCenter & at_for_metric(size_t) const {
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY return at_for_metric");
  }

  RfCenter & at_for_change(size_t) {
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY return at_for_change");
  }
    
  void append_zero(){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY append_zero");
  }

  void set_zero(size_t){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY set_zero");
  }


  void replace_with(size_t j, const RfCenter & t){
    (void)j;
    (void)t;
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY replace_with");
  }

  bool equals (size_t j, const RfCenter & t){
    (void)j;
    (void)t;
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY equals");
  }

  std::string string_for_sample(size_t){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY string_for_sample");
  }
  
};


template <typename AtomicType>
class VCenterData{
    
 public:

  size_t ndata;    
  const size_t dimension;
  std::vector<AtomicType> data;

  VCenterData (size_t dimension_):dimension(dimension_){}

  void append(const AtomicType * const datapoint){
    ndata += 1;
    data.insert(data.end(), datapoint, datapoint + dimension);
  }

  void add(size_t i, const AtomicType * const datapoint){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] += *(datapoint + d);
    }
  }

  void subtract(size_t i, const AtomicType * const datapoint){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] -= *(datapoint + d);
    }
  }
  
  void set_to_zero(size_t i){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] = 0;
    }
  }


  void scale(size_t i, const AtomicType * const datapoint, double alpha){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] = *(datapoint + d)*alpha;            
    }
  }
  
  void append_zero(){
    data.resize((ndata + 1)*dimension, 0);
    ++ndata;
  }
  
  void set_zero(size_t i){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] = 0;
    }
  }

    std::string string_for_sample(size_t i){
      std::stringstream ss;
      ss << "(";
      for (size_t d = 0; d < dimension; ++d){
        ss << " " << data[i*dimension + d] << " ";
      }
      ss << ")";
      return ss.str();
    }

  const AtomicType * const at_for_metric(size_t j) const {
    return data.data() + j*dimension;
  }

  AtomicType * const at_for_change(size_t j) {
    return data.data() + j*dimension;
  }


  void replace_with(size_t j, const AtomicType * const datapoint){
    for (size_t d = 0; d < dimension; ++d){
      data[j*dimension + d] = *(datapoint + d);
    }
  }

  bool equals(size_t i, const AtomicType * const datapoint) const{
    double abs_sum = 0.; 
    double sum_abs = 0.;
    for (size_t d = 0; d < dimension; ++d){

      abs_sum += std::abs(data[i*dimension + d] -  *(datapoint + d));
      sum_abs += std::abs(data[i*dimension + d]) +  std::abs(*(datapoint + d));
    }
    
    double relative_change = (abs_sum / sum_abs);
    return relative_change < 1e-7;
  }


  
  
  
};


template <typename TDataIn>
class VData{
  
  public:
    using AtomicType = typename TDataIn::AtomicType;
    typedef VCenterData<AtomicType> RefinementCenterData;

  private:
    typedef typename TDataIn::Sample Sample;
        
  public:
    typedef TDataIn DataIn;
    typedef typename TDataIn::InitBundle InitBundle;    
    
  private:
    size_t ndata;    
    const size_t dimension;
    std::vector<AtomicType> data;

    
  public: 
    VData (const DataIn & datain, bool as_empty): ndata(0), dimension(datain.dimension){
      if (as_empty == false){
        throw zentas::zentas_error("Currently, there is no implementation for constructing VData from DenseVectorDataUnrootedIn with data, due to the lack of any apparent need");
      }
    }
    
     size_t get_ndata() const{
      return ndata;
    }
    
    void append(const Sample & datapoint){
      ndata += 1;
      data.insert(data.end(), datapoint, datapoint + dimension);
    }

        
    const Sample at_for_metric(size_t j) const {
      return data.data() + j*dimension;
    }

    const Sample at_for_move(size_t j) const {
      return data.data() + j*dimension;
    }
        
    void remove_last(){
      if (ndata == 0){
        throw zentas::zentas_error("Request for remove_last, but ndata == 0");
      }
      ndata -= 1;
      data.resize(ndata*dimension);
    }
    
    void replace_with(size_t j, const Sample & datapoint){
      for (size_t d = 0; d < dimension; ++d){
        data[j*dimension + d] = *(datapoint + d);
      }
    }

    
    std::string string_for_sample(size_t i){
      std::stringstream ss;
      ss << "(";
      for (size_t d = 0; d < dimension; ++d){
        ss << " " << data[i*dimension + d] << " ";
      }
      ss << ")";
      return ss.str();
    }
    
    std::string get_string(){
      std::stringstream ss;
      for (size_t i = 0; i < ndata; ++i){
        ss << string_for_sample(i);
      }
      
      return ss.str();
    }
};



template <typename TSDataIn>
class SData{
  

  public:
    typedef typename TSDataIn::Sample Sample;
    typedef typename TSDataIn::AtomicType AtomicType;
    typedef TSDataIn DataIn;

  public:
    //typedef NoRefinementPossible<Sample> RefinementCenterData;

    
    typedef typename TSDataIn::InitBundle InitBundle;
    
    
    

   private: 
    size_t ndata;   
    size_t dimension;
    std::vector<AtomicType> data;
    std::vector<size_t> sizes;
    
   public: 

    
    SData (const DataIn & datain, bool as_empty): ndata(0), dimension(datain.get_max_size()){  //, mowri(true, false, ""){
      if (as_empty == false){
        throw zentas::zentas_error("Currently, there is no implementation for constructing SData from SDataIn with data, due to the lack of any apparent need");
      }
    }

     size_t get_ndata() const{
      return ndata;
    }
    
    void append(const Sample & s){
      ndata += 1;
      data.insert(data.end(), s.values, s.values + s.size);
      data.resize(ndata*dimension);
      sizes.push_back(s.size);
    }
    
    const Sample at_for_metric(size_t j) const {
      return Sample(sizes[j], data.data() + j*dimension);
    }

    const Sample at_for_move(size_t j) const {
      return Sample(sizes[j], data.data() + j*dimension);
    }
        
    void remove_last(){
      if (ndata == 0){
        throw zentas::zentas_error("Request for remove_last, but ndata == 0");
      }
      ndata -= 1;
      data.resize(ndata*dimension);
      sizes.resize(ndata);
    }
    
    void replace_with(size_t j, const Sample & s){
      for (size_t d = 0; d < s.size; ++d){
        data[j*dimension + d] = *(s.values + d);
      }
      sizes[j] = s.size;
    }
    
    
    std::string string_for_sample(size_t i){
      std::stringstream ss;
      
      ss << "{";
      for (size_t d = 0; d < sizes[i]; ++d){
        ss << data[i*dimension + d];
      }
      ss << "}";
      return ss.str();
    }
    
    std::string get_string(){
      std::stringstream ss;
      for (size_t j = 0; j < ndata; ++j){
        ss << string_for_sample(j) << " ";
      }
      return ss.str();
    }

};



template <typename TSDataIn>
class SparseVectorData{
  

  public:
    typedef typename TSDataIn::Sample Sample;
    typedef typename TSDataIn::AtomicType AtomicType;
    typedef TSDataIn DataIn;
    typedef typename TSDataIn::InitBundle InitBundle;

  public:
    typedef SparseRefinementCenterData<AtomicType> RefinementCenterData;

   
   private: 
    size_t ndata;
    size_t dimension; //note that this dimension is the max non-sparsity, NOT the true vector-space dimension. 
    std::vector<AtomicType> data;
    std::vector<size_t> indices_s;
    std::vector<size_t> sizes;
    
   public: 
    SparseVectorData (const DataIn & datain, bool as_empty): ndata(0), dimension(datain.get_max_size()){
      if (as_empty == false){
        throw zentas::zentas_error("Currently, there is no implementation for constructing SData from SDataIn with data, due to the lack of any apparent need");
      }
    }

     size_t get_ndata() const{
      return ndata;
    }
    
    void append(const Sample & s){
      ndata += 1;
      data.insert(data.end(), s.values, s.values + s.size);
      data.resize(ndata*dimension);
      indices_s.insert(indices_s.end(), s.indices, s.indices + s.size);
      indices_s.resize(ndata*dimension);
      sizes.push_back(s.size);
    }
    
    const Sample at_for_metric(size_t j) const {
      return Sample(sizes[j], data.data() + j*dimension, indices_s.data() + j*dimension);
    }

    const Sample at_for_move(size_t j) const {
      return Sample(sizes[j], data.data() + j*dimension, indices_s.data() + j*dimension);
    }
        
    void remove_last(){
      if (ndata == 0){
        throw zentas::zentas_error("Request for remove_last, but ndata == 0");
      }
      ndata -= 1;
      data.resize(ndata*dimension);
      indices_s.resize(ndata*dimension);
      sizes.resize(ndata);
    }
    
    void replace_with(size_t j, const Sample & s){
      for (size_t d = 0; d < s.size; ++d){
        data[j*dimension + d] = *(s.values + d);
        indices_s[j*dimension + d] = *(s.indices + d);
      }
      sizes[j] = s.size;
    }


    std::string string_for_sample(size_t i){
      std::stringstream ss;
      for (size_t d = 0; d < sizes[i]; ++d){
        ss << "(" << data[i*dimension + d] << ")" << data[i*dimension + d];
      }
      return ss.str();
    }
        
    std::string get_string(){
      std::stringstream ss;
      for (size_t j = 0; j < ndata; ++j){
        ss << j << "   ";
        ss << string_for_sample(j);
        ss << "\n";
      }
      
      return ss.str();
    }
};











/* **************************************************
 * : Rooted data structures (non big memory copies) :
 * ************************************************* */

template <typename TDataIn>
class BaseDataRooted{
  
  private:
    typedef typename TDataIn::Sample Sample;
    typedef typename std::remove_const<typename std::remove_pointer<Sample>::type>::type NP_Sample;
  
  public:
    typedef TDataIn DataIn;
    typedef typename TDataIn::InitBundle InitBundle;
 
  protected:
    size_t ndata;    
    const size_t dimension;
    std::vector<size_t> IDs;
    const  TDataIn * const  ptr_datain;
  
  public: 
    BaseDataRooted (const DataIn & datain, bool as_empty): ndata(0), dimension(datain.dimension), ptr_datain(&datain){
      if (as_empty == false){
        throw zentas::zentas_error("Currently, there is no implementation for constructing VData from DenseVectorDataUnrootedIn with data, due to the lack of any apparent need");
      }
    }
    
    size_t get_ndata() const{
      return ndata;
    }
    void append(size_t i){
      ndata += 1;
      IDs.push_back(i);
    }

    size_t at_for_move(size_t j) const {
      return IDs[j];
    }
    void remove_last(){
      if (ndata == 0){
        throw zentas::zentas_error("Request for remove_last, but ndata == 0");
      }
      ndata -= 1;
      IDs.resize(ndata);
    }
    void replace_with(size_t j, size_t i){
      IDs[j] = i;
    }

    std::string string_for_sample(size_t i){
      (void)i;
      return "string_for_sample not implemented for vdatarooted";
    }    
    
    std::string get_string(){
      return "get_string not implemented for vdatarooted";
    }
};


template <typename TDataIn>
class VDataRooted : public  BaseDataRooted <TDataIn>{

  public:
    using AtomicType = typename TDataIn::AtomicType;
    typedef VCenterData<AtomicType> RefinementCenterData;

        
  public:
    using BaseDataRooted<TDataIn>::ptr_datain;
    using BaseDataRooted<TDataIn>::IDs;
    using BaseDataRooted<TDataIn>::dimension;

  public:
    typedef typename TDataIn::Sample Sample;
    
    typedef TDataIn DataIn;

    VDataRooted (const DataIn & datain, bool as_empty) : BaseDataRooted<TDataIn> (datain, as_empty) {}

    const Sample at_for_metric(size_t j) const {
      return ptr_datain->data + IDs[j]*dimension;
    }  
};



template <typename TDataIn>
class SDataRooted : public BaseDataRooted<TDataIn> {

    
  public: 
    using BaseDataRooted<TDataIn>::ptr_datain;
    using BaseDataRooted<TDataIn>::IDs;
    using BaseDataRooted<TDataIn>::dimension;

  public:
    typedef typename TDataIn::Sample Sample;
    typedef TDataIn DataIn;

    
 public:
 
    SDataRooted (const DataIn & datain, bool as_empty) : BaseDataRooted<TDataIn> (datain, as_empty) {}
    const Sample at_for_metric(size_t j) const {
      return Sample(ptr_datain->get_size(IDs[j]), ptr_datain->get_data(IDs[j]));
    }
  
};


template <typename TDataIn>
class SparseVectorDataRooted : public BaseDataRooted<TDataIn>{


  public: 
    using BaseDataRooted<TDataIn>::ptr_datain;
    using BaseDataRooted<TDataIn>::IDs;
    using BaseDataRooted<TDataIn>::dimension;

  public:
    typedef typename TDataIn::Sample Sample;
    typedef TDataIn DataIn;
    using AtomicType = typename TDataIn::AtomicType;

  public:
    
    typedef SparseRefinementCenterData<AtomicType> RefinementCenterData;
    
    SparseVectorDataRooted (const DataIn & datain, bool as_empty) : BaseDataRooted<TDataIn> (datain, as_empty) {}
    
    const Sample at_for_metric(size_t j) const {
      return Sample(ptr_datain->get_size(IDs[j]), ptr_datain->get_data(IDs[j]), ptr_datain->get_indices_s(IDs[j]));
    }
    
};



/* ********************************************
 * : general functions for manipulating TData :
 * ******************************************** */

template <class TData>
void replace_with_last(TData & data, size_t j){
  if (data.get_ndata() == 0){
    throw zentas::zentas_error("Request for replace_with_last, but ndata == 0");
  }
  
  else if (j == data.get_ndata() - 1){
    throw zentas::zentas_error("Request for replace_with_last(ndata - 1) which is valid but non-sensical");
  }
  
  data.replace_with(j, data.at_for_move(data.get_ndata() - 1));
}

template <class TData>
void swap(TData & data_1, size_t j1, TData & data_2, size_t j2){
  data_2.append(data_1.at_for_move(j1));
  data_1.replace_with(j1, data_2.at_for_move(j2));
  replace_with_last(data_2, j2);
  data_2.remove_last();
}



} //namespace nszen

#endif
