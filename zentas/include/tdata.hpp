// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

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
#include <unordered_map>
#include "zentaserror.hpp"
#include "sparsevectorrfcenter.hpp"

#include "tdatain.hpp"

namespace nszen{  



template <typename TAtomic>
class SparseRefinementCenterData{

  public:  
  using RfCenter = SparseVectorRfCenter<TAtomic>;

  private:
  std::vector<RfCenter> cdata;

  public:  
  
  SparseRefinementCenterData(size_t dimension){(void)dimension;} //maybe use dimension to initialise size of maps

  
  void add(size_t i, const SparseVectorSample<TAtomic> & svs){
    for (size_t d = 0; d < svs.size; ++d){
      if (cdata[i].count(svs.indices[d]) == 1){
        cdata[i][svs.indices[d]] += svs.values[d];
      }
      else{
        cdata[i][svs.indices[d]] = svs.values[d];
      }
      if (cdata[i].at(svs.indices[d]) == 0){
        cdata[i].erase(svs.indices[d]);
      }
    }
  }

  void subtract(size_t i, const SparseVectorSample<TAtomic> & svs){
    for (size_t d = 0; d < svs.size; ++d){
      if (cdata[i].count(svs.indices[d]) == 1){
        cdata[i][svs.indices[d]] -= svs.values[d];
      }
      else{
        cdata[i][svs.indices[d]] = -svs.values[d];
      }
      
      if (cdata[i].at(svs.indices[d]) == 0){
        cdata[i].erase(svs.indices[d]);
      }
    }
    
  }
  
  void set_to_zero(size_t){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY set to zero");
  }

  void set_as_scaled(size_t i, const RfCenter & rfc, double alpha){
    cdata[i] = rfc;
    for(auto iter = cdata[i].begin(); iter != cdata[i].end(); ++iter){
      iter->second *= alpha;
    }
  }

  const RfCenter & at_for_metric(size_t i) const {
    return cdata[i];
  }

  RfCenter * at_for_change(size_t) {
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY return at_for_change");
  }
    
  void append_zero(){
    cdata.emplace_back(std::unordered_map<size_t, TAtomic> {});
  }

  void set_zero(size_t){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY set_zero");
  }


  void replace_with(size_t j, const RfCenter & t){
    cdata[j] = t;
  }

  bool equals (size_t j, const RfCenter & t){
    (void)j;
    (void)t;
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY equals");
  }

  std::string string_for_sample(size_t){
    throw zentas::zentas_error("SparseRefinementCenterData cannot CURRENTLY string_for_sample");
  }
  
  void set_sum_abs(size_t i, double & sum_abs){
    sum_abs = 0;
    for(auto iter = cdata[i].cbegin(); iter != cdata[i].cend(); ++iter){
      sum_abs += std::abs(iter->second);
    }
  }
};


template <typename TAtomic>
class VCenterData{
    
 public:
  
  size_t ndata;    
  size_t dimension;
  std::vector<TAtomic> data;

  VCenterData (size_t dimension_):ndata(0), dimension(dimension_){
  }


  void add(size_t i, const TAtomic * const datapoint){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] += *(datapoint + d);
    }
  }

  void subtract(size_t i, const TAtomic * const datapoint){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] -= *(datapoint + d);
    }
  }
  
  void set_to_zero(size_t i){
    for (size_t d = 0; d < dimension; ++d){
      data[i*dimension + d] = 0;
    }
  }


  void set_as_scaled(size_t i, const TAtomic * const datapoint, double alpha){
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

  const TAtomic * const at_for_metric(size_t j) const {
    return data.data() + j*dimension;
  }

  TAtomic * const at_for_change(size_t j) {
    return data.data() + j*dimension;
  }


  void replace_with(size_t j, const TAtomic * const datapoint){
    for (size_t d = 0; d < dimension; ++d){
      data[j*dimension + d] = *(datapoint + d);
    }
  }

  bool equals(size_t i, const TAtomic * const datapoint) const{
    double abs_sum = 0.; 
    double sum_abs = 0.;
    for (size_t d = 0; d < dimension; ++d){
      abs_sum += std::abs(data[i*dimension + d] -  *(datapoint + d));
      sum_abs += std::abs(data[i*dimension + d]) +  std::abs(*(datapoint + d));
    }
    
    double relative_change = (abs_sum / sum_abs);
    return relative_change < 1e-7;
  }
  
  
  void set_sum_abs(size_t i, double & sum_abs){
    sum_abs = 0;
    for (size_t d = 0; d < dimension; ++d){
      sum_abs += std::abs(data[i*dimension + d]);
    }
  }
  
};




template <typename TDataIn>
class VData{

  // does this imply default move constructor copy constructor are generated?
  public:
    VData () = default; 
        
  public:
    using AtomicType = typename TDataIn::AtomicType;
    using RefinementCenterData = VCenterData<AtomicType>;
        
  public:
    typedef TDataIn DataIn;


    // how to make sure that the returned pointer is not used after this object is deleted?  
    //ConstLengthInitBundle self_as_tdatain;
    ConstLengthInitBundle<AtomicType> get_as_datain_ib(){
      return {ndata, dimension, data.data()};
    }

     

  private:
    size_t ndata;    
    size_t dimension;
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
    
    void append(const AtomicType * const datapoint){
      ndata += 1;
      data.insert(data.end(), datapoint, datapoint + dimension);
    }

        
    const AtomicType * const at_for_metric(size_t j) const {
      return data.data() + j*dimension;
    }

    const AtomicType * const at_for_move(size_t j) const {
      return data.data() + j*dimension;
    }
        
    void remove_last(){
      if (ndata == 0){
        throw zentas::zentas_error("Request for remove_last, but ndata == 0");
      }
      ndata -= 1;
      data.resize(ndata*dimension);
    }
    
    void replace_with(size_t j, const AtomicType * const datapoint){
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
    SData() = default;


  public:
    using Sample = typename TSDataIn::Sample;
    using AtomicType = typename TSDataIn::AtomicType;
    typedef TSDataIn DataIn;


   private: 
    size_t ndata;   
    size_t dimension;
    std::vector<AtomicType> data;
    std::vector<size_t> sizes;
    
   public:     
    SData (const DataIn & datain, bool as_empty): ndata(0), dimension(datain.get_max_size()){
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
    SparseVectorData() = default;


  

  public:
    using AtomicType = typename TSDataIn::AtomicType ;
    using DataIn = TSDataIn ;
    using RefinementCenterData = SparseRefinementCenterData<AtomicType> ;

    // how to make sure that the returned pointer is not used after this object is deleted?  
    SparseVectorDataInitBundle<AtomicType> get_as_datain_ib(){
      return {ndata, sizes.data(), data.data(), indices_s.data()};
    }

    
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
    
    void append(const SparseVectorSample<AtomicType> & s){
      ndata += 1;
      data.insert(data.end(), s.values, s.values + s.size);
      data.resize(ndata*dimension);
      indices_s.insert(indices_s.end(), s.indices, s.indices + s.size);
      indices_s.resize(ndata*dimension);
      sizes.push_back(s.size);
    }
    
    SparseVectorSample<AtomicType> at_for_metric(size_t j) const {
      return SparseVectorSample<AtomicType>(sizes[j], data.data() + j*dimension, indices_s.data() + j*dimension);
    }

    SparseVectorSample<AtomicType> at_for_move(size_t j) const {
      return SparseVectorSample<AtomicType>(sizes[j], data.data() + j*dimension, indices_s.data() + j*dimension);
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
    
    void replace_with(size_t j, const SparseVectorSample<AtomicType> & s){
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



  public:
    BaseDataRooted() = default;




    
  public:
    typedef TDataIn DataIn;
 
  protected:
    size_t ndata;    
    size_t dimension;
    std::vector<size_t> IDs;
    const  TDataIn *  ptr_datain;
  
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
    VDataRooted() = default;


    
  public:
    using AtomicType = typename TDataIn::AtomicType;
    typedef VCenterData<AtomicType> RefinementCenterData;

    ConstLengthInitBundle<AtomicType> get_as_datain_ib(){
      throw zentas::zentas_error("the function get_as_datain_ib needs to be implemented for rooted (dense) data (data needs to be made from indices) ");
    }

        
    using BaseDataRooted<TDataIn>::ptr_datain;
    using BaseDataRooted<TDataIn>::IDs;
    using BaseDataRooted<TDataIn>::dimension;

  public:
    using Sample = typename TDataIn::Sample;
    
    typedef TDataIn DataIn;

    VDataRooted (const DataIn & datain, bool as_empty) : BaseDataRooted<TDataIn> (datain, as_empty) {}

    const AtomicType * const at_for_metric(size_t j) const {
      return ptr_datain->data + IDs[j]*dimension;
    }  
};



template <typename TDataIn>
class SDataRooted : public BaseDataRooted<TDataIn> {


  public:
    SDataRooted() = default;
    
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
    SparseVectorDataRooted() = default;

  public:
    using BaseDataRooted<TDataIn>::ptr_datain;
    using BaseDataRooted<TDataIn>::IDs;
    using BaseDataRooted<TDataIn>::dimension;

  public:
    typedef TDataIn DataIn;
    using AtomicType = typename TDataIn::AtomicType;

    SparseVectorDataInitBundle<AtomicType> get_as_datain_ib(){
      throw zentas::zentas_error("the function get_as_datain_ib needs to be implemented for rooted (sparse) data (data needs to be made from indices) ");
    }


  public:
    
    typedef SparseRefinementCenterData<AtomicType> RefinementCenterData;
    
    SparseVectorDataRooted (const DataIn & datain, bool as_empty) : BaseDataRooted<TDataIn> (datain, as_empty) {}
    
    SparseVectorSample<AtomicType> at_for_metric(size_t j) const {
      return SparseVectorSample<AtomicType>(ptr_datain->get_size(IDs[j]), ptr_datain->get_data(IDs[j]), ptr_datain->get_indices_s(IDs[j]));
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
