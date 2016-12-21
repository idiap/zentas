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

#ifndef ZENTAS_CLARANS_HPP
#define ZENTAS_CLARANS_HPP



#include "baseclarans.hpp"



namespace nszen{

 
template <class TMetric, class TData>
class ClaransL0 : public BaseClarans<TMetric, TData> {

  typedef typename TData::DataIn DataIn;  
  
  public:
    
    
    ClaransL0(const BaseClustererInitBundle<DataIn, TMetric> & ib, const BaseClaransInitBundle & clib): 
    BaseClarans<TMetric, TData> (ib, clib) {} 
    
    //ClaransL0(size_t K, const DataIn & datain, const size_t * const center_indices_init, size_t max_proposals, double maxtime): BaseClarans<TMetric, TData>(K, datain, center_indices_init, max_proposals, maxtime){}
    
   private: 
   
    /* There are no custom quantities to
     * monitor for a Sample other than those 
     * defined in BaseClarans, use functions 
     * defined in BaseClarans */
    virtual inline void put_sample_custom_in_cluster(size_t i, size_t k_nearest, const double * const distances) final override{
      put_nearest_2_infos_margin_in_cluster(i, k_nearest, distances);
    }
    
    virtual inline void reset_sample_custom(size_t k, size_t j, size_t nearest_center, const double * const distances) final override{
      reset_sample_nearest_2_infos_margin(k, j, nearest_center, distances);
    }

    virtual inline void custom_append(size_t k_new, size_t k, size_t j) final override{
      nearest_2_infos_margin_append(k_new, k, j);
    }
    
    virtual inline void custom_replace_with_last(size_t k, size_t j) final override{
      nearest_2_infos_margin_replace_with_last(k, j);
    }

    virtual inline void custom_replace_with(size_t k1, size_t j1, size_t k2, size_t j2) final override{
      nearest_2_infos_margin_replace_with(k1, j1, k2, j2);
    }

    
    virtual inline void custom_remove_last(size_t k) final override{
      nearest_2_infos_margin_remove_last(k);
    }

    virtual void update_sample_info() override final{
      basic_clarans_update_sample_info();
    }

    
    
    /* We do not record any 
     * center-center variables */
    virtual void set_center_center_info() final override {}
    
    virtual void update_center_center_info() final override {}

    virtual void custom_acceptance_call() override final {}
    
    
    /* We use the basic evaluation 
     * and update defined in BaseClarans */
    //bool evaluate_proposal(size_t k1, size_t k2, size_t j2) override final{
      //return basic_clarans_evaluate_proposal(k1, k2, j2);
    //}
    
    virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial) override final{
      
      //quelch warning
      serial = ~serial;
      //"TODO" : implement a parallel version of get_delta_E_l0
      return get_delta_E_l0(k1, k2, j2);
    }
    

    virtual void custom_info_test() override final {
      nearest_2_infos_margin_test();
    }

    virtual void put_sample_in_cluster(size_t i) override final {
      //BaseClusterer<TMetric, TData>::base_put_sample_in_cluster(i); 
      base_put_sample_in_cluster(i);
    }

    //virtual bool update_centers() override final{
      //return update_centers_greedy();
    //}

  
  
   
   public: 
    using BaseClarans<TMetric, TData>::put_nearest_2_infos_margin_in_cluster;   
    using BaseClarans<TMetric, TData>::reset_sample_nearest_2_infos_margin;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_append;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_replace_with_last;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_replace_with;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_remove_last;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_test;
    using BaseClarans<TMetric, TData>::basic_clarans_update_sample_info;
    //using BaseClarans<TMetric, TData>::basic_clarans_evaluate_proposal;
    using BaseClarans<TMetric, TData>::update_centers_greedy;
    using BaseClarans<TMetric, TData>::get_delta_E_l0;
    using BaseClusterer<TMetric, TData>::base_put_sample_in_cluster;
   
};

} //namespace nszen
 

#endif
