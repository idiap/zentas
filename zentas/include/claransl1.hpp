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
#ifndef ZENTAS_CLARANSL1_HPP
#define ZENTAS_CLARANSL1_HPP

#include "baseclarans.hpp"



namespace nszen{

 
template <class TMetric, class TData>
class ClaransL1 : public BaseClarans<TMetric, TData> {

  typedef typename TData::DataIn DataIn;
   
  private:
  
    std::unique_ptr<double []> up_dists_centers_old_k_to;
    double * const dists_centers_old_k_to;
    
    std::unique_ptr<double []> up_dists_centers_new_k_to;
    double * const dists_centers_new_k_to;    

    std::vector<XNearestInfo> center_nearest_center;


  public:
    
    
    ClaransL1(const BaseClustererInitBundle<DataIn,  TMetric> & ib, const BaseClaransInitBundle & clib): 
    BaseClarans<TMetric, TData> (ib, clib),
    up_dists_centers_old_k_to(new double [ib.K]), dists_centers_old_k_to(up_dists_centers_old_k_to.get()), up_dists_centers_new_k_to(new double [ib.K]), dists_centers_new_k_to(up_dists_centers_new_k_to.get()){}
        
  private: 

    /* ********************************************************************************************************
     * *************** no additional properties of samples, so default operators defined **********************
     * ********************************************************************************************************/
    virtual inline void put_sample_custom_in_cluster(size_t i, size_t k_first_nearest, const double * const distances) final override{
      put_nearest_2_infos_margin_in_cluster(i, k_first_nearest, distances);
    }
    
    virtual inline void reset_sample_custom(size_t k, size_t j, size_t nearest_center, const double * const distances) final override{
      reset_sample_nearest_2_infos_margin(k, j, nearest_center, distances);
    }

    virtual inline void custom_append(size_t k_to, size_t k, size_t j) final override{
      nearest_2_infos_margin_append(k_to, k, j);
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

    virtual void custom_info_test() override final {
      nearest_2_infos_margin_test();
    }

    virtual void put_sample_in_cluster(size_t i) override final {
      base_put_sample_in_cluster(i);
    }

    //virtual bool update_centers() override final{
      //return update_centers_greedy();
    //}

    /* ************************************************************
     * *************** center - center stuff **********************
     * ************************************************************/
    
     
    

    /* set center_nearest_center[k] for k in [0,K) */ 
    inline virtual void set_center_center_info() override final {
      set_center_center_info_l1(center_nearest_center);
    }
    

    
    
    virtual void update_center_center_info() override final {
      update_center_center_info_l1(center_nearest_center, dists_centers_old_k_to, dists_centers_new_k_to);
    }
    
          
    virtual void center_center_info_test() override final{
      center_center_info_test_l1(center_nearest_center);
    }
    
    
    
    virtual void custom_acceptance_call() override final {
      for (size_t k = 0; k < K; ++k){
        set_center_center_distance(k, k_to, dists_centers_old_k_to[k]);
        set_center_sample_distance(k, k_from, j_from, dists_centers_new_k_to[k]);
      }
    }
    
    

    /* *******************************************************************
     * *************** the 2 core clarans functions **********************
     * *******************************************************************/    


    virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial) override final{
      return get_delta_E_l1(k1, k2, j2, center_nearest_center[k1].d_x, serial);
    }


    virtual void update_sample_info() override final{
      update_sample_info_l1(dists_centers_old_k_to, dists_centers_new_k_to);
    }
    
    //virtual bool evaluate_proposal(size_t k1, size_t k2, size_t j2) override final{      
      //return evaluate_proposal_l1(k1, k2, j2, center_nearest_center[k1].d_x);
    //}

    

    
  public:

    using BaseClusterer<TMetric, TData>::K;
    using BaseClusterer<TMetric, TData>::set_center_sample_distance;
    using BaseClusterer<TMetric, TData>::set_center_center_distance;
    using BaseClusterer<TMetric, TData>::base_put_sample_in_cluster;    
    using BaseClarans<TMetric, TData>::k_to;
    using BaseClarans<TMetric, TData>::k_from;
    using BaseClarans<TMetric, TData>::j_from;
    using BaseClarans<TMetric, TData>::put_nearest_2_infos_margin_in_cluster;   
    using BaseClarans<TMetric, TData>::reset_sample_nearest_2_infos_margin;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_append;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_replace_with_last;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_replace_with;    
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_remove_last;
    using BaseClarans<TMetric, TData>::nearest_2_infos_margin_test;
    using BaseClarans<TMetric, TData>::update_sample_info_l1;
    using BaseClarans<TMetric, TData>::get_delta_E_l1;
    using BaseClarans<TMetric, TData>::center_center_info_test_l1;
    using BaseClarans<TMetric, TData>::update_center_center_info_l1;    
    using BaseClarans<TMetric, TData>::update_centers_greedy;      
    using BaseClarans<TMetric, TData>::set_center_center_info_l1;
  

};

} //namespace nszen
 

#endif


