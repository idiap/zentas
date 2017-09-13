// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_CLARANSL23_HPP
#define ZENTAS_CLARANSL23_HPP

#include "baseclarans.hpp"
  
namespace nszen{


class BaseClaransL23 : public BaseClarans {

   
  private:
  
    std::unique_ptr<double []> up_dists_centers_old_k_to;
    double * const dists_centers_old_k_to;
    
    std::unique_ptr<double []> up_cc;
    double * const cc; 

    std::unique_ptr<double []> up_d_min_cc;
    double * const d_min_cc;
    
    std::unique_ptr<size_t []> up_a_min_cc;
    size_t * const a_min_cc;
    
    


  public:
  
    BaseClaransL23(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb): BaseClarans(sb, eb),  
//    BaseClarans<TMetric, TData> (ib, clib), 
    up_dists_centers_old_k_to(new double [sb.K]), dists_centers_old_k_to(up_dists_centers_old_k_to.get()), up_cc(new double [sb.K*sb.K]), cc(up_cc.get()), up_d_min_cc(new double [sb.K*sb.K]), d_min_cc(up_d_min_cc.get()), up_a_min_cc(new size_t [sb.K*sb.K]), a_min_cc(up_a_min_cc.get()) {}
      
     double * get_cc(){
      return cc;
    }
    
     double get_d_min_cc(size_t k){
      return d_min_cc[k];
    }
    
    
  private: 

    /* ********************************************************************************************************
     * *************** no additional properties of samples, so default operators defined **********************
     * ********************************************************************************************************/
    virtual void put_sample_custom_in_cluster(size_t i, size_t k_first_nearest, const double * const distances) final override{
      put_nearest_2_infos_margin_in_cluster(i, k_first_nearest, distances);
    }
    
    
    virtual void custom_append(size_t k_to, size_t k, size_t j) final override{
      nearest_2_infos_margin_append(k_to, k, j);
    }
    
    virtual void custom_replace_with_last(size_t k, size_t j) final override{
      nearest_2_infos_margin_replace_with_last(k, j);
    }
    
    virtual void custom_replace_with(size_t k1, size_t j1, size_t k2, size_t j2) final override{
      nearest_2_infos_margin_replace_with(k1, j1, k2, j2);
    }
    
    virtual void custom_remove_last(size_t k) final override{
      nearest_2_infos_margin_remove_last(k);
    }

    virtual void custom_info_test() override final {
      nearest_2_infos_margin_test();
    }



    /* ************************************************************
     * *************** center - center stuff **********************
     * ************************************************************/
     
    /* set cc[k,kp] for all k,kp, and set min_cc */ 
     virtual void set_center_center_info() override final {
      set_center_center_info_l2(cc, d_min_cc, a_min_cc);
    }

    
    virtual void update_center_center_info() override final {
      //updating cc is handled in custom_acceptance_call, nothing to do here.
    }
      
          
    virtual void center_center_info_test() override final{
      for (size_t k1 = 0; k1 < K; ++k1){
        for (size_t k2 = 0; k2 < K; ++k2){
          double distance;
          set_center_center_distance_nothreshold(k1, k2, distance);
          if (std::abs(distance - cc[k1*K + k2]) > 1e-7*(std::abs(distance) + std::abs(cc[k1*K + k2]))){
            throw zentas::zentas_error("cc not correct in center_center_info_test");
          }
        }
      }  
    }
    
    
    virtual void custom_acceptance_call() override final {
      acceptance_call_l2(cc, dists_centers_old_k_to, d_min_cc, a_min_cc);
    }

    virtual void put_sample_in_cluster(size_t i) override final {
      triangular_put_sample_in_cluster(i, cc);
    }
    
    virtual void update_sample_info() override final{
      update_sample_info_l23(dists_centers_old_k_to, cc);
    }
};

} //namespace nszen
 

#endif
  
  

