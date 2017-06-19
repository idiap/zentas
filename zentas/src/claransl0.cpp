#include "claransl0.hpp"

namespace nszen{

    
ClaransL0::ClaransL0(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb): BaseClarans (sb, eb) {} 
    
void ClaransL0::put_sample_custom_in_cluster(size_t i, size_t k_nearest, const double * const distances){
  put_nearest_2_infos_margin_in_cluster(i, k_nearest, distances);
}
    
void ClaransL0::reset_sample_custom(size_t k, size_t j, size_t nearest_center, const double * const distances){
  reset_sample_nearest_2_infos_margin(k, j, nearest_center, distances);
}

void ClaransL0::custom_append(size_t k_new, size_t k, size_t j){
  nearest_2_infos_margin_append(k_new, k, j);
}

void ClaransL0::custom_replace_with_last(size_t k, size_t j){
  nearest_2_infos_margin_replace_with_last(k, j);
}

void ClaransL0::custom_replace_with(size_t k1, size_t j1, size_t k2, size_t j2){
  nearest_2_infos_margin_replace_with(k1, j1, k2, j2);
}

void ClaransL0::custom_remove_last(size_t k){
  nearest_2_infos_margin_remove_last(k);
}

void ClaransL0::update_sample_info(){
  basic_clarans_update_sample_info();
}

/* We do not record any 
 * center-center variables */
void ClaransL0::set_center_center_info() {}

void ClaransL0::update_center_center_info() {}

void ClaransL0::custom_acceptance_call() {}

double ClaransL0::get_delta_E(size_t k1, size_t k2, size_t j2, bool serial){
  
  (void)serial;
  /* I am not going to implement a parallel version of get_delta_E_l0, l0 is not important */
  return get_delta_E_l0(k1, k2, j2);
}


void ClaransL0::custom_info_test() {
  nearest_2_infos_margin_test();
}

void ClaransL0::put_sample_in_cluster(size_t i) {
  base_put_sample_in_cluster(i);
}

}
