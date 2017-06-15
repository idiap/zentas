#include "voronoil0.hpp"

namespace nszen{

    
VoronoiL0::VoronoiL0(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb): SkeletonClusterer(sb) {(void)eb;}
      
void VoronoiL0::set_redistribute_order(std::vector<size_t> & redistribute_order) {
  std::iota(redistribute_order.begin(), redistribute_order.end(), 0);      
}

void VoronoiL0::put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) {

  (void)k1;
  (void)k2;
  (void)d2;
  (void)e2;
}

    
std::string VoronoiL0::get_round_summary()  {
  std::stringstream ss;
  ss << get_base_summary_string();
  return ss.str();
}

void VoronoiL0::update_sample_info() {
  double adistance;
  double min_distance;
  size_t min_k;
  for (size_t k = 0; k < K; ++k){
    for (size_t j = 0; j < get_ndata(k); ++j){
      min_distance = std::numeric_limits<double>::max();
      for (size_t kp = 0; kp < K; ++kp){
        set_center_sample_distance(kp, k, j, min_distance, adistance);
        if (adistance < min_distance){
          min_distance = adistance;
          min_k = kp;
        }
      }
      reset_nearest_info(k,j, min_k, min_distance, f_energy(min_distance));
    }
  }
}


bool VoronoiL0::update_centers() {
  
  bool modified = false;
  double E_old; 
  double E_prop;
  double E_prop_best;
  size_t j_prop_best;
  
  for (size_t k = 0; k < K; ++k){

    E_old = get_cluster_energy(k);
    
    // compute energies with all other centers
    E_prop_best = std::numeric_limits<double>::max();
    j_prop_best = 0;
    
    for (size_t j_prop = 0; j_prop < get_ndata(k); ++j_prop){
      E_prop = get_e1(k, j_prop); 
      for (size_t j = 0; j < get_ndata(k); ++j){
        E_prop += f_energy(get_sample_sample_distance_nothreshold(k, j_prop, j));
      }
      if (E_prop < E_prop_best){
        E_prop_best = E_prop;
        j_prop_best = j_prop;
      }
    }
    
    if (E_old > E_prop_best){
      modified = true;
      swap_center_with_sample(k, j_prop_best);
    }
  }
  
  return modified;
}


}
