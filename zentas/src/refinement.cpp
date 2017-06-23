#include "skeletonclusterer.hpp"

#include <tuple>
#include <memory>
namespace nszen{


void SkeletonClusterer::initialise_refinement(){
  

  /* put centers into clusters, initialise rf_center and old_rf_center */
  for (size_t k = 0; k < K; ++k){
    
    put_sample_in_cluster(center_IDs[k]);
    append_zero_to_rf_center_data();
    append_zero_to_old_rf_center_data();
  }
  
  /* do whatever else needs to be done to have refinement variables ready, also clear-up unnec. mem */
  custom_initialise_refinement();

  //TODO : verify if test
  
}



void SkeletonClusterer::run_refinement(){
 
  while (halt_refinement() == false){

    /* ************** *
    * UPDATE CENTERS *
    * ****************/
    set_t_ncalcs_update_centers_start();
    bool modified_centers = rf_update_centers();

    if (modified_centers == false){
      update_t_ncalcs_center_end();
      mowri << rf_get_round_summary() << zentas::Endl;
      break;
    }

    /* ************************* *
    * UPDATE CENTER CENTER INFO *
    * ***************************/
    rf_update_center_center_info();
    if (with_tests == true) rf_post_center_update_test();
    update_t_ncalcs_center_end();

    for (auto & x : to_leave_cluster){
      x.resize(0);
    }
    
    /* ****************** *
    * UPDATE SAMPLE INFO *
    * ********************/
    rf_update_sample_info();
    if (with_tests == true) rf_post_sample_update_test();
    update_t_ncalcs_sample_update_end();


    
    /* ************ *
    * REDISTRIBUTE *
    * **************/
    rf_redistribute();
    if (with_tests == true) rf_post_redistribute_test();        
    update_t_redistibute_end();

    
    /* ************************* *
    * UPDATE CLUSTER STATISTICS *
    * ***************************/
    rf_update_all_cluster_statistics();
    if (with_tests == true) rf_post_update_statistics_test();
    update_t_update_all_cluster_stats_end();


    if (with_tests == true) mowri << zentas::Endl;
    mowri << rf_get_round_summary() << zentas::Endl;
    ++round;
  }

}

bool SkeletonClusterer::halt_refinement(){
  /* TODO : make refinement stop criterion, pass in via constructor.  */
  return false;
}


bool SkeletonClusterer::rf_update_centers(){
  
  bool has_changed = false;
  double sum_abs;
  max_delta_C = 0;
  for (size_t k = 0; k < K; ++k){
    set_old_rf_center_data(k);
    /* metric dependent function, room for optimisation (l2 quadratic could track changing membership) */
    set_rf_center_data(k);
    set_delta_rf_new_and_old(k, std::numeric_limits<double>::max(), delta_C[k]);
    max_delta_C = std::max(max_delta_C, delta_C[k]);
    set_sum_abs_rf_new_and_old(k, sum_abs);    
    cluster_has_changed[k] = (delta_C[k] > 1e-7*sum_abs);
    has_changed = cluster_has_changed[k] == true ? true : has_changed;
  }

  return has_changed;
}


void SkeletonClusterer::rf_set_2_smallest(size_t k, size_t j, size_t & min_k1, double & min_d1, size_t & min_k2, double & min_d2){
  
  min_k1 = 0;
  min_k2 = 0;
  min_d1 = std::numeric_limits<double>::max();
  min_d2 = std::numeric_limits<double>::max();
  
  double worker_dist;
  for (size_t kp = 0; kp < K; ++kp){
    set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
    if (worker_dist < min_d2){
      if (worker_dist < min_d1){
        min_d2 = min_d1;
        min_k2 = min_k1;
        min_d1 = worker_dist;
        min_k1 = kp;
      }
      else{
        min_d2 = worker_dist;
        min_k2 = kp;        
      }
    }
  }
}

//Hamerly's algorithm. 
void SkeletonClusterer::rf_update_sample_info_standard(){

  size_t min_k1, min_k2;
  double min_d1, min_d2;
  
  for (size_t k = 0; k < K; ++k){
    for (size_t j = 0; j < get_ndata(k); ++j){
      rf_set_2_smallest(k, j, min_k1, min_d1, min_k2, min_d2);
      reset_nearest_info(k,j, min_k1, min_d1, f_energy(min_d1));
      lower_2[k][j] = min_d2;
      v_b[k][j] = min_k2;
    }
  }
}


        
//Hamerly's algorithm. 
void SkeletonClusterer::rf_update_sample_info_hamerly(){

  size_t min_k1, min_k2;
  double min_d1, min_d2;
  
  for (size_t k = 0; k < K; ++k){
    for (size_t j = 0; j < get_ndata(k); ++j){
      //update bounds.
      lower_2[k][j] -= max_delta_C;
      //reset nearest_info (we don't have bound u like Hamerly, auto set distance to prev nearest)
      set_rf_center_sample_distance_nothreshold(k, k, j, min_d1);
      reset_nearest_info(k, j, k, min_d1, f_energy(min_d1));
      
      if (lower_2[k][j] < nearest_1_infos[k][j].d_x){
        rf_set_2_smallest(k, j, min_k1, min_d1, min_k2, min_d2);
        reset_nearest_info(k,j, min_k1, min_d1, f_energy(min_d1));
        lower_2[k][j] = min_d2;
        v_b[k][j] = min_k2;
      }
    }
  }
}

//exponion-like algorithm. 
void SkeletonClusterer::rf_update_sample_info_exponion(){
  
  size_t min_k1, min_k2;
  double min_d1, min_d2;
  double dist_b;
  double worker_dist;  
  
  // compl_exp_rad refers to a new radius, 
  // which like the original exponion radius, 
  // guarantees nearest and second nearest 
  // (in ball around previous nearest).  
  double max_compl_exp_rad;

  // the mimimum distance to a center other than `self'. 
  double min_cc;
  
  // sort in ascending order by first (double) argument
  auto tupsorter = [](std::tuple<double, size_t> & x, std::tuple<double, size_t> & y){ return std::get<0>(x) < std::get<0>(y);};

  // how many samples fail the Hamerly test
  size_t n_to_reconsider;
  
  size_t max_cluster_size = 0;
  for (size_t k = 0; k < K; ++k){
    max_cluster_size = std::max(max_cluster_size, get_ndata(k));
  }
  
  std::vector<double> compl_exp_rad(max_cluster_size); 
  std::vector<size_t> to_reconsider(max_cluster_size);
  std::vector<double> dist_to_b(max_cluster_size);
  std::vector<double> dist_to_a(max_cluster_size);
  std::vector<std::tuple<double, size_t>> exprad_index (max_cluster_size); 
  std::vector<size_t> insertion_index (max_cluster_size);


  //version using unique ptrs, no memory zeroing.
  //std::unique_ptr<double []> pt_compl_exp_rad(new double[max_cluster_size]); 
  //auto compl_exp_rad = pt_compl_exp_rad.get();
  
  //std::unique_ptr<size_t []> pt_to_reconsider(new size_t[max_cluster_size]);
  //auto to_reconsider = pt_to_reconsider.get();
  
  //std::unique_ptr<double []> pt_dist_to_b(new double [max_cluster_size]);
  //auto dist_to_b = pt_dist_to_b.get();
  
  //std::unique_ptr<double []> pt_dist_to_a(new double [max_cluster_size]);
  //auto dist_to_a = pt_dist_to_a.get();
  
  //std::unique_ptr< std::tuple<double, size_t> [] > pt_exprad_index (new std::tuple<double, size_t> [max_cluster_size]); 
  //auto exprad_index = pt_exprad_index.get();
  
  //std::unique_ptr<size_t []> pt_insertion_index (new size_t [max_cluster_size]);
  //auto insertion_index = pt_insertion_index.get();



  // not all inter-centroid distances are required. 
  // these vectors only store those which are required.
  // (ie there exists a sample with a sufficiently large compl_exp_rad).  
  size_t n_cc_required;
  std::vector<std::tuple<double, size_t>> cc_required_tuple(K);
  std::vector<size_t> cc_required_indices(K);
  std::vector<double> cc_required_distances(K);

  

  for (size_t k = 0; k < K; ++k){

    n_to_reconsider = 0;
    n_cc_required = 0;  

    // determine which samples fail the outer hamerly test (to `reconsider') 
    // update those that pass the outer hamerly test as nec,
    // coupte compl_exp_rad for those that fail the outer hamerly test.
    n_to_reconsider = 0;
    max_compl_exp_rad = 0;
    
    for (size_t j = 0; j < get_ndata(k); ++j){
      lower_2[k][j] -= max_delta_C;
      set_rf_center_sample_distance_nothreshold(k, k, j, min_d1);
      
      if (lower_2[k][j] < min_d1){
        set_rf_center_sample_distance_nothreshold(v_b[k][j], k, j, dist_b);
        to_reconsider[n_to_reconsider] = j;
        dist_to_b[n_to_reconsider] = dist_b;
        dist_to_a[n_to_reconsider] = min_d1;
        compl_exp_rad[n_to_reconsider] = (1 + 1e-6)*(dist_b < min_d1 ? 2*min_d1 : min_d1 + dist_b);
        max_compl_exp_rad = std::max(max_compl_exp_rad, compl_exp_rad[n_to_reconsider]);        
        ++n_to_reconsider;
      }
      
      else{
        reset_nearest_info(k, j, k, min_d1, f_energy(min_d1));
      }
    }
    
    n_cc_required = 0;
    min_cc = std::numeric_limits<double>::max();
    for (size_t kp = 0; kp < K; ++kp){
      //threshold of max_compl_exp_rad
      set_rf_center_center_distance(k, kp, max_compl_exp_rad, worker_dist);
      if (kp != k && worker_dist < max_compl_exp_rad){
        cc_required_tuple[n_cc_required] = std::make_tuple(worker_dist, kp);
        ++n_cc_required;
        min_cc = std::min(min_cc, worker_dist);
      }
    }
    cc_required_tuple[n_cc_required] = std::make_tuple(std::numeric_limits<double>::max(), std::numeric_limits<size_t>::max());
    ++n_cc_required;    
    
    std::sort(cc_required_tuple.begin(), cc_required_tuple.begin() + n_cc_required, tupsorter);
    cc_required_distances.resize(n_cc_required);
    cc_required_indices.resize(n_cc_required);
    for (size_t jp = 0; jp < n_cc_required; ++jp){
      cc_required_distances[jp] = std::get<0>(cc_required_tuple[jp]);
      cc_required_indices[jp] = std::get<1>(cc_required_tuple[jp]);
    }
    
    
    
    

    for (size_t ji = 0; ji < n_to_reconsider; ++ji){
      // min of exponion radius as per paper and new compl_exp_rad
      std::get<0>(exprad_index[ji]) = std::min((1 + 1e-6)*(2*dist_to_a[ji] + min_cc), compl_exp_rad[ji]);
      std::get<1>(exprad_index[ji]) = ji;
    }
    
    
    std::sort (exprad_index.begin(), exprad_index.begin() + n_to_reconsider, tupsorter);



    size_t i_to_recon = 0;
    size_t i_cc = 0;
    while (i_to_recon < n_to_reconsider){
      
      double exponion_radius = std::get<0>(exprad_index[i_to_recon]);
      size_t ji = std::get<1>(exprad_index[i_to_recon]);
      while (i_to_recon < n_to_reconsider && cc_required_distances[i_cc] > exponion_radius){
        insertion_index[ji] = i_cc;
        ++i_to_recon;
        exponion_radius = std::get<0>(exprad_index[i_to_recon]);
        ji = std::get<1>(exprad_index[i_to_recon]);
      }
      ++i_cc;
    }

    for (size_t ji = 0; ji < n_to_reconsider; ++ji){
      size_t j = to_reconsider[ji];
      size_t upto = insertion_index[ji];
      
      min_k1 = k;
      min_k2 = v_b[k][j];
      min_d1 = dist_to_a[ji];
      min_d2 = dist_to_b[ji];      
      if (min_d1 > min_d2){
        std::swap(min_d1, min_d2);
        std::swap(min_k1, min_k2);
      }

      for (size_t kpi = 0; kpi < upto; ++kpi){
        size_t kp = cc_required_indices[kpi];
                
        if (kp != v_b[k][j]){
        set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
        if (worker_dist < min_d2){
          if (worker_dist < min_d1){
            min_d2 = min_d1;
            min_k2 = min_k1;
            min_d1 = worker_dist;
            min_k1 = kp;
          }
          else{
            min_d2 = worker_dist;
            min_k2 = kp;        
          }
        }
        }
      }
      reset_nearest_info(k,j, min_k1, min_d1, f_energy(min_d1));
      lower_2[k][j] = min_d2;
      v_b[k][j] = min_k2;
    }
  }
}



void SkeletonClusterer::rf_update_sample_info(){
  rf_update_sample_info_exponion();
}


void SkeletonClusterer::rf_redistribute(){


  size_t k_new;
  size_t j;
  size_t j_new;
  
  for (size_t k = 0; k < K; ++k){
    for (size_t ji = to_leave_cluster[k].size(); ji-- > 0;){
      j = to_leave_cluster[k][ji];
      k_new = nearest_1_infos[k][j].a_x;
      j_new = get_ndata(k_new);
      cluster_has_changed[k] = true;
      cluster_has_changed[k_new] = true;

      nearest_1_infos[k_new].push_back(nearest_1_infos[k][j]);
      sample_IDs[k_new].push_back(sample_IDs[k][j]);
      append_across(k_new, k, j);
      rf_custom_append(k_new, k, j);
    

      //this is `remove_with_tail_pull'
      if (j != get_ndata(k) - 1){
        nearest_1_infos[k][j] = *(nearest_1_infos[k].end() - 1);
        sample_IDs[k][j] = *(sample_IDs[k].end() - 1);
        replace_with_last_element(k,j);
        rf_custom_replace_with_last(k,j);
      }
      nearest_1_infos[k].pop_back();
      sample_IDs[k].pop_back();
      remove_last(k);
      rf_custom_remove_last(k);
    }
  }


}

void SkeletonClusterer::rf_custom_append(size_t k_new, size_t k, size_t j){
  lower_2[k_new].push_back(lower_2[k][j]);
  v_b[k_new].push_back(v_b[k][j]);
}

void SkeletonClusterer::rf_custom_replace_with_last(size_t k, size_t j){
  lower_2[k][j] = lower_2[k].back();
  v_b[k][j] = v_b[k].back();
}

void SkeletonClusterer::rf_custom_remove_last(size_t k){
  lower_2[k].pop_back();
  v_b[k].pop_back();
}


void SkeletonClusterer::rf_update_center_center_info(){
  //TODO
}




std::string SkeletonClusterer::rf_get_round_summary(){
  return get_round_summary();
}



void SkeletonClusterer::rf_post_center_update_test(){
  std::cout << "TODO : rf_post_center_update_test" << std::endl;
}

void SkeletonClusterer::rf_post_sample_update_test(){
  std::cout << "TODO : rf_post_sample_update_test" << std::endl;
}

void SkeletonClusterer::rf_post_redistribute_test(){
  std::cout << "TODO : rf_post_redistribute_test" << std::endl;
}

void SkeletonClusterer::rf_post_update_statistics_test(){
  std::cout << "TODO : rf_post_update_statistics_test" << std::endl;
}

// TODO : any custom cluster stats here ?
void SkeletonClusterer::rf_update_all_cluster_statistics(){
  old_E_total = E_total;
  E_total = 0;
  for (size_t k = 0; k < K; ++k){
    if (cluster_has_changed[k] == true){ 
      cluster_energies[k] = 0;
      for (size_t j = 0; j < get_ndata(k); ++j){
        cluster_energies[k] += nearest_1_infos[k][j].e_x;
      }
      cluster_mean_energies[k] = cluster_energies[k] / static_cast<double> (get_ndata(k));
      cluster_has_changed[k] = false;
    }
    E_total += cluster_energies[k];
  }
}
    
}

