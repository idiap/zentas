#include "skeletonclusterer.hpp"

#include <tuple>
#include <memory>
namespace nszen{


enum class RfAlg {Exponion = 0, Yinyang = 1};


void SkeletonClusterer::custom_swap(size_t k1, size_t k2){
  std::swap(lower_2[k1], lower_2[k2]);
  std::swap(upper_1[k1], upper_1[k2]);
  std::swap(v_b[k1], v_b[k2]);
}


void SkeletonClusterer::rf_swap(size_t k1, size_t k2){
  swap_data(k1, k2);
  std::swap(nearest_1_infos[k1], nearest_1_infos[k2]);
  std::swap(sample_IDs[k1], sample_IDs[k2]);  
  custom_swap(k1, k2); //alg specific. 
}



void SkeletonClusterer::initialise_refinement(){

  RfAlg rf_alg = RfAlg::Exponion;
  
  if (rf_alg == RfAlg::Exponion){
    rf_n_groups = K/10;
  }

  else if (rf_alg == RfAlg::Yinyang){
    // K/10 if memory will not exceed 2 GB, 
    double max_ngbs = 2.0;
    // otherwise (1e9 / (ndata*sizeof(double)))
    rf_n_groups = std::min<size_t>(K/10,  static_cast<size_t>(max_ngbs*1e9 / ((sizeof(size_t) + sizeof(double))*ndata)));
  }
  
  else{
    //error.
  }

  /* put centers into clusters, initialise rf_center and old_rf_center */
  for (size_t k = 0; k < K; ++k){
    put_sample_in_cluster(center_IDs[k]);
    append_zero_to_rf_center_data();
    append_zero_to_rf_sum_data();    
    append_zero_to_old_rf_center_data();
  }

  custom_initialise_refinement();
  
  
  // cluster the centers
  rf_groups = get_subclustered_centers_labels(rf_n_groups);
  if (rf_groups.size() != K){
    throw zentas::zentas_error("groups.size() != K just after center clustering, this is strange (wrong), logic error");
  }
  
  if (*std::max_element(rf_groups.begin(), rf_groups.end()) >= rf_n_groups){
    throw zentas::zentas_error("max group index too large, logic error");
  }

  n_in_group.resize(rf_n_groups, 0);
  for (size_t k = 0; k < K; ++k){
    ++n_in_group[rf_groups[k]];
  }

  cum_in_group.resize(rf_n_groups + 1, 0);
  for (size_t g = 0; g < rf_n_groups; ++g){
    cum_in_group[g+1] = cum_in_group[g] + n_in_group[g];
  }

  // sort the clusters by group
  size_t p_in = 0;
  size_t p_out = 0;
  
  for (size_t g = 0; g < rf_n_groups; ++g){

    if (p_in != cum_in_group[g]){
      throw zentas::zentas_error("This does not make sense, p_in should be auto here");
    }
    
    p_out = cum_in_group[g+1];
    
    while (p_in < cum_in_group[g] || p_out < K){
      if (rf_groups[p_in] == g){
        ++p_in;
      }
      
      else if (rf_groups[p_out] != g){
        ++p_out;
      }
      
      else {
        rf_swap(p_in, p_out);
        std::swap(rf_groups[p_in], rf_groups[p_out]);
      }
    }
  }
  
    /* do whatever else needs to be done to have refinement variables ready, also clear-up unnec. mem */
  custom_rf_clear_initmem();

  
  max_delta_group.resize(rf_n_groups);
  l_gC.resize(K*rf_n_groups, std::numeric_limits<double>::max());
  for (size_t k = 0; k < K; ++k){
    //initialise sum
    for (size_t j = 0; j < get_ndata(k); ++j){
      rf_increment_sum(k, j);
    }
    
    //initialise group bounds (for centers)
    double adist;
    for (size_t g = 0; g < rf_n_groups; ++g){
      for (size_t kp = cum_in_group[g]; kp < cum_in_group[g+1]; ++kp){
        set_rf_center_center_distance(k, kp, l_gC[k*rf_n_groups + g], adist);
        l_gC[k*rf_n_groups + g] = std::min(l_gC[k*rf_n_groups + g], adist);
      } 
    }    
  }
  
  rf_round = 1;
  
  // maintaining historical upper bounds since time immemorial
  rf_cum_u_delta_gC.resize(rf_n_groups, 0);  
  rf_cum_u_delta_gC.reserve(200*rf_n_groups); //expect memory for at least 200 rounds.   

  // initialise group bounds (for data)
  rf_l_groups.resize(ndata*rf_n_groups);
  rf_t_groups.resize(ndata*rf_n_groups);
  
  size_t running_ID = 0;
  rf_glt_ID.resize(K);
  for (size_t k = 0; k < K; ++k){
    rf_glt_ID[k].resize(get_ndata(k));
    for (size_t j = 0; j < get_ndata(k); ++j){
      rf_glt_ID[k][j] = running_ID;
      ++running_ID;
    }
  }
}
  


void SkeletonClusterer::rf_tighten_nearest(size_t k){
  double min_d1;
  for (size_t j = 0; j < get_ndata(k); ++j){
    set_rf_center_sample_distance_nothreshold(k, k, j, min_d1);
    reset_nearest_info(k,j, k, min_d1, f_energy(min_d1));
    upper_1[k][j] = min_d1;
  }
}

void SkeletonClusterer::run_refinement(){
 
  while (halt_refinement() == false){

    /* ************** *
    * UPDATE CENTERS *
    * ****************/
    set_t_ncalcs_update_centers_start();
    bool modified_centers = rf_update_centers();

    // no change, prepare to break. 
    if (modified_centers == false){
      
      for (size_t k = 0; k < K; ++k){
        rf_tighten_nearest(k);
      }
      
      rf_update_energies();
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
    if (is_rf_correct_d1_round() == true){
      rf_update_energies();
      if (with_tests == true) rf_post_update_statistics_test();
    }
    else{
      E_total = -1;    
    }
    update_t_update_all_cluster_stats_end();


    if (with_tests == true) mowri << zentas::Endl;
    mowri << rf_get_round_summary() << zentas::Endl;
    ++round;
    ++rf_round;
  }
  


}

bool SkeletonClusterer::is_rf_correct_d1_round(){
  return (round < 5  || round%40 == 0);
}


bool SkeletonClusterer::is_rf_tighten_cluster_radius_round(){
  return (round < 5 || round % 15 == 0);
}


bool SkeletonClusterer::halt_refinement(){
  /* TODO : make refinement stop criterion, pass in via constructor.  */
  return false;//round > 3;
}


bool SkeletonClusterer::rf_update_centers(){

  bool has_changed = false;
  double sum_abs;
  
  for (size_t g = 0; g < rf_n_groups; ++g){
    max_delta_group[g] = 0;
    for (size_t k = cum_in_group[g]; k < cum_in_group[g+1]; ++k){
      set_old_rf_center_data(k);
      set_rf_center_data(k);
      set_delta_rf_new_and_old(k, std::numeric_limits<double>::max(), delta_C[k]);
      max_delta_group[g] = std::max(max_delta_group[g], delta_C[k]);
      set_sum_abs_rf_new_and_old(k, sum_abs);    
      cluster_has_changed[k] = (delta_C[k] > 1e-7*sum_abs);
      has_changed = cluster_has_changed[k] == true ? true : has_changed;
    }
  }

  rf_cum_u_delta_gC.resize(rf_n_groups*(rf_round + 1), 0);  
  for (size_t rf_or = 0; rf_or < rf_round; ++rf_or){
    for (size_t g = 0; g < rf_n_groups; ++g){
      rf_cum_u_delta_gC[rf_or*rf_n_groups + g] += max_delta_group[g];
    }
  }

  return has_changed;
}




// An exponion-LIKE algorithm
void SkeletonClusterer::rf_update_sample_info_exponion(){
  
  size_t min_k1, min_k2;
  double min_d1, min_d2;
  double worker_dist;  
  
  // the mimimum distance to a center other than `self'. 
  double min_cc;
  
  // sort in ascending order by first (double) argument
  auto tupsorter = [](std::tuple<double, size_t> & x, std::tuple<double, size_t> & y){ return std::get<0>(x) < std::get<0>(y);};

  //  not all inter-centroid distances are required. 
  //  these vectors only store those which are required.
  //  ie there exists a sample with a sufficiently large compl_exp_rad.  
  size_t n_cc_required;
  std::vector<double> cc_all(K);
  std::vector<std::tuple<double, size_t>> cc_required_tuple(K);
  std::vector<size_t> cc_required_indices(K);
  std::vector<double> cc_required_distances(K);


  for (size_t k = 0; k < K; ++k){
    
    // update upper bound on distance to furthest member
    if (is_rf_tighten_cluster_radius_round() == true){
      u1_C[k] = 0;
      for (size_t j = 0; j < get_ndata(k); ++j){
        u1_C[k] = std::max(u1_C[k], upper_1[k][j]);
      }
    }
    u1_C[k] += delta_C[k];
    
    // update lower bound on distance to nearest member of each group 
    for (size_t g = 0; g < rf_n_groups; ++g){
      l_gC[k*rf_n_groups + g] -= max_delta_group[g];
      l_gC[k*rf_n_groups + g] -= delta_C[k];
    }
 
    min_cc = std::numeric_limits<double>::max();
    std::vector<size_t> groups_in_vicinity;

    // in own group, set distance to all centers other than self
    for (size_t kp = cum_in_group[rf_groups[k]]; kp < cum_in_group[rf_groups[k]+1]; ++kp){
      if (kp != k){
        set_rf_center_center_distance(k, kp, std::numeric_limits<double>::max(), cc_all[kp]);
        min_cc = std::min(min_cc, cc_all[kp]);
      }
    }
    groups_in_vicinity.push_back(rf_groups[k]);
    
    // we want to determine which centers are within the maximum exponion radius : 
    // 2*(max distance to an element of k) + (min distance to another center from k) (**)
    // we have (1) a LOWER bound on the distance on the centers in each group to center k, and
    //         (2) a UPPER bound on (**) of 2*u1_C[k] + (min_cc so far) + 
    
    for (size_t g = 0; g < rf_n_groups; ++g){
      if ((g != rf_groups[k]) && (l_gC[k*rf_n_groups + g] < 2*u1_C[k] + min_cc)){ // do not need to add s to make max_delta_possible_second correct.
        groups_in_vicinity.push_back(g);
        l_gC[k*rf_n_groups + g] = std::numeric_limits<double>::max();
        for (size_t kp = cum_in_group[g]; kp < cum_in_group[g+1]; ++kp){
          set_rf_center_center_distance(k, kp, std::numeric_limits<double>::max(), cc_all[kp]);
          min_cc = std::min(min_cc, cc_all[kp]);
          l_gC[k*rf_n_groups + g] = std::min(l_gC[k*rf_n_groups + g], cc_all[kp]);
        }
      }
    }
    
    n_cc_required = 0;  
    double max_delta_possible_second = 0;
    for (auto & g : groups_in_vicinity) {
      for (size_t kp = cum_in_group[g]; kp < cum_in_group[g+1]; ++kp) {
        if ((kp != k) && (cc_all[kp] < 2*u1_C[k] + min_cc)){ // do not need to add max_delta_C to make max_delta_possible_second correct.
          cc_required_tuple[n_cc_required] = std::make_tuple(cc_all[kp], kp);
          ++n_cc_required;
          max_delta_possible_second = std::max(max_delta_possible_second, delta_C[kp]);
        }
      }
    }
    
    cc_required_tuple[n_cc_required] = std::make_tuple(std::numeric_limits<double>::max(), std::numeric_limits<size_t>::max());
    ++n_cc_required;
    
    
    std::sort(cc_required_tuple.begin(), cc_required_tuple.begin() + n_cc_required, tupsorter);
    
    

    for (size_t jp = 0; jp < n_cc_required; ++jp){
      cc_required_distances[jp] = std::get<0>(cc_required_tuple[jp]);
      cc_required_indices[jp] = std::get<1>(cc_required_tuple[jp]);
    }



    // TODO : declare vectors outside of k-loop. 
    // all unique groups in cc_required_indices, in order of appearance in cc_required_indices.
    std::vector<size_t> cc_required_groups;
    cc_required_groups.push_back(rf_groups[cc_required_indices[0]]); 
    
    // how many unique groups have appeared at and preceding index jp in cc_required_indices. 
    std::vector<size_t> cc_n_groups_required(n_cc_required, std::numeric_limits<size_t>::max());
    cc_n_groups_required[0] = 1;

    for (size_t jp = 1; jp < n_cc_required; ++jp){
      auto new_g = rf_groups[cc_required_indices[jp]];
      cc_n_groups_required[jp] = cc_n_groups_required[jp - 1];
      if (std::count(cc_required_groups.begin(), cc_required_groups.end(), new_g) == 0){
        cc_required_groups.push_back(new_g);
        ++cc_n_groups_required[jp];
      }
    }
    
    
    double compl_exp_rad_x;
    double exprad;
    size_t kp;

    
    if (is_rf_correct_d1_round() == true){
      rf_tighten_nearest(k);
    }
    
    else{
      if (delta_C[k] > 0){
        for (size_t j = 0; j < get_ndata(k); ++j){
          upper_1[k][j] += delta_C[k];
        }
      }
    }

    for (size_t j = 0; j < get_ndata(k); ++j){


      size_t ID = sample_IDs[k][j]; //TODO, starting_IDs //DOES NOT HELP: rf_glt_ID[k][j];//
      lower_2[k][j] -= max_delta_possible_second; // this is better than max_delta_C. but is it correct? Yes, but TODO the name is wrong
      if (lower_2[k][j] < upper_1[k][j]){

        set_rf_center_sample_distance_nothreshold(k, k, j, min_d1); 
        reset_nearest_info(k, j, k, min_d1, f_energy(min_d1)); 
      
        if (lower_2[k][j] < upper_1[k][j]){
        
          set_rf_center_sample_distance_nothreshold(v_b[k][j], k, j, min_d2);
          compl_exp_rad_x = (min_d2 < min_d1 ? 2*min_d1 : min_d1 + min_d2);
          exprad = (1 + 1e-6)*std::min(2*min_d1 + min_cc, compl_exp_rad_x);
          min_k1 = k;
          min_k2 = v_b[k][j];
          
          if (min_d1 > min_d2){
            std::swap(min_d1, min_d2);
            std::swap(min_k1, min_k2);
          }
    
          //get insertion index:
          size_t insertion_index = 0;
          while (exprad > cc_required_distances[insertion_index]){
            ++insertion_index;
          }
          
          bool simplexpo = (true);//insertion_index == 0); // || cc_n_groups_required[insertion_index-1] <= 10);
          
          if (simplexpo == false){
            
            for (size_t groupi = 0; groupi < cc_n_groups_required[insertion_index-1]; ++groupi){
              size_t group = cc_required_groups[groupi];

              // if the group bound is less than min_d2, then nearest or second nearest might be in the group. 
              if (rf_l_groups[ID*rf_n_groups + group] - rf_cum_u_delta_gC[rf_n_groups*rf_t_groups[ID*rf_n_groups + group] + group] < min_d2){
                rf_t_groups[ID*rf_n_groups + group] = rf_round;
                rf_l_groups[ID*rf_n_groups + group] = std::numeric_limits<double>::max();
                
                for (size_t kp = cum_in_group[group]; kp < cum_in_group[group+1]; ++kp){              
                  set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);

                  if (kp != v_b[k][j] && kp != k){
                    if (worker_dist < min_d1){
                      min_d2 = min_d1;
                      min_k2 = min_k1;
                      min_d1 = worker_dist;
                      min_k1 = kp;
                    }
                    
                    else if (worker_dist < min_d2){
                      min_d2 = worker_dist;
                      min_k2 = kp;        
                    }
                  }
                
                  rf_l_groups[ID*rf_n_groups + group] = std::min(rf_l_groups[ID*rf_n_groups + group], worker_dist);
                }
              }
            }
          }
          
          
          else if (simplexpo == true){
            for (size_t kpi = 0; kpi < insertion_index; ++kpi){
              kp = cc_required_indices[kpi];
              if (kp != v_b[k][j]){
                
                set_rf_center_sample_distance(kp, k, j, min_d2, worker_dist);
                
                if (worker_dist < min_d1){
                  min_d2 = min_d1;
                  min_k2 = min_k1;
                  min_d1 = worker_dist;
                  min_k1 = kp;
                }
                
                else if (worker_dist < min_d2){
                  min_d2 = worker_dist;
                  min_k2 = kp;        
                }
              }
            }
          }
  
          reset_nearest_info(k,j, min_k1, min_d1, f_energy(min_d1));        
          lower_2[k][j] = min_d2;
          v_b[k][j] = min_k2;
          upper_1[k][j] = min_d1;
          
          if (min_d1 > u1_C[min_k1]){
            u1_C[min_k1] = min_d1;
          }
        }
      }
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

      rf_cumulative_correction(k_new, k, j); //correct sums or whatever. 
      
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
  upper_1[k_new].push_back(upper_1[k][j]);
  rf_glt_ID[k_new].push_back(rf_glt_ID[k][j]);
    
}

void SkeletonClusterer::rf_custom_replace_with_last(size_t k, size_t j){  
  lower_2[k][j] = lower_2[k].back();
  v_b[k][j] = v_b[k].back();
  upper_1[k][j] = upper_1[k].back();
  rf_glt_ID[k][j] = rf_glt_ID[k].back();

}

void SkeletonClusterer::rf_custom_remove_last(size_t k){
  lower_2[k].pop_back();
  v_b[k].pop_back();
  upper_1[k].pop_back();
  rf_glt_ID[k].pop_back();  
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
void SkeletonClusterer::rf_update_energies(){
  old_E_total = E_total;
  E_total = 0;
  for (size_t k = 0; k < K; ++k){
    cluster_energies[k] = 0;
    for (size_t j = 0; j < get_ndata(k); ++j){
      cluster_energies[k] += nearest_1_infos[k][j].e_x;
    }
    cluster_mean_energies[k] = cluster_energies[k] / static_cast<double> (get_ndata(k));
    E_total += cluster_energies[k];
  }
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

//Standard algorithm. 
void SkeletonClusterer::rf_update_sample_info_standard(){

  size_t min_k1, min_k2;
  double min_d1, min_d2;
  
  for (size_t k = 0; k < K; ++k){
    for (size_t j = 0; j < get_ndata(k); ++j){
      rf_set_2_smallest(k, j, min_k1, min_d1, min_k2, min_d2);
      reset_nearest_info(k,j, min_k1, min_d1, f_energy(min_d1));
    }
  }
}

void SkeletonClusterer::set_rf_center_sample_distance_nothreshold(size_t k, size_t k1, size_t j1, double & distance){
  set_rf_center_sample_distance(k, k1, j1, std::numeric_limits<double>::max(), distance);
}

    
}
