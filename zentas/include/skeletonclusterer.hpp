#ifndef ZENTAS_SKELETONCLUSTERER_HPP
#define ZENTAS_SKELETONCLUSTERER_HPP


#include "outputwriter.hpp"
#include "zentaserror.hpp"
#include "tenergyfunc.hpp"
#include "energyinit.hpp"
#include "initialisation.hpp"
#include "zentasinfo.hpp"

#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <random>
#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>


namespace nszen{


/* k-means++ helper class */
class P2Bundle{

  private:
  
  size_t ndata;
  
  /* Note on data layout here:
   * I've also tried with d1 and d2 separate, 
   * slower (as d1 and d2 always used together). */
  std::unique_ptr <std::array<double, 2 > [] > d12_;    
  std::unique_ptr <size_t []> k_1_;
  std::unique_ptr <size_t []> k_2_;
  std::unique_ptr <size_t []> ori_;
  
  public:
  
  size_t & k_1(size_t i){
    return k_1_[i];
  }
  
  size_t & k_2(size_t i){
    return k_2_[i];
  }
  
  double & d_1(size_t i ){
    return std::get<0>(d12_[i]);
  }
  
  double & d_2(size_t i){
    return std::get<1>(d12_[i]);
  }
  
  size_t & ori(size_t i){
    return ori_[i];
  }
  
  size_t get_ndata() const{
    return ndata;
  }
  
  void initialise_data(size_t n);
  
  P2Bundle() = default;  
  P2Bundle(size_t n);
  P2Bundle(const std::vector<size_t> & ori);
  
};

/* Note on data layout here:
 * The choice to store (a,d,e) contiguously initially 
 * for clearer code. However, I have implemented a version 
 * with a,d and e stored as separate vectors, and there 
 * is no speed-up */
struct XNearestInfo{
  // the assigned center of the data point
  size_t a_x;
  // the distance to the assigned center
  double d_x;
  // the energy of d_x
  double e_x;
  
  void reset(size_t new_a_x, double new_d_x, double new_e_x);
  void reset(XNearestInfo & nearest_x_infos);
  XNearestInfo(size_t a_x, double d_x, double e_x);
  std::string get_string();
  
};


class SkeletonClustererInitBundle{
  public:
    size_t K;
    size_t ndata;
    std::chrono::time_point<std::chrono::high_resolution_clock> bigbang;
    const size_t * const center_indices_init_predefined;
    std::string init_method;
    double max_time;
    double min_mE;
    size_t max_rounds;
    size_t nthreads;
    size_t seed;
    std::string energy;
    bool with_tests;
    size_t * const indices_final;
    size_t * const labels;
    const EnergyInitialiser * ptr_energy_initialiser;
    
  SkeletonClustererInitBundle(size_t K_, size_t nd_, std::chrono::time_point<std::chrono::high_resolution_clock> bb_, const size_t * const center_indices_init_predefined_, const std::string & init_method_, double max_time_, double min_mE_, size_t max_rounds_, size_t nthreads_, size_t seed_, const std::string & energy_,  bool with_tests_, size_t * const indices_final_, size_t * const labels_, const EnergyInitialiser * ptr_energy_initialiser_);
  
  
};


class SkeletonClusterer{

  public:  

  
  SkeletonClusterer(const SkeletonClustererInitBundle & sb);


  /* TODO certain variables should be private. */
  protected: 
  
  std::vector<std::vector<size_t>> aq2p_original_indices;
  std::vector<P2Bundle> aq2p_p2buns;
  zentas::outputwriting::OutputWriter mowri;    
  const size_t K;    
  std::chrono::time_point<std::chrono::high_resolution_clock> bigbang;
  const size_t ndata;
  std::function<double(double)> f_energy;    
  std::string initialisation_method;
  size_t * center_IDs;
  std::vector<std::vector<XNearestInfo>> nearest_1_infos;
  std::vector<std::vector<size_t>> sample_IDs;
  std::vector<std::vector<size_t>> to_leave_cluster;
  std::vector<bool> cluster_has_changed;
  std::vector<double> cluster_energies;
  std::vector<double> cluster_mean_energies;
  double E_total;
  double old_E_total;
  size_t round;
  std::vector<size_t> v_center_indices_init;
  size_t * const center_indices_init;      
  size_t time_prehistory = 0;
  size_t time_in_update_centers = 0;
  size_t time_in_update_sample_info = 0;
  size_t time_in_redistribute = 0;
  size_t time_initialising = 0;
  size_t time_in_update_all_cluster_statistics = 0;
  size_t time_total = 0;
  size_t time_to_initialise_centers = 0;
  std::chrono::time_point<std::chrono::high_resolution_clock> tstart;
  std::chrono::time_point<std::chrono::high_resolution_clock> tstart_initialise_centers;
  size_t ncalcs_in_update_centers = 0;
  size_t ncalcs_in_update_sample_info = 0;
  size_t ncalcs_initialising = 0;
  size_t ncalcs_total = 0;
  /* passed in through constructor */
  std::chrono::time_point<std::chrono::high_resolution_clock> 
  t_update_centers_start, 
  t_update_centers_end, 
  t_update_sample_info_end, 
  t_redistribute_end, 
  t_update_all_cluster_statistics_end;
  size_t 
  ncalcs_update_centers_start, 
  ncalcs_update_centers_end, 
  ncalcs_update_sample_info_end;

  size_t max_time_micros;
  double min_mE;
  size_t * const labels;
  std::mutex mutex0;
  size_t nthreads;
  size_t nthreads_fl;
  size_t max_rounds;
  std::string energy;
  bool with_tests;
  std::uniform_int_distribution<size_t> dis;
  /* generates [0,1] uniform, will use to sample centers */
  std::uniform_real_distribution<double> dis_uni01;          
  std::default_random_engine gen;
  std::vector<double> kmoo_cc;
  P2Bundle kmoo_p2bun;
  bool km_is_exhausted;
 
   /* *****************
   * metric virtuals *
   * ***************** */
  virtual double get_rel_calccosts() = 0;
  virtual size_t get_ncalcs() = 0;
  virtual void set_cc_pp_distance(size_t k1, size_t k2) = 0;
  virtual void set_center_sample_pp_distance(size_t k, size_t bin, size_t i, double & adistance) = 0;
  virtual void set_center_sample_distance(size_t k, size_t k1, size_t j1, double threshold, double & distance) = 0;
  virtual void set_center_sampleID_distance(size_t k, size_t i, double threshold, double & distance) = 0;
  virtual void set_sampleID_sampleID_distance(size_t i1, size_t i2, double threshold, double & distance) = 0;
  virtual void set_sample_sample_distance(size_t k1, size_t j1, size_t k2, size_t j2, double threshold, double & adistance) = 0;
  virtual void set_center_center_distance(size_t k1, size_t k2, double threshold, double & adistance) = 0;
  virtual bool get_do_refinement() = 0;

  /* ***************
   * data virtuals *
   * *************** */
  virtual std::string string_for_sample(size_t k, size_t j) = 0;
  virtual std::string string_from_ID(size_t i) = 0;
  virtual std::string string_for_center(size_t k) = 0;
  virtual void append_from_ID(size_t k, size_t i) = 0;
  virtual void append_to_centers_from_ID(size_t i) = 0;
  virtual void add_empty_cluster() = 0;
  virtual void append_pp_from_ID(size_t i) = 0;
  virtual void append_pp_from_bin(size_t bin, size_t j) = 0;
  virtual void append_aq2p_p2buns(size_t bin, size_t i) = 0;
  virtual size_t get_ndata(size_t k) = 0;
  virtual void swap_center_data(size_t k, size_t j) = 0;
  virtual void replace_with_last_element(size_t k, size_t j) = 0;
  virtual void remove_last(size_t k) = 0;
  virtual void  print_centers() = 0;
  virtual void append_across(size_t k_to, size_t k_from, size_t j_from) = 0;
  virtual void replace_with(size_t k_to, size_t j_to, size_t k_from, size_t j_from) = 0;
  virtual size_t get_p2buns_dt_ndata(unsigned bini) = 0;
  virtual void reset_p2buns_dt(unsigned n_bins) = 0;
  virtual void kmoo_finish_with() = 0;
  virtual void centers_replace_with_sample(size_t k1, size_t k2, size_t j2) = 0;
  virtual void swap_data(size_t k1, size_t k2) = 0;


  /* general rule : functions with suffix 
   * 'basic' will not touch to_leave_cluster */
  void reset_nearest_info_basic(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest);
  void reset_nearest_info(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest);
  size_t get_time_in_update_centers();    
  size_t get_time_in_update_sample_info();    
  double get_time_remaining();
  size_t get_start(size_t ti, size_t nthreads, size_t j_A, size_t j_Z);
  size_t get_end(size_t ti, size_t nthreads, size_t j_A, size_t j_Z);
  size_t get_sample_from(std::vector<double> & v_cum_nearest_energies);
  std::string get_base_summary_string();    


  public:
  
  virtual void set_center_center_info() = 0;
  virtual void update_center_center_info() = 0;
  void signal_cluster_change(size_t k);
  size_t draw_j_uniform(size_t k);
  void kmoo_prepare();
  void print_ndatas();
  void default_initialise_with_kmeanspp();
  void triangular_kmeanspp_aq2(size_t n_bins);
  void triangular_kmeanspp();
  void kmpp_inner(size_t i, size_t k, double a_distance, P2Bundle & p2bun);
  size_t get_nthreads();
  double get_nthreads_fl();
  void triangular_kmeanspp_after_initial(size_t aq2p_bin, size_t k0, size_t k1, bool from_full_data);
  void test_parameters_to_tkai(size_t k0, size_t k1, size_t ndata_1, size_t ndata_2);
  void set_center_sample_distance_nothreshold(size_t k, size_t k1, size_t j1, double & distance);
  void set_rf_center_sample_distance_nothreshold(size_t k, size_t k1, size_t j1, double & distance);
  void set_sampleID_sampleID_distance_nothreshold(size_t i1, size_t i2, double & distance);
  double get_center_sample_distance_nothreshold(size_t k, size_t k1, size_t j1);
  double get_sample_sample_distance_nothreshold(size_t k , size_t j1, size_t j2);
  void set_center_center_distance_nothreshold(size_t k1, size_t k2, double & adistance);
  void output_halt_kmedoids_reason();
  bool halt_kmedoids();

  void set_t_ncalcs_update_centers_start();
  void update_t_ncalcs_center_end();
  void update_t_ncalcs_sample_update_end();
  void update_t_redistibute_end();
  void update_t_update_all_cluster_stats_end();
  
  void run_kmedoids();    
  
  void core_kmedoids_loops();  
  void populate_labels();
  void go();
      

  /* TODO certain variables should be private. */
  protected:
  
  virtual void put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) = 0;
  virtual void initialise_with_kmeanspp() = 0;
  virtual void put_sample_custom_in_cluster(size_t, size_t, const double * const) {}
  /* when implementing this, one must guarantee that cluster_has_changed is correctly set in this function */
  virtual void reset_sample_custom(size_t, size_t, size_t, const double * const) {}
  virtual bool update_centers() = 0;       
  virtual void update_sample_info() = 0;
  virtual void custom_append(size_t, size_t, size_t) {}
  virtual void custom_replace_with_last(size_t, size_t) {}
  virtual void custom_replace_with(size_t, size_t, size_t, size_t) {};
  virtual void custom_remove_last(size_t) {};
  virtual void increment_custom_cluster_statistics(size_t, size_t) {};
  virtual void set_normalised_custom_cluster_statistics(size_t k) = 0;
  virtual void set_to_zero_custom_cluster_statistics(size_t k) = 0;
  virtual std::string get_round_summary() = 0;    

  private:
  ////////////////////////////
  ////////// tests ///////////
  ////////////////////////////
  void post_initialisation_test();
  virtual void center_center_info_test();
  void post_center_update_test();
  void post_sample_update_test();
  void post_redistribute_test();    
  void post_update_statistics_test();
  void to_leave_cluster_test();
  virtual void custom_cluster_statistics_test();
  void cluster_statistics_test();    
  void fundamental_triangle_inequality_test();
  virtual void custom_info_test();
  void info_tests();
  virtual void custom_ndata_test();
  void ndata_tests();    
  void as_assigned_test();
  void injective_ID_test();
  void post_initialise_centers_test();



  protected:

  /* single threaded redistribution. mulithreading here 
   * will be very tricky. requires that to_leave_cluster
   * be reliably set. centers are unchanged by this
   * function. this function simply moves samples between 
   * clusters. additional complexity required to maintain 
   * randomness (immigrant samples inserted into random indices)
   * */
  void redistribute();
  void set_all_cluster_statistics();
  void update_all_cluster_statistics();
  virtual void put_sample_in_cluster(size_t i) = 0;
  
  /* redistribute_order is the order in which clusters `send away samples' */
  /* this function : set the vector to be indices {0,...K-1} in a specific order */
  /* clarans is best when k_to is first (done in virtual of baseclarans) */
  virtual void set_redistribute_order(std::vector<size_t> & redistribute_order) = 0;
  
  /* set energy of cluster k, and set as custom cluster statistics */
  void set_cluster_statistics(size_t k);

  size_t get_a1(size_t k, size_t j){
    return nearest_1_infos[k][j].a_x;
  }
  
  double get_d1(size_t k, size_t j){
    return nearest_1_infos[k][j].d_x;
  }
  
  double get_e1(size_t k, size_t j){
    return nearest_1_infos[k][j].e_x;
  }
  
  double get_e1_tail(size_t k){
    return nearest_1_infos[k].back().e_x;
  }

  double get_E_total(){
    return E_total;
  }
  
  double get_cluster_energy(size_t k){
    return cluster_energies[k];
  }

  /* remove the j'th sample from cluster k, and 
   * (if it is not the last element) fill the 
   * hole with the tail (size drops by 1) */
  void remove_with_tail_pull(size_t k, size_t j);
  void initialise_all() ;    
  void put_samples_in_clusters();
  void pll_put_samples_in_cluster(size_t i_a, size_t i_z, std::vector<size_t> & non_center_IDs);
  void final_push_into_cluster_post_kmeanspp(size_t i, size_t k1, size_t k2, double d1, double d2);


  public:
  
  void final_push_into_cluster_basic(size_t i, size_t nearest_center, double min_distance);
  void reset_multiple_sample_infos(size_t k_to, size_t j_a, size_t j_z);
  void base_put_sample_in_cluster(size_t i) ;    
  void triangular_put_sample_in_cluster(size_t i, const double * const cc);
  void swap_center_with_sample(size_t k, size_t j);
  void move_center_into_its_own_cluster(size_t k);    
  size_t draw_k_uniform();    
  size_t draw_k_prop_ndata();
  void reset_sample_infos_basic(size_t k, size_t j);    
  void reset_sample_infos(size_t k, size_t j);
  void final_push_into_cluster(size_t i, size_t nearest_center, double min_distance, const double * const distances); 
  void overwrite_center_with_sample(size_t k1, size_t k2, size_t j2);

  void populate_afk_mc2();
  void initialise_center_indices();
  

  public:

  void initialise_refinement();
  void run_refinement();






  /* refinement (k-means etc.) functions defined in refinement.cpp */

  std::vector<std::vector<double>> lower_2;
  std::vector<std::vector<size_t>> v_b; // second nearest, when set.
  std::vector<std::vector<double>> upper_1;
  std::vector<double> delta_C;
  std::vector<double> u1_C; // upper bound on furthest cluster member.   
  size_t rf_n_groups;
  std::vector<double> l_gC; //of size K * rf_n_groups.
  std::vector<size_t> n_in_group;
  std::vector<size_t> cum_in_group;
  std::vector<double> max_delta_group; //of size rf_n_groups.
  std::vector<size_t> rf_groups;
  std::vector<std::vector<size_t>> rf_glt_ID; 
  std::vector<double> rf_l_groups; // lower bounds on distances to groups. size ndata*rf_n_groups. mapped from rf_glt_ID. 
  std::vector<size_t> rf_t_groups; // times when lower bounds set. size ndata*rf_n_groups.  mapped from rf_glt_ID.  
  std::vector<double> rf_cum_u_delta_gC; // upper bounds on distances moved since start, by group (of size T*rf_n_groups) 
  size_t rf_round;
  
  bool halt_refinement();
  bool rf_update_centers();
  bool is_rf_correct_d1_round();
  void rf_set_2_smallest(size_t k, size_t j, size_t & min_k1, double & min_d1, size_t & min_k2, double & min_d2);
  void rf_update_center_center_info();
  void rf_update_sample_info();
  void rf_update_sample_info_standard();  
  void rf_update_sample_info_hamerly();
  void rf_update_sample_info_exponion();
  void rf_swap(size_t k1, size_t k2);
  virtual void custom_initialise_refinement() {}
  virtual void custom_rf_clear_initmem() {}
  std::string rf_get_round_summary();
  void rf_redistribute();
  void rf_update_energies();
  void rf_tighten_nearest(size_t k);
  bool is_rf_tighten_cluster_radius_round();
  void custom_swap(size_t k1, size_t k2);
  void rf_custom_append(size_t k_new, size_t k, size_t j);
  void rf_custom_replace_with_last(size_t k, size_t j);
  void rf_custom_remove_last(size_t k);
  void rf_post_center_update_test();
  void rf_post_sample_update_test();
  void rf_post_redistribute_test();
  void rf_post_update_statistics_test();




  public:
  
  /* data */
  virtual void zero_refinement_sum(size_t k) { (void)k; throw zentas::zentas_error("zero_refinement_sum not possible"); };
  virtual void add_to_refinement_sum(size_t k, size_t j) {(void)k; (void)j; throw zentas::zentas_error("add_to_refinement_sum not possible"); };
  virtual void subtract_from_refinement_sum(size_t k, size_t j) {(void)k; (void)j;  throw zentas::zentas_error("subtract_from_refinement_sum not possible"); };
  virtual void append_zero_to_rf_center_data(){throw zentas::zentas_error("virtual function append_zero_to_rf_center_data not possible");  }
  virtual void append_zero_to_rf_sum_data(){throw zentas::zentas_error("virtual function append_zero_to_rf_sum_data not possible");  }
  //virtual void a_to_rf_sum_data(){throw zentas::zentas_error("virtual function append_zero_to_rf_sum_data not possible");  }
  virtual void append_zero_to_old_rf_center_data(){throw zentas::zentas_error("virtual function append_zero_to_rf_center_data not possible");  }
  virtual void set_old_rf_center_data(size_t k){(void)k; throw zentas::zentas_error("virtual function append_zero_to_rf_center_data not possible");  }
  virtual void set_delta_rf_new_and_old(size_t k, double threshold, double &) {(void)k; (void)threshold; throw zentas::zentas_error("virtual function get_delta_rf_new_and_old not possible");}
  virtual void set_sum_abs_rf_new_and_old(size_t k, double &) {(void)k; throw zentas::zentas_error("virtual function get_sum_abs_rf_new_and_old not possible");}
  virtual void set_rf_center_data(size_t k){(void)k; throw zentas::zentas_error("virtual function set_rf_center_data not possible"); }
  virtual std::string string_for_rf_center(size_t k){(void)k; throw zentas::zentas_error("strinf_for_rf_center not possible");}

  /* metric */
  virtual void set_rf_center_sample_distance(size_t k, size_t k1, size_t j1, double threshold, double & distance)
  {(void)k; (void)k1; (void)j1; (void)threshold; (void)distance; throw zentas::zentas_error("virtual function set_rf_center_sample_distance not possible"); }
  virtual void set_rf_center_center_distance(size_t k1, size_t k2, double threshold, double & adistance)
  {(void)k1; (void)k2; (void)threshold; (void)adistance; throw zentas::zentas_error("virtual function set_rf_center_center_distance not possible"); }
  virtual std::vector<size_t> get_subclustered_centers_labels(size_t sub_K) {(void)sub_K; throw zentas::zentas_error("cluster_centers not possible"); };
  virtual void rf_cumulative_correction(size_t k_new, size_t k, size_t j) {(void)k_new, (void)k; (void)j; throw zentas::zentas_error("rf_cumulative_correction not possible"); };
  virtual void rf_increment_sum(size_t k, size_t j) {(void)k; (void)j; throw zentas::zentas_error("rf_increment_sum not possible"); };

};

}    
    

#endif
