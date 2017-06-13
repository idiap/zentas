#ifndef ZENTAS_SKELETONCLUSTERER_HPP
#define ZENTAS_SKELETONCLUSTERER_HPP


#include "outputwriter.hpp"
#include "zentaserror.hpp"
#include "tenergyfunc.hpp"
#include "energyinit.hpp"
#include "initialisation.hpp"

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
    /* also tried with d1 and d2 separate, slower (as always used together). */
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

/* The choice to store (a,d,e) contiguously for a sample was more for ease of coding than performance. 
 * However, I have implemented a version with a,d and e stored as separate vectors, with no speed-up */
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


class SkeletonClusterer{

  public:  
    
    SkeletonClusterer(size_t K, std::chrono::time_point<std::chrono::high_resolution_clock> bigbang, const size_t * const center_indices_init_, size_t ndata_in, std::string init_method, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t max_rounds, std::string energy, bool with_tests, size_t random_sd, const EnergyInitialiser & energy_initialiser);
    


  protected: //certain should be private after BaseClusterer, maybe even here.

    /* generic */
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

    bool in_refinement {false};
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




    size_t get_time_in_update_centers();    
    size_t get_time_in_update_sample_info();    
    double get_time_remaining();

    size_t get_start(size_t ti, size_t nthreads, size_t j_A, size_t j_Z);
    size_t get_end(size_t ti, size_t nthreads, size_t j_A, size_t j_Z);
    size_t get_sample_from(std::vector<double> & v_cum_nearest_energies);
    std::string get_base_summary_string();    
    
    /* metric virtuals */
    virtual double get_rel_calccosts() = 0;
    virtual size_t get_ncalcs() = 0;
    virtual void set_cc_pp_distance(size_t k1, size_t k2) = 0;
    virtual void set_center_sample_pp_distance(size_t k, size_t bin, size_t i, double & adistance) = 0;
    virtual void set_center_sample_distance(size_t k, size_t k1, size_t j1, double threshold, double & distance) = 0;
    virtual void set_center_sampleID_distance(size_t k, size_t i, double threshold, double & distance) = 0;
    virtual void set_sampleID_sampleID_distance(size_t i1, size_t i2, double threshold, double & distance) = 0;
    virtual void set_sample_sample_distance(size_t k1, size_t j1, size_t k2, size_t j2, double threshold, double & adistance) = 0;
    virtual void set_center_center_distance(size_t k1, size_t k2, double threshold, double & adistance) = 0;
    virtual void set_center_sample_distance(size_t k, size_t k1, size_t j1, double & distance) = 0;

    /* data virtuals */
    virtual std::string string_for_sample(size_t k, size_t j) = 0;
    virtual std::string string_for_center(size_t k) = 0;
    virtual void append_from_ID(size_t k, size_t i) = 0;
    virtual void append_pp_from_ID(size_t i) = 0;
    virtual void append_pp_from_bin(size_t bin, size_t j) = 0;
    virtual void append_aq2p_p2buns(size_t bin, size_t i) = 0;
    virtual size_t get_ndata(size_t k) = 0;
    virtual void swap_center_data(size_t k, size_t j) = 0;
    
    /* rule : functions with suffix 'basic' will not touch to_leave_cluster */
    void reset_nearest_info_basic(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest);
    void reset_nearest_info(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest);



};

}    
    

#endif
