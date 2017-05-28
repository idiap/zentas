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
#ifndef ZENTAS_BASECLUSTERER_HPP
#define ZENTAS_BASECLUSTERER_HPP

#include "tdata.hpp"
#include "tdatain.hpp"
#include "tmetric.hpp"
#include "tenergyfunc.hpp"

#include "outputwriter.hpp"
#include "initialisation.hpp"

#include <algorithm>
#include <type_traits>
#include <chrono>
#include <random>


#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>



namespace nszen{

class P2Bundle{

  public:

    size_t ndata;
    
    
    /* also tried with d1 and d2 separate, slower (as always used together). */
    std::unique_ptr <std::array<double, 2 > [] > d12_;    
    std::unique_ptr <size_t []> k_1_;
    std::unique_ptr <size_t []> k_2_;
    std::unique_ptr <size_t []> ori_;
    
    
    inline size_t & k_1(size_t i){
      return k_1_[i];
    }
    
    inline size_t & k_2(size_t i){
      return k_2_[i];
    }
    
    inline double & d_1(size_t i ){
      return std::get<0>(d12_[i]);
    }

    inline double & d_2(size_t i){
      return std::get<1>(d12_[i]);
    }
    
    inline size_t & ori(size_t i){
      return ori_[i];
    }
    
    inline size_t get_ndata() const{
      return ndata;
    }
    
    void initialise_data(size_t n){
      ndata = n;
      d12_.reset(new std::array<double, 2> [n]);
      k_1_.reset(new size_t[n]);
      k_2_.reset(new size_t[n]);
      ori_.reset(new size_t[n]);
      
      for (unsigned i = 0; i < n; ++i){
        d_1(i) = std::numeric_limits<double>::max();
        d_2(i) = std::numeric_limits<double>::max();
      }
      for (unsigned i = 0; i < n; ++i){
        k_1(i) = 0; /* note that in the first round of kmeanspp, the lookup in cc requires this  */ 
      }
    }
    
    P2Bundle() = default;
    
    P2Bundle(size_t n){
      initialise_data(n);
    }
    
    /* TODO : how to use preceding initialised to initialise here */
    P2Bundle(const std::vector<size_t> & ori){
      initialise_data(ori.size());
      for (size_t i = 0; i < ndata; ++i){
        ori_[i] = ori[i];
      }
    }

};



template <typename TData>
class KmooBundle {
  
  public:
  
  std::vector<double> cc;
  P2Bundle p2bun;
  TData c_dt;
  
  typedef typename TData::DataIn DataIn;

  KmooBundle(size_t K, size_t ndata, const DataIn * ptr_datain):cc(K*K), p2bun(ndata), c_dt(*ptr_datain, true){
    
  }
  
  KmooBundle(const DataIn * ptr_datain):c_dt(*ptr_datain, true) {}
  
  void reset(size_t K, size_t ndata){
    cc.resize(K*K);
    p2bun = P2Bundle(ndata);
  }
  
};



 
 
struct EnergyInitialiser{

  private:
    double critical_radius;
    double exponent_coeff;
    
  public: 
    EnergyInitialiser():critical_radius(0), exponent_coeff(0) {}
    EnergyInitialiser(double critical_radius, double exponent_coeff):critical_radius(critical_radius), exponent_coeff(exponent_coeff) {}
    
    double get_critical_radius () const{
      return critical_radius;
    }
    
    double get_exponent_coeff() const{
      return exponent_coeff;
    }
   
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
  
  inline void reset(size_t new_a_x, double new_d_x, double new_e_x){
    a_x = new_a_x;
    d_x = new_d_x;
    e_x = new_e_x;
  }
  
  inline void reset(XNearestInfo & nearest_x_infos){
    a_x = nearest_x_infos.a_x;
    d_x = nearest_x_infos.d_x;
    e_x = nearest_x_infos.e_x;
  }

  
  XNearestInfo(size_t a_x, double d_x, double e_x):a_x(a_x), d_x(d_x), e_x(e_x){}
  
  std::string get_string(){
    return  std::to_string(a_x) + "\t " + std::to_string(d_x) + "\t " + std::to_string(e_x);
  }
};


template <class TDataIn, class TMetric>
struct BaseClustererInitBundle{

  typedef typename TMetric::Initializer TMetricInitializer;

  size_t K;
  const TDataIn * const ptr_datain;
  const size_t * const center_indices_init;
  std::string initialisation_method;
  size_t seed;
  double maxtime;
  double minmE;
  size_t * const indices_final;
  size_t * const labels;
  size_t nthreads;
  size_t maxrounds;
  std::string energy;
  const TMetricInitializer * ptr_metric_initializer;
  const EnergyInitialiser * ptr_energy_initialiser;

  
  BaseClustererInitBundle(size_t K, const TDataIn * const ptr_datain, const size_t * const center_indices_init, std::string initialisation_method_, size_t seed, double maxtime, double minmE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t maxrounds, std::string energy, const TMetricInitializer & metric_initializer, const EnergyInitialiser & energy_initialiser): 
  
  K(K), ptr_datain(ptr_datain), center_indices_init(center_indices_init), initialisation_method(initialisation_method_), seed(seed), maxtime(maxtime), minmE(minmE), indices_final(indices_final), labels(labels), nthreads(nthreads), maxrounds(maxrounds), energy(energy), ptr_metric_initializer(&metric_initializer), ptr_energy_initialiser(&energy_initialiser) {}

};

template <class TMetric, class TData>
class BaseClusterer{
  
  /* A comment to put somewhere : when two clusters are equally near to a sample, the cluster to which it is assigned is unpredictable. One might have a statement that it goes to the one with lower index for example, but this would mean that if a center is at equal distance to the nearest it cannot be eliminated. This is less efficient, and well, I don't want to go back and change vast swathes of code. 
   * 
   * As a result of this unpredictability, l0, l1, l2 can differ: the important point is that random proposals depend on cluster sizes, and so while the energy and proposal evaluation do not care whether a sample is assigned to a cluster at equal distance to another, the subsequent proposals may differ. 
   * 
   * */
  
  public:
  
    zentas::outputwriting::OutputWriter mowri;

    typedef typename TMetric::Initializer TMetricInitializer;
    typedef typename TData::DataIn DataIn;
    const size_t K;
    const size_t ndata;
    std::function<double(double)> f_energy;
    std::string initialisation_method;

  private:
    TData centers_data;
    size_t * center_IDs;
    std::vector<TData> cluster_datas; 
    std::vector<std::vector<XNearestInfo>> nearest_1_infos;
    std::vector<std::vector<size_t>> sample_IDs;
    std::vector<std::vector<size_t>> to_leave_cluster;
    std::vector<bool> cluster_has_changed;
    const DataIn * const ptr_datain;
    TMetric metric;
    
    std::vector<double> cluster_energies;
    std::vector<double> cluster_mean_energies;
    double E_total;
    double old_E_total;
    size_t round;
    
    std::vector<size_t> v_center_indices_init;
    size_t * const center_indices_init;    
    
    
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
    


  protected:
    std::uniform_int_distribution<size_t> dis;

    /* generates [0,1] uniform, will use to sample centers */
    std::uniform_real_distribution<double> dis_uni01;      
    
    std::default_random_engine gen;
  
  private:  
    size_t maxtime_micros;
    double minmE;
    
    size_t * const labels;
    
    std::mutex mutex0;
    
    size_t nthreads;
    size_t nthreads_fl;

    size_t maxrounds;
    
    std::string energy;

    std::unique_ptr<KmooBundle<TData>> up_kmoo_bundle;
    
    

/////////////////////////////////////
/////////// public //////////////////
/////////////////////////////////////
  public:
    BaseClusterer(size_t K, /* number of clusters */
     const DataIn & datain, /* initialisating data */
     const size_t * const center_indices_init_, /* The K sample indices to initialise with (if initialisation_method is from_indices_init)*/
     std::string initialisation_method, /* */
     size_t seed, /* starting seed for generating random numbers, uses a custom "std::default_random_engine"  */
     double maxtime, /* will stop in go ( ) at first opportunity after maxtime */
     double minmE, /* will stop in go ( ) at first opportunity if mE is less than minmE */
     size_t * const indices_final, /* the K final indices (to populate) */
     size_t * const labels, /* the assigned cluster of the datain.ndata samples */
     size_t nthreads,
     size_t maxrounds, 
     std::string energy, 
     const TMetricInitializer & metric_initializer, 
     const EnergyInitialiser & energy_initialiser): 
    
    mowri(true, false, ""),
    K(K), ndata(datain.get_ndata()), initialisation_method(initialisation_method),
    centers_data(datain, true), nearest_1_infos(K), sample_IDs(K), to_leave_cluster(K), cluster_has_changed(K, true), ptr_datain(& datain), 
    metric(datain, nthreads, metric_initializer), 
    cluster_energies(K,0), cluster_mean_energies(K), E_total(std::numeric_limits<double>::max()), old_E_total(0),
    round(0), v_center_indices_init(K), center_indices_init(v_center_indices_init.data()), gen(seed), maxtime_micros(static_cast<size_t>(maxtime*1000000.)), minmE(minmE), labels(labels), nthreads(nthreads), nthreads_fl(static_cast<double> (nthreads)), maxrounds(maxrounds), energy(energy), up_kmoo_bundle(new KmooBundle<TData>(ptr_datain))

     {
       
      
      

       
       
       
      if (energy.compare("identity") == 0){
        f_energy = nszen::Identity(); 
      } 
       
      else if (energy.compare("quadratic") == 0){
        f_energy = nszen::Quadratic();
      }
      
      else if (energy.compare("cubic") == 0){
        f_energy = nszen::Cubic();
      }
      
      else if (energy.compare("squarepotential") == 0){
        f_energy = nszen::SquarePotential(energy_initialiser.get_critical_radius());
      }
      
      else if (energy.compare("log") == 0){
        f_energy = nszen::Log();
      }
      
      else if (energy.compare("exp") == 0){
        f_energy = nszen::Exponential(energy_initialiser.get_exponent_coeff());
      }
      
      else if (energy.compare("sqrt") == 0){
        f_energy = nszen::SquareRoot();
      }
      
      else{
        throw std::runtime_error(std::string("Unrecognised energy function, ") + energy);
      }
      


      /* initialisation from indices. */
      if (initialisation_method == "from_indices_init"){
        populate_from_indices_init(center_indices_init_, center_indices_init, K, ndata);
      }
      
      else{
        //will do in go ( )
      }

      center_IDs = indices_final;

    }







  
    BaseClusterer(const BaseClustererInitBundle<DataIn, TMetric> & ib): BaseClusterer(ib.K, *(ib.ptr_datain), ib.center_indices_init, ib.initialisation_method, ib.seed, ib.maxtime, ib.minmE, ib.indices_final, ib.labels, ib.nthreads, ib.maxrounds, ib.energy, *(ib.ptr_metric_initializer), *(ib.ptr_energy_initialiser)) {}

    size_t get_time_in_update_centers(){
      return time_in_update_centers;
    }
    
    size_t get_time_in_update_sample_info(){
      return time_in_update_sample_info;
    }
    
    inline double get_time_remaining(){
      auto t1 = std::chrono::high_resolution_clock::now();
      time_total = std::chrono::duration_cast<std::chrono::microseconds>(t1 - tstart).count();
        
      
      if (maxtime_micros > time_total){
        return maxtime_micros - time_total;
      }
      else{
        return -1;
      }
    }


    inline void final_push_into_cluster_basic(size_t i, size_t nearest_center, double min_distance){
      cluster_datas[nearest_center].append(ptr_datain->at_for_move(i));
      nearest_1_infos[nearest_center].emplace_back(nearest_center, min_distance, f_energy(min_distance));
      sample_IDs[nearest_center].push_back(i);      
    }

    //TODO : move to private section
    inline void final_push_into_cluster(size_t i, size_t nearest_center, double min_distance, const double * const distances){
      //get your lock on, time for polyphonics. 
      std::lock_guard<std::mutex> lockraii(mutex0);
      
      /* the common part of pushing into a cluster */
      final_push_into_cluster_basic(i, nearest_center, min_distance);
      
      /* as usual: the final parameter up_distances.get() guarantees correctness only of lowest and second lowest, other values may exceed */
      //TODO : there may be computation in here which should not be under the lock. 
      put_sample_custom_in_cluster(i, nearest_center, distances);

    }
    
    inline void final_push_into_cluster_post_kmeanspp(size_t i, size_t k1, size_t k2, double d1, double d2){
      std::lock_guard<std::mutex> lockraii(mutex0);
      final_push_into_cluster_basic(i, k1, d1);
      put_nearest_2_infos_margin_in_cluster_post_kmeanspp(k1, k2, d2, f_energy(d2));
    }
    


    void reset_multiple_sample_infos(size_t k_to, size_t j_a, size_t j_z){
      for (size_t j = j_a; j < j_z; ++j){
        reset_sample_infos(k_to, j);
      }
    }
    
    size_t get_start(size_t ti, size_t nthreads, size_t j_A, size_t j_Z){
      size_t n_js = j_Z - j_A;
      double t_fl = static_cast<double>(ti);
      size_t j_a = j_A + (t_fl / static_cast<double> (nthreads)) * n_js;
      return j_a;
    }

    size_t get_end(size_t ti, size_t nthreads, size_t j_A, size_t j_Z){
      size_t n_js = j_Z - j_A;
      double t_fl = static_cast<double>(ti);
      size_t j_z = j_A + ( ( t_fl + 1.) / static_cast<double> (nthreads)) * n_js;
      return j_z;
    }
    



    
    inline void base_put_sample_in_cluster(size_t i) {
      
      
      std::unique_ptr<double []> up_distances (new double [K]);
      double min_distance = std::numeric_limits<double>::max();
      size_t nearest_center = 0;
      double second_min_distance = std::numeric_limits<double>::max();
      
      for (size_t k = 0; k < K; ++k){
        set_center_sampleID_distance(k, i, second_min_distance, up_distances[k]);
        if (up_distances[k] <= second_min_distance){
          if (up_distances[k] < min_distance){
            second_min_distance = min_distance;
            nearest_center = k;
            min_distance = up_distances[k];
          }
          else{
            second_min_distance = up_distances[k];
          }
        }
      }
      final_push_into_cluster(i, nearest_center, min_distance, up_distances.get());
    }      
    
     inline void triangular_put_sample_in_cluster(size_t i, const double * const cc) {
      std::unique_ptr<double []> up_distances (new double [K]);
      double min_distance = std::numeric_limits<double>::max();
      double second_min_distance = std::numeric_limits<double>::max();

      size_t nearest_center = 0;
      
      for (size_t k = 0; k < K; ++k){
        up_distances[k] = std::numeric_limits<double>::max();
        if (cc[nearest_center*K + k] < second_min_distance + min_distance){
          set_center_sampleID_distance(k, i, second_min_distance, up_distances[k]);
          if (up_distances[k] < second_min_distance){
            if (up_distances[k] < min_distance){
              second_min_distance = min_distance;
              nearest_center = k;
              min_distance = up_distances[k];
            }
            else{
              second_min_distance = up_distances[k];
            }
          }
        }
      }
      final_push_into_cluster(i, nearest_center, min_distance, up_distances.get());
    }


    void test_parameters_to_tkai(size_t k0, size_t k1, size_t ndata_1, size_t ndata_2){
      if (k0 >= K || k1 < k0 || k1 > K){
        std::stringstream errm;
        errm << "invalid k0, k1 in kmeanspp : ";
        errm << "k0 = " << k0 << "     k1 = " << k1 << "     K = " << K << ".";
        throw std::runtime_error(errm.str());
      }
 
      if (ndata_1 != ndata_2){
        throw std::runtime_error("ndata_1 != ndata_2 in test_parameters_to_tkai");        
      }
    }


        inline size_t get_sample_from(std::vector<double> & v_cum_nearest_energies){
          
          if (v_cum_nearest_energies.back() == 0){
            throw std::runtime_error("exhausted in get_sample_from (kmeans++). in particular, the cumulative energy is zero. not sure what to do in this situation, could just return indices uniformly at random. but being cautious and throwing an error. please change it here if this is annoying. ");  
            /* alternatively, return uniformly with this : 
             * return dis(gen)%(v_cum_nearest_energies.size() - 1); */
          }
          
          return std::distance(
          v_cum_nearest_energies.begin(), 
          std::lower_bound(v_cum_nearest_energies.begin(), v_cum_nearest_energies.end(), dis_uni01(gen)*v_cum_nearest_energies.back())) - 1;
        }


            inline void kmpp_inner(size_t i, size_t k, double a_distance, P2Bundle & p2bun){
            
              if (a_distance < p2bun.d_2(i)){
                if (a_distance < p2bun.d_1(i)){
                  p2bun.d_2(i) = p2bun.d_1(i);
                  p2bun.k_2(i) = p2bun.k_1(i);
                  p2bun.d_1(i) = a_distance;
                  p2bun.k_1(i) = k;
                }
                else{
                  p2bun.d_2(i) = a_distance;
                  p2bun.k_2(i) = k;
                }
              }
            }



    
    
    /* where 
     * c_dt will store the centers
     * c_ind stores the indices (as per original ptr_datain).
     * cc stores the inter-ceter distances.
     * p2bun stores a1,d1,a2,d2,ori of data used for kmeanspp.
     * p2bun_d2 (const) is the data used for kmeanpp.
     * get centers [k0, k1).
     *  */

    void triangular_kmeanspp_after_initial(TData & c_dt, size_t * const c_ind, double * const cc, P2Bundle & p2bun, const TData * const ptr_p2bun_dt, size_t k0, size_t k1, bool from_full_data){
      
      if (from_full_data != (ptr_p2bun_dt == nullptr)){
        throw std::runtime_error("logic error in triangular_kmeanspp_after_initial : from_full_data != (ptr_p2bun_dt == nullptr)");
      }
      
      size_t kmeanspp_ndata = from_full_data ? ndata : ptr_p2bun_dt->get_ndata();
      
      if (kmeanspp_ndata -1 < k1 - k0){
        std::stringstream ss;
        ss << "triangular_kmeanspp_after_initial, attempting to find more centers than data - 1 : ";
        ss << "kmeanspp_ndata = " << kmeanspp_ndata << "   and   k1 - k0 = " << k1 - k0; 
        throw std::runtime_error(ss.str());        
      }
      
      test_parameters_to_tkai(k0, k1, p2bun.get_ndata(), kmeanspp_ndata);
      
      double a_distance;

      /* the cumulative sampling distribution */
      std::vector<double> v_cum_nearest_energies (kmeanspp_ndata + 1);      
      v_cum_nearest_energies[0] = 0.;
      
      if (k0 == 0){
        for (size_t i = 0; i < kmeanspp_ndata; ++i){
          v_cum_nearest_energies[i+1] = v_cum_nearest_energies[i] + 1.;
        }
      }

      else{
        for (size_t i = 0; i < kmeanspp_ndata; ++i){
          v_cum_nearest_energies[i+1] = v_cum_nearest_energies[i] + f_energy(p2bun.d_1(i));
        }
      }

      /* depending on whether we're using all the data or the data in ptr_p2bun_dt, distance calculations are performed differently. */      
      std::function<void(size_t, size_t)> set_distance_kk;
      std::function<void(size_t, size_t)> set_distance_ik;
      std::function<void(size_t)> update_c_ind_c_dt;

      if (from_full_data == true) {
  
        set_distance_kk = [this, &c_ind, &cc](size_t k, size_t kp){
          set_sampleID_sampleID_distance(c_ind[k], c_ind[kp], cc[k*K + kp]);  
        };
        
        set_distance_ik = [this, &a_distance, &p2bun, &c_ind](size_t i, size_t k){
          set_sampleID_sampleID_distance(c_ind[k], i, p2bun.d_2(i), a_distance);
        };
        
        update_c_ind_c_dt = [this, &c_ind, &c_dt, &v_cum_nearest_energies](size_t k){
          c_ind[k] = get_sample_from(v_cum_nearest_energies); 
          c_dt.append(ptr_datain->at_for_move(c_ind[k]));
        };
      }
      
      else{

        set_distance_kk = [this, &c_ind, &cc, &c_dt](size_t k, size_t kp){
          metric.set_distance(c_dt.at_for_metric(k), c_dt.at_for_metric(kp), cc[k*K + kp]);
        };
        
        set_distance_ik = [this, &a_distance, &p2bun, &c_ind, &c_dt, &ptr_p2bun_dt](size_t i, size_t k){
          metric.set_distance(c_dt.at_for_metric(k), ptr_p2bun_dt->at_for_metric(i), p2bun.d_2(i), a_distance);
        };
        
        update_c_ind_c_dt = [this, &c_ind, &c_dt, &v_cum_nearest_energies, &ptr_p2bun_dt, &p2bun](size_t k){
          
          size_t sami = get_sample_from(v_cum_nearest_energies);
          c_ind[k] = p2bun.ori(sami);
          c_dt.append(ptr_p2bun_dt->at_for_move(sami));
        };
      
      }
      
      
      
      /* k-means++, at last */
      for (size_t k = k0; k < k1; ++k){
        
        update_c_ind_c_dt(k);
 
        
        /* update cc */
        cc[k*K + k] = 0.;
        for (size_t kp = 0; kp < k; ++kp){
          set_distance_kk(k, kp);
          cc[kp*K + k] = cc[k*K + kp];
        }
 
        /* update nearest distances, second nearest distances, and centers */
        for (size_t i = 0; i < kmeanspp_ndata; ++i){
          /* note that if k0 = 0, k_1(i) is 0 so this is fine. */
          if (cc[p2bun.k_1(i)*K + k] < p2bun.d_1(i) + p2bun.d_2(i)){
            set_distance_ik(i,k);
            kmpp_inner(i, k, a_distance, p2bun);
          }
          v_cum_nearest_energies[i+1] = v_cum_nearest_energies[i] + f_energy(p2bun.d_1(i));
        }
      }
    }


    

    /* all parameters passed in are to be set. TODO : parallelisation */
    inline void triangular_kmeanspp(size_t * const centers, double * const cc, P2Bundle & p2bun, TData & c_dt){
      
      for (size_t i = 0; i < ndata; ++i){
        p2bun.ori(i) = i;
      }
      
      size_t k0 = 0;
      size_t k1 = K;
      triangular_kmeanspp_after_initial(c_dt, centers, cc, p2bun, nullptr, k0, k1, true);
    }




    
        
    inline void triangular_kmeanspp_aq2(size_t * const centers, double * const cc, P2Bundle & p2bun, TData & c_dt, size_t n_bins){
       
       
      /* experiments so far show that multithreading does not help here, can hurt. what's weird is that even if nthreads = 1 in 
       * the pll version, it's sig slower than the serial version.  */
      bool multithread_kmpp = false;
      
      double a_distance;

      /* non_tail_k will be how many k's for which only 1/n_bins of the data is used.
       * it will be (n_bins - 1)/n_bins, rounded down to the nearest multiple of n_bins.
       * an analysis suggests that the tail should be sqrt(K/nbins), so this is conservative.
       * we also ensure that the tail is at least of length min(K, 50). 
       * */
      size_t non_tail_k;
      non_tail_k = static_cast<size_t>(static_cast<double>(K)*(1. - 1./static_cast<double>(n_bins)));
      size_t minimal_tail_k = std::min<size_t>(K, 50);
      size_t maximal_non_tail_k = K - minimal_tail_k;
      non_tail_k = std::min(non_tail_k, maximal_non_tail_k);
                  
            
      /* tail_k is how many k's have use all the data (at the end) */
      size_t tail_k = K - non_tail_k;
      
      /* the data will be copied randomly into bins */
      std::vector<TData> p2buns_dt (n_bins, {*ptr_datain, true});
      std::vector<std::vector<size_t>> original_indices(n_bins);
      std::vector<P2Bundle> p2buns(n_bins);

      
      std::vector<size_t> bin_indices (n_bins);
      for (size_t bin = 0; bin < n_bins; ++bin){
        bin_indices[bin] = bin;
      }
      for (size_t i = 0; i < ndata; ++i){
        if (i%n_bins == 0){
          /* shuffle the indices */
          for (size_t bin = 0; bin < n_bins; ++bin){
            std::swap(bin_indices[bin], bin_indices[bin + dis(gen)%(n_bins - bin)]);
          }
        } 
        size_t bin = bin_indices[i%n_bins];
        p2buns_dt[bin].append(ptr_datain->at_for_move(i));
        original_indices[bin].push_back(i);
      }
      
      for (size_t bin = 0; bin < n_bins; ++bin){
        p2buns[bin] = P2Bundle(original_indices[bin]);
      }

        auto update_nearest_info = [&p2buns, &p2buns_dt, &c_dt, &cc, &a_distance, this](size_t bin, size_t k_start, size_t k_end, size_t i0, size_t i1)
        {
          for (size_t i = i0; i < i1; ++i){
            for (size_t k = k_start; k < k_end; ++k){
              if (cc[p2buns[bin].k_1(i)*K + k] < p2buns[bin].d_1(i) + p2buns[bin].d_2(i)){
                metric.set_distance(c_dt.at_for_metric(k), p2buns_dt[bin].at_for_metric(i), p2buns[bin].d_2(i), a_distance);
                kmpp_inner(i, k, a_distance, p2buns[bin]);
              }
            }
          }
        };
      

      for (size_t bin = 0; bin < n_bins; ++bin){
        size_t k0 = (bin*non_tail_k)/n_bins;
        size_t k1 = ((bin + 1)*non_tail_k)/n_bins;        
        

        /* update nearest info from 0 to k0 of this bin*/        
        if (multithread_kmpp == false){
          update_nearest_info(bin, 0, k0, 0, p2buns[bin].get_ndata());
        }

        else{
          std::vector<std::thread> threads;
          for (size_t ti = 0; ti < get_nthreads(); ++ti){
            threads.emplace_back([this, ti, bin, &p2buns, k0, &update_nearest_info] () {
              update_nearest_info(
              bin, 0, k0, 
              get_start(ti, get_nthreads(), 0, p2buns[bin].get_ndata()),
              get_end(ti, get_nthreads(), 0, p2buns[bin].get_ndata()));});
          }
          
          for (auto & t : threads){
            t.join();
          }
        }
                
        
        /* we don't even try to parallelize to center by center kmeans++ part, too many syncs. */
        triangular_kmeanspp_after_initial(c_dt, centers, cc, p2buns[bin], &p2buns_dt[bin], k0, k1, false);
      }
      
      /* update all until the tail k's (*&^*) */
      for (size_t bin = 0; bin < n_bins; ++bin){
        size_t k1 = ((bin + 1)*non_tail_k)/n_bins;

        
        if (multithread_kmpp == false){
          update_nearest_info(bin, k1, non_tail_k, 0, p2buns[bin].get_ndata());
        }
        
        else{
          std::vector<std::thread> threads;
          for (size_t ti = 0; ti < get_nthreads(); ++ti){
            threads.emplace_back([this, ti, bin, &p2buns, k1, non_tail_k, & update_nearest_info] () {update_nearest_info(
              bin, k1, non_tail_k, 
              get_start(ti, get_nthreads(), 0, p2buns[bin].get_ndata()),
              get_end(ti, get_nthreads(), 0, p2buns[bin].get_ndata()));});
          }
          
          for (auto & t : threads){
            t.join();
          }
        }
      }
      
      
      //get p2bun up to date. 
      for (size_t bin = 0; bin < n_bins; ++bin){
        for (size_t i = 0; i < p2buns[bin].get_ndata(); ++i){
          size_t index0 = p2buns[bin].ori(i);
          p2bun.d_1(index0) = p2buns[bin].d_1(i);
          p2bun.d_2(index0) = p2buns[bin].d_2(i);          
          p2bun.k_1(index0) = p2buns[bin].k_1(i);
          p2bun.k_2(index0) = p2buns[bin].k_2(i);
          p2bun.ori(index0) = index0;
        }
      }
      
      //set the tail k's
      triangular_kmeanspp_after_initial(c_dt, centers, cc, p2bun, nullptr, K - tail_k, K, true);

    }
    
        
    std::string get_base_summary_string(){
      std::ostringstream out;
      out  << round << "\tmE: " <<  std::setprecision(7) << E_total / static_cast<double>(ndata) << std::setprecision(5) <<  
      "\tatime: "  <<  time_to_initialise_centers/1000  << 
      "  itime: "  <<  time_initialising/1000  <<  
      "  ctime: "  <<  time_in_update_centers/1000  <<  
      "\tutime: "  <<  time_in_update_sample_info/1000  <<  
      "\trtime: "  <<  time_in_redistribute/1000  <<  
      "\tttime: "  <<  time_total/1000  <<   
      "\tlg2 nc(c): "  << std::log2(ncalcs_in_update_centers) <<  
      "\tlg2 nc: "  <<  std::log2(ncalcs_total)  <<  
      "\t pc: " << static_cast<double> (metric.get_rel_calccosts())  ;
      
      //   <<  "\tstime: "  <<  time_in_update_all_cluster_statistics/1000
      
      return out.str();
    
      //return "*" + std::to_string(round) + "\t E: " + std::to_string(E_total) + "\t itime: " + std::to_string( time_initialising/1000 ) + "\t ctime: " + std::to_string( time_in_update_centers/1000 ) + "\t utime: " + std::to_string( time_in_update_sample_info/1000 ) + "\t rtime: " + std::to_string( time_in_redistribute/1000 ) + "\t ttime: " + std::to_string( time_total/1000 )  + "\t stime: " + std::to_string( time_in_update_all_cluster_statistics/1000 ) +  "\t log2 nc(c): " + std::to_string(std::log2(ncalcs_in_update_centers) )+ "\t log2 nc: " + std::to_string( std::log2(ncalcs_total) );
    }
        
    inline size_t get_a1(size_t k, size_t j) const{
      return nearest_1_infos[k][j].a_x;
    }
    
    inline double get_d1(size_t k, size_t j){
      return nearest_1_infos[k][j].d_x;
    }
    
    inline double get_e1(size_t k, size_t j){
      return nearest_1_infos[k][j].e_x;
    }
    
    inline double get_e1_tail(size_t k){
      return nearest_1_infos[k].back().e_x;
    }
    
    inline size_t get_ndata(size_t k){
      return cluster_datas[k].get_ndata();
    }
    
    inline size_t get_nthreads(){
      return nthreads;
    }
    
    inline double get_nthreads_fl(){
      return nthreads_fl;
    }
    
    
    
    /* rule : functions with suffix 'basic' will not touch to_leave_cluster */
    inline void reset_nearest_info_basic(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest){
      nearest_1_infos[k][j].reset(k_nearest, d_nearest, e_nearest);      
      cluster_has_changed[k] = true;
    }
    
    inline void reset_nearest_info(size_t k, size_t j, size_t k_nearest, double d_nearest, double e_nearest){
      reset_nearest_info_basic(k, j, k_nearest, d_nearest, e_nearest);      
      if (k != k_nearest){
        std::lock_guard<std::mutex> lock(mutex0);
        to_leave_cluster[k].push_back(j);
      }
    }
    
    
    
    inline void swap_center_with_sample(size_t k, size_t j){
      cluster_has_changed[k] = true;
      center_IDs[k] = sample_IDs[k][j]; //just added.
      
      nszen::swap<TData>(centers_data, k, cluster_datas[k], j);
      
      
      reset_sample_infos_basic(k, j);
    }
          
    
    inline double get_E_total(){
      return E_total;
    }
    
    inline double get_cluster_energy(size_t k){
      return cluster_energies[k];
    }
    
    void move_center_into_its_own_cluster(size_t k){
      //move it from the original data (more likely to have a cache miss, but for now less code = good)
      //note that we use the base_version of the function, at the overloaded version (for clarans l23) uses cc which might not
      //be reliable at this point.
      base_put_sample_in_cluster(center_IDs[k]);
      center_IDs[k] = 911911; 
      cluster_has_changed[k] = true;
    }
    
    void overwrite_center_with_sample(size_t k1, size_t k2, size_t j2){
      centers_data.replace_with(k1, cluster_datas[k2].at_for_move(j2));
      center_IDs[k1] = sample_IDs[k2][j2];
      cluster_has_changed[k1] = true;
    }

    //remove the j'th sample from cluster k, and (if it is not the last element) fill the hole with the tail (size drops by 1) 
    inline void remove_with_tail_pull(size_t k, size_t j){
      if (j != get_ndata(k) - 1){
        nearest_1_infos[k][j] = *(nearest_1_infos[k].end() - 1);
        sample_IDs[k][j] = *(sample_IDs[k].end() - 1);
        replace_with_last<TData> (cluster_datas[k], j);
        custom_replace_with_last(k,j);
      }
      
      nearest_1_infos[k].pop_back();
      sample_IDs[k].pop_back();
      cluster_datas[k].remove_last();
      custom_remove_last(k);
      
      cluster_has_changed[k] = true;            
    }

    size_t draw_k_uniform(){
      return dis(gen)%K;
    }
    
    size_t draw_k_prop_ndata(){
      std::unique_ptr<size_t []> up_cum_ndatas (new size_t [K]);
      size_t cum_ndata = 0;
      for (size_t k = 0; k < K; ++k){
        cum_ndata += get_ndata(k);
        up_cum_ndatas[k] = cum_ndata;
      }
      if (up_cum_ndatas[K-1] != (ndata - K)){
        mowri << "Weird : cum_ndatas[K-1] != ndata : up_cum_ndatas[K-1] = " << up_cum_ndatas[K-1] << " and ndata = " << ndata << zentas::Endl;
        throw std::logic_error("(see above)");
      }
      size_t i = dis(gen)%(ndata - K);
      size_t k_below = 0;
      while (i >= up_cum_ndatas[k_below]){
        ++k_below;
      }
 
      return k_below;
    }
    
    inline void signal_cluster_change(size_t k){
      cluster_has_changed[k] = true;
    }
    
    size_t draw_j_uniform(size_t k){
      return dis(gen)%get_ndata(k);
    }
    
    //set the distance from center k to the j1'th element of cluster k1.
    inline void set_center_sample_distance(size_t k, size_t k1, size_t j1, double & distance) {
      metric.set_distance(centers_data.at_for_metric(k), cluster_datas[k1].at_for_metric(j1), distance);
    }

    inline void set_center_sample_distance(size_t k, size_t k1, size_t j1, double threshold, double & distance) {
      metric.set_distance(centers_data.at_for_metric(k), cluster_datas[k1].at_for_metric(j1), threshold, distance);
    }
        
    inline void set_center_sampleID_distance(size_t k, size_t i, double threshold, double & distance) {
      metric.set_distance(centers_data.at_for_metric(k), ptr_datain->at_for_metric(i), threshold, distance);
    }

    inline void set_sampleID_sampleID_distance(size_t i1, size_t i2, double threshold, double & distance) {
      metric.set_distance(ptr_datain->at_for_metric(i1), ptr_datain->at_for_metric(i2), threshold, distance);
    }

    inline void set_sampleID_sampleID_distance(size_t i1, size_t i2, double & distance) {
      metric.set_distance(ptr_datain->at_for_metric(i1), ptr_datain->at_for_metric(i2), distance);
    }
        
        
    inline double get_center_sample_distance(size_t k, size_t k1, size_t j1) {
      double distance;
      set_center_sample_distance(k, k1, j1, distance);
      return distance;
    }
    
    inline double get_sample_sample_distance(size_t k , size_t j1, size_t j2) {
      double adistance;
      metric.set_distance(cluster_datas[k].at_for_metric(j1), cluster_datas[k].at_for_metric(j2), adistance);
      return adistance;
    }
    
    inline void set_sample_sample_distance(size_t k1, size_t j1, size_t k2, size_t j2, double threshold, double & adistance) {
      metric.set_distance(cluster_datas[k1].at_for_metric(j1), cluster_datas[k2].at_for_metric(j2), threshold, adistance);
    }
    
    
    inline void set_center_center_distance(size_t k1, size_t k2, double threshold, double & adistance) {
      metric.set_distance(centers_data.at_for_metric(k1), centers_data.at_for_metric(k2), threshold, adistance);
    }
    inline void set_center_center_distance(size_t k1, size_t k2, double & adistance) {
      metric.set_distance(centers_data.at_for_metric(k1), centers_data.at_for_metric(k2), adistance);
    }
    

    
    void print_centers(){
      mowri << centers_data.get_string() << zentas::Flush;
    }

    
    void print_ndatas(){
      for (size_t k = 0; k < K; ++k){
        mowri << get_ndata(k) << " ";
      }
      mowri << zentas::Endl;
    }
    
    /* for inter-center distances etc. if relevant, otherwise just make {}. */
    virtual void set_center_center_info() = 0;
    virtual void update_center_center_info() = 0;
    

    
    void default_initialise_with_kmeanspp(){

      unsigned n_bins;

      if (initialisation_method == "kmeans++"){
        n_bins = 1;
      }
      
      else{
        std::string prefix = "kmeans++-";
        n_bins = extract_INT(initialisation_method, prefix.size());
        if (n_bins == 0){
          throw std::runtime_error("n_bins passed to kmeans++ should be positive.");
        }
      }

      
      up_kmoo_bundle->reset(K, ndata);
      
      
      if (n_bins == 1){
        triangular_kmeanspp(center_indices_init, up_kmoo_bundle->cc.data(), up_kmoo_bundle->p2bun, up_kmoo_bundle->c_dt);
      }
      
      else{
        triangular_kmeanspp_aq2(center_indices_init, up_kmoo_bundle->cc.data(), up_kmoo_bundle->p2bun, up_kmoo_bundle->c_dt, n_bins);
      }
      
      /* TODO here: rearrange everything, so that center_indices_init are in order. what a palava. */
      std::vector<std::array<size_t, 2>> vi (K);
      for (size_t k = 0; k < K; ++k){
        std::get<0>(vi[k]) = center_indices_init[k];
        std::get<1>(vi[k]) = k;
      }
      
      auto fg = std::less<size_t>();
      std::sort(vi.begin(), vi.end(), [&fg](std::array<size_t, 2> & lhs, std::array<size_t, 2> & rhs){ return fg(std::get<0>(lhs), std::get<0>(rhs)); });
      //for (size_t k = 0; k < K; ++k){
        //std::cout << "k: " << k << "  i:" << std::get<0>(vi[k]) << " kori:" << std::get<1>(vi[k]) << std::endl;
      //}
      
      std::vector<size_t> old_to_new (K);
      for (unsigned k = 0; k < K; ++k){
        center_indices_init[k] = std::get<0>(vi[k]);
        old_to_new[std::get<1>(vi[k])] = k;
      }
      
      for (unsigned i = 0; i < ndata; ++ i){
        up_kmoo_bundle->p2bun.k_1(i) = old_to_new[up_kmoo_bundle->p2bun.k_1(i)];
        up_kmoo_bundle->p2bun.k_2(i) = old_to_new[up_kmoo_bundle->p2bun.k_2(i)];
      }
      
      //std::cout << "done the " << K << " centers" << std::endl;
      //std::abort();
      
      
       
    }
    
    void go(){

      bool with_tests = false;
      if (with_tests == true){
        mowri << "\n\nCOMPILED WITH TESTS ENABLED : WILL BE SLOW" <<zentas::Endl;
      }


      /*TODO: prevent code duplication of this string (in pyzentas.pyx and here) */      
      mowri << 
R"(
(The prevent output to terminal, set capture_output to false)

The output below contains the following round-by-round statistics
------------------------------------------------------------------------------------------------------------------------------------
first column  : round (where a round is defined by center change)
mE            : mean energy over samples
itime         : time [in milliseconds] taken for initialisation (first assignments)
ctime         : time spent in center update. For clarans, this is the time spent evaluating proposals. 
                For Voronoi, this is time in updating assignments
utime         : time spent in updating. For clarans, this is the time spent determining the nearest and second nearest 
                centers following a swap. For voronoi, this is time spent determining medoids.
rtime         : time spent implementing the center move. This cost involves moving data between clusters while maintaining 
                it in a random order. If rooted = True, this can be expected be higher that when rooted = False,
                (with a correspondlingly lower ctime for rooted = True)
ttime         : total time.
lg2 nc(c)     : log based 2 of the number of distance calculations made in center update (corresponding to ctime)
lg2 nc        : log based 2 of total number of distance calculations made.
pc            : distance calculations can be terminated early if it can be established that they are above some threshold. 
                This is particularly effective when the Levenshtein or normalised Levenshtein distance is used.
                pc measures how effective the early stopping in a distance computation is. For details see the appendix of our paper
nprops        : for clarans, the number of rejected proposals before one is accepted.       
------------------------------------------------------------------------------------------------------------------------------------
)";


      





      tstart_initialise_centers = std::chrono::high_resolution_clock::now(); 
      
      /* initialisation from indices. */
      if (initialisation_method == "from_indices_init"){
        //already done in constructor
      }
      
      /* initialisation uniformly. */
      else if (initialisation_method == "uniform"){
        populate_uniformly(center_indices_init, K, ndata, dis, gen);
      }
      
      /* initialisation uniformly according to Bachem et al. */
      else if (initialisation_method.substr(0,8) == "afk-mc2-"){
        populate_afk_mc2<TMetric, DataIn>(initialisation_method, center_indices_init, *ptr_datain, metric, K, ndata, mowri, gen, dis, f_energy);
      }
      
      /* initialisation with kmeans++ TODO : in dev mode */
      else if (initialisation_method == "kmeans++" || initialisation_method.substr(0, 9) == "kmeans++-"){
        initialise_with_kmeanspp();
      }
      
      else {
        std::stringstream vims_ss;
        vims_ss << "The valid strings for initialisation_method are [from_indices_init, uniform, kmeans++-INT, afk-mc2-INT (for some positive INT)]";
        throw std::runtime_error(vims_ss.str());
      }
  
      std::sort(v_center_indices_init.begin(), v_center_indices_init.end());
      std::copy(center_indices_init, center_indices_init + K, center_IDs);

      for (size_t k = 0; k < K; ++k){
        centers_data.append(ptr_datain->at_for_move(center_indices_init[k]));
        /* here we make an empty cluster, using datain to determine nec. `shape' parameters */
        cluster_datas.emplace_back(*ptr_datain, true);
      }
  
      /* if the indices are from the user or in debug mode, we run a test that initialising indices look fine. */
      if (with_tests || initialisation_method == "from_indices_init"){
        post_initialise_centers_test();
      }
      
      auto t_endit = std::chrono::high_resolution_clock::now();
      time_to_initialise_centers = std::chrono::duration_cast<std::chrono::microseconds>(t_endit - tstart_initialise_centers).count();

      time_total = time_initialising;
      
      
      
      
      
      


      
      std::chrono::time_point<std::chrono::high_resolution_clock> t0;
      std::chrono::time_point<std::chrono::high_resolution_clock> t1;
      std::chrono::time_point<std::chrono::high_resolution_clock> t2;
      std::chrono::time_point<std::chrono::high_resolution_clock> t3;
      std::chrono::time_point<std::chrono::high_resolution_clock> t4;
      tstart = std::chrono::high_resolution_clock::now();
      size_t ncalcs0;
      size_t ncalcs1;
      size_t ncalcs2;
      
  
      //initialisation
      set_center_center_info();
      put_samples_in_clusters();                  
      set_all_cluster_statistics();
      
      /* (CHECK POINT) all assignments and cluster statistics must be correct. Tests to pass : all */
      if (with_tests == true){
        post_initialisation_test();
      }
      

      //mowri << "round :" << -1 << "\t tenergy : " <<  E_total << " \t itime : " << time_initialising << " \t ctime : " << time_in_update_centers << " \t utime : " << time_in_update_sample_info << " \t rtime : " << time_in_redistribute << " \t ttime : " << time_total << zentas::Endl;        
  
      t1 = std::chrono::high_resolution_clock::now();
      time_initialising = std::chrono::duration_cast<std::chrono::microseconds>(t1 - tstart).count();
      ncalcs_initialising = metric.get_ncalcs();
      
      round_summary();
      
      //mowri << time_total << "   " << maxtime MM zentas::Endl;
      while ((time_total < maxtime_micros) && (round < maxrounds) && (E_total / static_cast<double>(ndata)) >= minmE){ // (E_total != old_E_total) && 

        t0 = std::chrono::high_resolution_clock::now();
        ncalcs0 = metric.get_ncalcs();
        
        bool modified_centers = update_centers(); 
        
        if (modified_centers == true){
          update_center_center_info();
          /* (CHECK POINT) the above should swap center data with cluster data. It must set cluster_has_changed[k] to true if cluster k changes
           * It may change a_x etc. Tests to pass : injective_ID_test. */
          if (with_tests == true){
            post_center_update_test();
          }
        }
             
        t1 = std::chrono::high_resolution_clock::now();
        ncalcs1 = metric.get_ncalcs();
        time_in_update_centers += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        ncalcs_in_update_centers += ncalcs1 - ncalcs0;

          

        
        if (modified_centers == false){
          time_total = std::chrono::duration_cast<std::chrono::microseconds>(t1 - tstart).count();
          ncalcs_total = metric.get_ncalcs();
          round_summary();
          
          break;
        }
        
        else{
          
          
          for (auto & x : to_leave_cluster){
            x.resize(0);
          }
          
          update_sample_info(); 
          /* (CHECK POINT) all sample info (a_x etc) must be correct, although a_x can disagree with k. to_leave_cluster must be correctly set */
          if (with_tests == true){
            post_sample_update_test();
          }
          
          t2 = std::chrono::high_resolution_clock::now();
          ncalcs2 = metric.get_ncalcs();
          time_in_update_sample_info += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
          ncalcs_in_update_sample_info += ncalcs2 - ncalcs1;
          
          redistribute();                
          /* (CHECK POINT) a_x must correspond to k. Tests to pass : as_assigned_test. */
          if (with_tests == true){
            post_redistribute_test();
          }
          
     
          t3 = std::chrono::high_resolution_clock::now();
          time_in_redistribute += std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
          
          
          
          
          update_all_cluster_statistics();
          /* (CHECK POINT) all cluster statistics (cluster_energies etc) must be correct. Tests to pass : all */
          if (with_tests == true){
            post_update_statistics_test();
          }
  
          t4 = std::chrono::high_resolution_clock::now();
          time_in_update_all_cluster_statistics += std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();  
          
          time_total = std::chrono::duration_cast<std::chrono::microseconds>(t4 - tstart).count();
          ncalcs_total = metric.get_ncalcs();
          round_summary();
          
          
          ++round;        
        }
      }
      
      if (time_total >= maxtime_micros){
        mowri << "exceeded maxtime (" << maxtime_micros/1000. << ")" << zentas::Endl;
      }
      
      if (round >= maxrounds){
        mowri << "exceeded maxrounds (" << maxrounds << ")" << zentas::Endl;
      }
      
      if (E_total / static_cast<double>(ndata) < minmE){
        mowri << "mE below termination minmE (" << minmE << ")" << zentas::Endl;        
      }
      
      
      //populate labels : 
      for (size_t k = 0; k < K; ++k){
        labels[center_IDs[k]] = k;
        for (size_t j = 0; j < get_ndata(k); ++j){
          labels[sample_IDs[k][j]] = k;
        }
      }
    }

    

    /* basic, because has no effect on to_leave_cluster */
    void reset_sample_infos_basic(size_t k, size_t j){
      
      std::unique_ptr<double []> up_distances (new double [K]);
      double min_distance = std::numeric_limits<double>::max();
      size_t nearest_center = 0;
      double second_min_distance = std::numeric_limits<double>::max();

      for (size_t kp = 0; kp < K; ++kp){
        set_center_sample_distance(kp, k, j, second_min_distance, up_distances[kp]);
        if (up_distances[kp] < second_min_distance){
          if (up_distances[kp] < min_distance){
            nearest_center = kp;
            min_distance = up_distances[kp];
          }
          else{
            second_min_distance = up_distances[kp];
          }
        }
      }
      
      reset_nearest_info_basic(k,j, nearest_center, min_distance, f_energy(min_distance));
      reset_sample_custom(k,j, nearest_center, up_distances.get());
    }
    
    void reset_sample_infos(size_t k, size_t j){
      //mowri << k << "--" << j << " " << zentas::Flush;
      reset_sample_infos_basic(k,j);
      if (k != nearest_1_infos[k][j].a_x){
        std::lock_guard<std::mutex> lock (mutex0);
        to_leave_cluster[k].push_back(j); //haha!
      }
    }
    


/////////////////////////////////////
/////////// private /////////////////
/////////////////////////////////////    
  private:

    virtual void put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) = 0;

    virtual void initialise_with_kmeanspp() = 0;
        
    /* called from put_sample_in_clusters (virtual inline functions are not nonsensical) */    
    //virtual inline void put_sample_custom_in_cluster(size_t i, size_t a_x, const double * const distances) = 0;
    virtual inline void put_sample_custom_in_cluster(size_t, size_t, const double * const) {}

    
    /* when implementing this, one must guarantee that cluster_has_changed is correctly set in this function */
    /* lowest and second lowest values in distances are reliable */
    //virtual inline void reset_sample_custom(size_t k, size_t j, size_t nearest_center, const double * const distances) = 0;
    virtual inline void reset_sample_custom(size_t, size_t, size_t, const double * const) {}
    
    virtual bool update_centers() = 0; 
    virtual void update_sample_info() = 0;
      
    /* functions called from within redistribute */
    //virtual inline void custom_append(size_t k_new, size_t k, size_t j) = 0;
    virtual inline void custom_append(size_t, size_t, size_t) {}
    
    //virtual inline void custom_replace_with_last(size_t k, size_t j) = 0;
    virtual inline void custom_replace_with_last(size_t, size_t) {}
    
    //virtual inline void custom_replace_with(size_t k1, size_t j1, size_t k2, size_t j2) = 0;
    virtual inline void custom_replace_with(size_t, size_t, size_t, size_t) {};
    
    //virtual inline void custom_remove_last(size_t k) = 0;
    virtual inline void custom_remove_last(size_t) {};

    /* functions called from within set_cluster_statistics */
    //virtual inline void increment_custom_cluster_statistics(size_t k, size_t j) = 0;
    virtual inline void increment_custom_cluster_statistics(size_t, size_t) {};

    virtual void set_normalised_custom_cluster_statistics(size_t k) = 0;
    virtual void set_to_zero_custom_cluster_statistics(size_t k) = 0;

    virtual void round_summary() = 0;    


        
    void put_samples_in_clusters(){
    
      // get sorted center IDs.
      std::unique_ptr<size_t []> up_sorted_center_IDs (new size_t [K]);
      auto sorted_center_IDs = up_sorted_center_IDs.get();
      std::copy(center_IDs, center_IDs + K, sorted_center_IDs);
      std::sort(sorted_center_IDs, sorted_center_IDs + K);
      
      // make array with non center IDs
      std::vector<size_t> non_center_IDs (ndata - K);
      size_t non_center_number = 0;
      size_t center_number = 0;
      
      for (size_t ID = 0; ID < ndata; ++ID){
        if (center_number == K || ID != sorted_center_IDs[center_number]){
          non_center_IDs[non_center_number] = ID;
          ++non_center_number;
        }
        
        else{
          ++center_number;
        }
      }
      
      for (size_t i = ndata - K - 1 ; i > 0; --i){
        size_t swap_with = dis(gen)%(i+1);
        std::swap(non_center_IDs[i], non_center_IDs[swap_with]);
      }
      
      
      std::vector<std::thread> threads;      
      for (size_t ti = 0; ti < nthreads; ++ti){
        threads.emplace_back([this,ti,  &non_center_IDs] () { pll_put_samples_in_cluster(get_start(ti, get_nthreads(), 0, ndata - K), get_end(ti, get_nthreads(), 0, ndata - K), non_center_IDs); } );
      }
      for (auto & t : threads){
        t.join();
      }

    }
      
    void pll_put_samples_in_cluster(size_t i_a, size_t i_z, std::vector<size_t> & non_center_IDs){

      std::string prefix = "kmeans++";
      if (initialisation_method.substr(0, prefix.size()) == prefix){
        //use the kmeans++ initialisation bundle. TODO here. 
        for (size_t i = i_a; i < i_z; ++i){
          size_t nc_ID = non_center_IDs[i];
          size_t k1 = up_kmoo_bundle->p2bun.k_1(nc_ID); 
          size_t k2 = up_kmoo_bundle->p2bun.k_2(nc_ID);
          double d1 = up_kmoo_bundle->p2bun.d_1(nc_ID);
          double d2 = up_kmoo_bundle->p2bun.d_2(nc_ID);
          
          //size_t oripotter = up_kmoo_bundle->p2bun.ori(nc_ID);
          
          
          //std::cout << " (i)" << i << " (nc_ID)" << nc_ID << " (oripotter)" << oripotter << std::endl;
          final_push_into_cluster_post_kmeanspp(nc_ID, k1, k2, d1, d2);

        }
      }      

      else{        
        for (size_t i = i_a; i < i_z; ++i){
          put_sample_in_cluster(non_center_IDs[i]);
        }
      }
      
      up_kmoo_bundle.reset();
    }
  


    void put_samples_in_clusters_old_linear(){
    
      std::unique_ptr<size_t []> up_sorted_center_IDs (new size_t [K]);
      auto sorted_center_IDs = up_sorted_center_IDs.get();
      std::copy(center_IDs, center_IDs + K, sorted_center_IDs);
      std::sort(sorted_center_IDs, sorted_center_IDs + K);
      
      size_t sorted_ID_rank = 0;
      for (size_t i = 0; i < this->ndata; ++i){
        //check if data point i is a center, if so : it does not get added to cluster
        if (sorted_ID_rank < K && i == sorted_center_IDs[sorted_ID_rank]){
          ++sorted_ID_rank;
        }
        
        else{
          put_sample_in_cluster(i);
        }
      }
    }
        
    /* set energy of cluster k, and set as custom cluster statistics */
    void set_cluster_statistics(size_t k){
      cluster_energies[k] = 0;
      set_to_zero_custom_cluster_statistics(k);
      for (size_t j = 0; j < get_ndata(k); ++j){
        cluster_energies[k] += nearest_1_infos[k][j].e_x;
        increment_custom_cluster_statistics(k,j); // includes max
      }
      cluster_mean_energies[k] = cluster_energies[k] / static_cast<double> (get_ndata(k));
      set_normalised_custom_cluster_statistics(k);
    }
    
    void set_all_cluster_statistics(){
      E_total = 0;
      for (size_t k = 0; k < K; ++k){
        set_cluster_statistics(k);
        cluster_has_changed[k] = false;
        E_total += cluster_energies[k];
      }
    }
    
    void update_all_cluster_statistics(){
      old_E_total = E_total;
      E_total = 0;
      for (size_t k = 0; k < K; ++k){
        if (cluster_has_changed[k] == true){ 
          set_cluster_statistics(k);
          cluster_has_changed[k] = false;
        }
        
        //mowri << k << " ----- ( " << center_IDs[k] << " )  ---- energy : "  << cluster_energies[k] << zentas::Endl;
        
        //for (size_t j = 0; j < get_ndata(k); ++j){
          //mowri << sample_IDs[k][j] << "  (" <<nearest_1_infos[k][j].e_x << ")  " << zentas::Flush;
        //}
        //mowri << zentas::Endl;
        


        E_total += cluster_energies[k];
      }
      
      //mowri << "\n";
    }
    

    
    virtual void put_sample_in_cluster(size_t i) = 0;
    
    
    
    // For clarans, it is preferential to redistribute from k_to first (TODO : detailed explanation of what redistribute() does). This function allows clarans to set the order of redistribution. 
    virtual void set_redistribute_order(std::vector<size_t> & redistribute_order) = 0;
    

    
    
    // single threaded redistribution. TODO : make multithreaded.
    // requires that to_leave_cluster be reliably set. 
    // centers are unchanged by this function. this function simply moves samples between clusters.
    // additional complexity required to maintain randomness (samples inserted into random indices)
    void redistribute(){
      size_t k_new;
      size_t j;
      
      size_t j_new;
      
      std::vector<size_t> redistribute_order (K, 0);
      set_redistribute_order(redistribute_order);
      
      
      std::vector<size_t> redistribute_rank (K,0);
      for (size_t k = 0; k < K; ++k){
        redistribute_rank[redistribute_order[k]] = k;
      }
      
      
      bool target_is_a_mover;
      for (auto & k : redistribute_order){ //size_t k = 0; k < K; ++k){
        if (nthreads > 1){
          std::sort(to_leave_cluster[k].begin(), to_leave_cluster[k].end());
        }
        // it would be nice if one could auto over a vector in reverse
        for (size_t ji = to_leave_cluster[k].size(); ji-- > 0;){
          
          j = to_leave_cluster[k][ji];
          
          if (j != std::numeric_limits<size_t>::max()){
            k_new = nearest_1_infos[k][j].a_x;
            
            cluster_has_changed[k] = true;
            cluster_has_changed[k_new] = true;
            
            // now for the redistribution 
            // (1) make space somewhere in k_new (k,j) will be inserted at k_new, j_new:
            j_new = dis(gen)%(get_ndata(k_new) + 1);
            
            size_t insert_attempts = 0;
            while (j_new != get_ndata(k_new) && nearest_1_infos[k_new][j_new].a_x != k_new && insert_attempts < 10){
              j_new = dis(gen)%(get_ndata(k_new) + 1);
              ++insert_attempts;
            }
            
            //if (insert_attempts == 10){
              //mowri << "from : "  << k <<  " k_to " << redistribute_order[0] << " k_new : " << k_new << " j_new : " << j_new << zentas::Endl;
              ////throw std::runtime_error("insert attempts is 10. Wow");
            //}
            
            /* case 1 : k,j goes on the tail */
            if (j_new == get_ndata(k_new)){
              nearest_1_infos[k_new].push_back(nearest_1_infos[k][j]);
              sample_IDs[k_new].push_back(sample_IDs[k][j]);
              cluster_datas[k_new].append(cluster_datas[k].at_for_move(j));
              custom_append(k_new, k, j);
            }
            
            /* case 2 : k_new, j_new goes on the tail, then k,j goes in at k_new, j_new */
            else{
              
              if (nearest_1_infos[k_new][j_new].a_x != k_new){
                target_is_a_mover = true;
              }
              
              else{
                target_is_a_mover = false;
              }
              
              nearest_1_infos[k_new].push_back(nearest_1_infos[k_new][j_new]);
              sample_IDs[k_new].push_back(sample_IDs[k_new][j_new]);
              cluster_datas[k_new].append(cluster_datas[k_new].at_for_move(j_new));
              custom_append(k_new, k_new, j_new);
              
              nearest_1_infos[k_new][j_new] = nearest_1_infos[k][j];
              sample_IDs[k_new][j_new] = sample_IDs[k][j];
              cluster_datas[k_new].replace_with(j_new, cluster_datas[k].at_for_move(j));
              custom_replace_with(k_new, j_new, k,j);
              
              
              
              if (target_is_a_mover && redistribute_rank[k_new] > redistribute_rank[k]){ //quadratic : ouch. sort this out. 
                for (auto & x : to_leave_cluster[k_new]){
                  if (x == j_new){
                    x = std::numeric_limits<size_t>::max();
                    to_leave_cluster[k_new].push_back(get_ndata(k_new) - 1);
                    break;
                  }
                }
              }
            }
            
            //swap with tail.
            remove_with_tail_pull(k, j);

          } //j !=  max size_t
        }
      }
    }
    
    
////////////////////////////
////////// tests ///////////
////////////////////////////
    
    
    void post_initialisation_test(){
      mowri << " post initialisation test... " << zentas::Flush;
      mowri << " (injective) " << zentas::Flush;
      injective_ID_test();
      mowri << " (ndata) " << zentas::Flush;
      ndata_tests();
      mowri << " (info) " << zentas::Flush;
      info_tests();
      mowri << " (statistics) " << zentas::Flush;
      cluster_statistics_test();
      mowri << "initialised object looks ok." << zentas::Endl;
    }
    
    virtual void center_center_info_test(){}
    
    
    void post_center_update_test(){
      mowri << " post center update test... " << zentas::Flush;
      ndata_tests();
      center_center_info_test();
      mowri << "done." << zentas::Flush;
    }
    
    void post_sample_update_test(){
      mowri << " post sample update test... " << zentas::Flush;
      ndata_tests();
      info_tests();
      to_leave_cluster_test();
      mowri << "done." << zentas::Flush;
    }
    
    void post_redistribute_test(){
      mowri << " post redistribute test... " << zentas::Flush;
      ndata_tests();
      as_assigned_test();  
        mowri << "done." << zentas::Flush;

    }
    
    void post_update_statistics_test(){
      injective_ID_test();
      cluster_statistics_test();
    }
    
    void to_leave_cluster_test(){
      std::string errm ("error in leave_cluster_test");
      for (size_t k = 0; k < K; ++k){
        for (size_t j = 0; j < get_ndata(k); ++j){
          if (nearest_1_infos[k][j].a_x != k){
            bool is_listed = false;
            for (auto & x : to_leave_cluster[k]){
              if (x == j){
                is_listed = true;
              }
            }
            if (is_listed == false){
              mowri << "\nis_listed == false" << zentas::Endl;
              throw std::logic_error(errm);
            }
          }
        }
        
        
        for (auto & x : to_leave_cluster[k]){
          
          size_t nappears = 0;
          for (auto & xp : to_leave_cluster[k]){
            if (xp == x){
              ++nappears;
            }
          }
          if (nappears != 1){
            mowri << std::to_string(x) << " appears " << nappears << " times in the to leave list of cluster " << k << zentas::Endl;
          }
          
          if (nearest_1_infos[k][x].a_x == k){
            mowri << "\ncluster and size ( " << k << " : " << get_ndata(k) << ") " << zentas::Endl;
            mowri << "j : " << x << zentas::Endl;
            mowri << "size of to_leave_cluster[k] " << to_leave_cluster[k].size() << zentas::Endl;
            mowri << "members of to_leave_cluster[k]" << zentas::Endl;
            for (auto & x : to_leave_cluster[k]){
              mowri << x << " " << zentas::Flush;
            }
            mowri <<  "\n(a,d,e) : " << nearest_1_infos[k][x].get_string() << zentas::Endl;
            mowri << "\nto leave cluster but k is 'k'" << zentas::Endl;
            throw std::logic_error(errm);
          }
        }
      }
    }
  
    virtual void custom_cluster_statistics_test() = 0;

    void cluster_statistics_test(){
      std::string errs("cluster statistics test failed.");
      //energies.
      double E__0 = 0;
      
      for (size_t k = 0; k < K; ++k){
        double E__0k = 0;
        for (size_t j = 0; j < get_ndata(k); ++j){
          E__0k += nearest_1_infos[k][j].e_x;
        }
        if (E__0k != cluster_energies[k]){
          throw std::logic_error(errs);
        }
        E__0 += E__0k;
      }
      
      if (E__0 != get_E_total()){
        throw std::logic_error(errs);
      }
      
      custom_cluster_statistics_test();
      
    }
  
  
    virtual void custom_info_test() {}
    
    void info_tests(){
      std::string errm("info test failed. ");
      for (size_t k = 0; k < K; ++k){
        for (size_t j = 0; j < get_ndata(k); ++j){
          std::vector<double> distances(K);
          size_t k_first_nearest;
          double d_first_nearest = std::numeric_limits<double>::max();
          for (size_t k2 = 0; k2 < K; ++k2){
            set_center_sample_distance(k2, k, j, distances[k2]);
            if (distances[k2] < d_first_nearest){
              d_first_nearest = distances[k2];
              k_first_nearest = k2;
            }
          }
          double e_first_nearest = f_energy(d_first_nearest);
          
          //if (k_first_nearest != nearest_1_infos[k][j].a_x){
            //mowri << "\n" << k_first_nearest << "  " << nearest_1_infos[k][j].a_x << zentas::Endl;
            //mowri << d_first_nearest << "  " << nearest_1_infos[k][j].d_x << zentas::Endl;
            //throw std::logic_error(errm + " k_first_nearest != nearest_1_infos[k][j].a_x");
          //}
  
          if (d_first_nearest != nearest_1_infos[k][j].d_x){
            mowri << "k: " << k << "  j: " << j << zentas::Endl;
            mowri << "\n" << "k1 (comp) " << k_first_nearest << "    k1 (testing) " << nearest_1_infos[k][j].a_x << zentas::Endl;
            mowri << std::setprecision(20);
            mowri <<  "d1 (comp) " << d_first_nearest << "    d1 (testing) " << nearest_1_infos[k][j].d_x << zentas::Endl;
            
            mowri << "the first min(10, size) distances are : " << zentas::Endl;
            for (size_t jp = 0; jp < std::min<size_t>(10, nearest_1_infos[k].size()); ++jp){
              mowri << " " << nearest_1_infos[k][j].d_x << " ";
            }
            mowri << zentas::Endl;
            
            throw std::logic_error(errm + "d_first_nearest != nearest_1_infos[k][j].d_x");
          }
  
          if (e_first_nearest != nearest_1_infos[k][j].e_x){
            throw std::logic_error(errm + "e_first_nearest != nearest_1_infos[k][j].e_x");
          }          
        }      
      }
      custom_info_test();
    }
  
    virtual void custom_ndata_test(){}
    void ndata_tests(){
      std::string errm("ndata test failure");
      size_t ndata_internal = 0;
      for (size_t k = 0; k < K; ++k){
        if (get_ndata(k) != nearest_1_infos[k].size()){
          throw std::logic_error(errm);
        }
        ndata_internal += get_ndata(k);
      }
      
      custom_ndata_test();
      
      
      if (ndata_internal != ndata - K){
        throw std::logic_error(errm);
      }
    }

    void as_assigned_test(){
      for (size_t k = 0; k < K; ++k){
        for (size_t j = 0; j < get_ndata(k); ++j){
          if (nearest_1_infos[k][j].a_x != k){
            std::string errstring = "A sample in cluster " + std::to_string(k) + " has a_x " + std::to_string(nearest_1_infos[k][j].a_x);
            throw std::runtime_error(errstring);
          }
        }
      }
    }
    

    void injective_ID_test(){
      for (size_t i = 0; i < ndata; ++i){
        size_t count_i = 0;
        std::vector<size_t> as_center_ID;
        std::vector<size_t> in_cluster;
        std::vector<size_t> at_index;
        
        for (size_t k = 0; k < K; ++k){
          if (i == center_IDs[k]){
            ++count_i;
            as_center_ID.push_back(k);
          }
          for (size_t j = 0; j < get_ndata(k); ++j){
            if (i == sample_IDs[k][j]){
              ++count_i;
              in_cluster.push_back(k);
              at_index.push_back(j);
            }
          }
        }
        if (count_i != 1){
          std::string errstring = "Error caught. ID " + std::to_string(i) + " appears not exactly once. ";
          errstring = errstring + "It appears as centers of [";
          for (auto & x : as_center_ID){
            errstring = errstring + " " + std::to_string(x) + " ";
          }
          errstring = errstring + " ], as well as as member (k,j {size of k}) in [";
          for (size_t q = 0; q < in_cluster.size(); ++ q){
            errstring = errstring + " (" + std::to_string(in_cluster[q]) + "," + std::to_string(at_index[q]) + " {" + std::to_string(get_ndata(in_cluster[q]))+ "}) ";
          }
          
          errstring = errstring + "].";
          throw std::logic_error(errstring);
        }
      }
    }


    void post_initialise_centers_test(){
      
      for (size_t k = 0; k < K; ++k){
        if (center_indices_init[k] == center_indices_init[(k+1)%ndata]){
          std::stringstream errm_ss;
          errm_ss << "initialising indices must be unique, received " << center_indices_init[k] << " at least twice. The initialising indices are \n";
          for (size_t k = 0; k < K; ++k){
            errm_ss << " " << center_indices_init[k] << " ";
          }
          errm_ss << "\n";
          throw std::runtime_error(errm_ss.str());
        }
      }
    }


};






} //namespace nszen
 

#endif
