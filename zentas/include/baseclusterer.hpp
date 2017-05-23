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

#include <algorithm>
#include <type_traits>
#include <chrono>
#include <random>


#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>



namespace nszen{
 
 
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
  
    zentas::outputwriting::OutputWriter mowri; //(true, false, "");

    typedef typename TMetric::Initializer TMetricInitializer;
    typedef typename TData::DataIn DataIn;
    const size_t K;
    const size_t ndata;
    //const TEnergyFunc f_E_rock;
    std::function<double(double)> f_energy; //makes no difference if this or that.

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
    std::chrono::time_point<std::chrono::high_resolution_clock> tstart;
    

    size_t ncalcs_in_update_centers = 0;
    size_t ncalcs_in_update_sample_info = 0;
    size_t ncalcs_initialising = 0;
    size_t ncalcs_total = 0;


  protected:
    std::uniform_int_distribution<size_t> dis;
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
    K(K), ndata(datain.get_ndata()), //f_E_rock(),
    centers_data(datain, true), nearest_1_infos(K), sample_IDs(K), to_leave_cluster(K), cluster_has_changed(K, true), ptr_datain(& datain),
    metric(datain, nthreads, metric_initializer), 
    cluster_energies(K,0), cluster_mean_energies(K), E_total(std::numeric_limits<double>::max()), old_E_total(0),
    round(0), v_center_indices_init(K), center_indices_init(v_center_indices_init.data()), gen(seed), maxtime_micros(static_cast<size_t>(maxtime*1000000.)), minmE(minmE), labels(labels), nthreads(nthreads), nthreads_fl(static_cast<double> (nthreads)), maxrounds(maxrounds), energy(energy)

     {
       
       
      

       
       
       
      if (energy.compare("identity") == 0){
        f_energy = nszen::Identity(); //std::function<double(double)> ( [](double x){ return x*x; } );
      } 
       
      else if (energy.compare("quadratic") == 0){
        f_energy = nszen::Quadratic(); //std::function<double(double)> ( [](double x){ return x*x; } );
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
        for (size_t k = 0; k < K; ++k){
          /* confirm uniqueness and range. TODO : do I want to do this? performance? */
          if (center_indices_init_[k] >= ndata){
            throw std::runtime_error("initialising center index out of bounds in BaseClusterer constructor");
          }
          for (unsigned j = 0; j < k; ++j){
            if (center_indices_init_[k] == center_indices_init_[j]){
              throw std::runtime_error("initialising center index at j and k are the same in BaseClusterer constructor");
            }
          }
          
          v_center_indices_init[k] = center_indices_init_[k];
        }
      }
      
      
      /* initialisation from indices. */
      else if (initialisation_method == "uniform"){

        bool accepted;
        size_t proposed_i;
        for (size_t k = 0; k < K; ++k){
          accepted = false;
          while (accepted == false){
            accepted = true;
            proposed_i = dis(gen)%ndata;
            
            for (size_t k_m = 0; k_m < k; ++k_m){
              if (v_center_indices_init[k_m] == proposed_i){
                accepted = false;
              }
            }
          }
          if (accepted == true){
            v_center_indices_init[k] = proposed_i;
          }
        }
      }
      
      else if (initialisation_method.substr(0,8) == "afk-mc2-"){
        
        if (initialisation_method.size() < 9){
          std::stringstream errm_ss;
          errm_ss << "invalid initialisation_method " << initialisation_method << ". It is not of the form afk-mc2-INT, where INT is positive."; 
          throw std::runtime_error(errm_ss.str());
        }
        
        std::string digit_substring = initialisation_method.substr(8, initialisation_method.size() - 8);
        auto striter = digit_substring.begin(); 
        while (striter != digit_substring.end()){
          char x = *striter;
          if (std::isdigit(x) == false){
            std::stringstream errm_ss;
            errm_ss << "Unexpected character while attempting to extract integer from " << initialisation_method << " (" << digit_substring << ")" << ", `" << x << "'";
            throw std::runtime_error(errm_ss.str());
          }
          ++striter;
        }
        
        size_t chain_length = std::stoi(digit_substring);
        mowri << "chain length : " << chain_length << zentas::Endl;
        
        
        
        
        
        
        
        /* number of attempted moves before sampling */
        
        

        double adistance;
        
        
        std::vector<bool> is_center (ndata, false);
        /* the first center is chosen at random */
        size_t index0 = dis(gen)%ndata;
        v_center_indices_init[0] = index0;
        is_center[index0] = true;

        /* energies */
        std::vector<double> e0 (ndata);
        /* total unnormalised energy */
        double e0_sum = 0;        
        /* cumulative normalised energies */
        std::vector<double> E0 (ndata);        
        
        
        /* set raw energies and Z (E0_sum) */
        for (size_t i = 0; i < ndata; ++i){
          metric.set_distance(datain.at_for_metric(index0), datain.at_for_metric(i), adistance);
          e0[i] = f_energy(adistance);
          e0_sum += e0[i];
        }
        
        E0[0] = e0[0];
        for (size_t i = 1; i < ndata; ++i){
          E0[i]  = E0[i-1] + e0[i];
        }
        for (size_t i = 0; i < ndata; ++i){
          E0[i] /= e0_sum;
        }
        
        /* will now sample 2*ndata samples from (almost) q of Bachem et al. */ 
        size_t n_pseudo_samples = 2*ndata;
        std::vector<size_t> pseudo_sample (n_pseudo_samples);

        double d_ndata = static_cast<double>(ndata);        
        std::uniform_real_distribution<double> dis_uni01(0.0, 1.0);
        double rand_offset = dis_uni01(gen) / d_ndata;
        
        unsigned running_index = 0;
        for (size_t i = 0; i < ndata; ++i){
          /* the uniform distribution component (with replacement, not what Bachem et al do but, but good enough) */
          pseudo_sample[2*i] = i;
          /* the non-uniform distribution component. Again, the sampling is not independent as in Bachem et al, but good enough  */
          while (E0[running_index] < (static_cast<double>(i) + rand_offset)/d_ndata){
            ++ running_index;
          }
          pseudo_sample[2*i+1] = running_index; 
        }
        
        /* shuffle the samples */
        size_t swap_i;
        for (size_t i = 0; i < ndata; ++i){
          swap_i = dis(gen)%(ndata - i);
          std::swap(pseudo_sample[swap_i], pseudo_sample[ndata - 1 - i]);
        }
        
        //std::cout << "pseudo samples : \n";
        //for (size_t i = 0; i < 2*ndata; ++i){
          //std::cout << pseudo_sample[i] << " " << std::flush;
        //}

        size_t q_index = 0;
        /* x of Bachem et al */
        size_t sample_index_current;
        /* d_x of Bachem et al */
        double e_current;

        /* y of Bachem et al */
        size_t sample_index_proposal;
        /* d_y of Bachem et al */
        double e_proposal;
        


        
        for (size_t k = 1; k < K; ++k){
          
          do {
            q_index += 1;
            q_index %= n_pseudo_samples;
            sample_index_current = pseudo_sample[q_index];            
          } while (is_center[sample_index_current] == true);
          
          
          e_current = std::numeric_limits<double>::max();
          for (size_t kp = 0; kp < k; ++kp){
            metric.set_distance(datain.at_for_metric(sample_index_current), 
                                datain.at_for_metric(v_center_indices_init[kp]), 
                                adistance);
            e_current = std::min<double>(e_current, f_energy(adistance));
          }
          
          for (size_t m = 1; m < chain_length; ++m){
            
            do {
              q_index += 1;
              q_index %= n_pseudo_samples;
              sample_index_proposal = pseudo_sample[q_index];
            } while (is_center[sample_index_proposal] == true);
            

            e_proposal = std::numeric_limits<double>::max();
            for (size_t kp = 0; kp < k; ++kp){
              metric.set_distance(datain.at_for_metric(sample_index_proposal), 
                                  datain.at_for_metric(v_center_indices_init[kp]), 
                                  adistance);
              e_proposal = std::min<double>(e_proposal, f_energy(adistance));
            }
  
              
            if ((e_proposal/e_current)*((e0[sample_index_current]*2*ndata + e0_sum)/(e0[sample_index_proposal]*2*ndata + e0_sum))  >  dis_uni01(gen)){
              e_current = e_proposal;
              sample_index_current = sample_index_proposal;
            }
          }
          is_center[sample_index_current] = true;
          v_center_indices_init[k] = sample_index_current;          
        }
      }
      
      else {
        std::stringstream vims_ss;
        vims_ss << "The valid strings for initialisation_method are [from_indices_init, uniform, afk-mc2-INT (for some positive INT)]";
        throw std::runtime_error(vims_ss.str());
      }

      std::sort(v_center_indices_init.begin(), v_center_indices_init.end());

      center_IDs = indices_final;
      std::copy(center_indices_init, center_indices_init + K, center_IDs);

      for (size_t k = 0; k < K; ++k){
        centers_data.append(datain.at_for_move(center_indices_init[k]));
        // here we make an empty cluster, using datain to determine nec. `shape' parameters
        cluster_datas.emplace_back(datain, true);
      }
      
      for (size_t k = 0; k < K; ++k){
        if (center_indices_init[k] == center_indices_init[(k+1)%ndata]){
          std::stringstream errm_ss;
          errm_ss << "initialising indices must be unique, received " << center_indices_init[k] << " at least twice\n";
          throw std::runtime_error(errm_ss.str());
        }
      }
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
        //mowri << maxtime_micros - time_total << zentas::Endl;
        return maxtime_micros - time_total;
      }
      else{
        return -1;
      }
    }


    //TODO : move to private section
    inline void final_push_into_cluster(size_t i, size_t nearest_center, double min_distance, const double * const distances){
      //get your lock on, time for polyphonics. 
      std::lock_guard<std::mutex> lockraii(mutex0);
      cluster_datas[nearest_center].append(ptr_datain->at_for_move(i));
      nearest_1_infos[nearest_center].emplace_back(nearest_center, min_distance, f_energy(min_distance));
      sample_IDs[nearest_center].push_back(i);
      /* as usual: the final parameter up_distances.get() guarantees correctness only of lowest and second lowest, other values may exceed */
      
      //TODO : there may be computation in here which should not be under the lock. 
      put_sample_custom_in_cluster(i, nearest_center, distances);
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
    
        
    std::string get_base_summary_string(){
      std::ostringstream out;
      out  << round << "\tmE: " <<  std::setprecision(7) << E_total / static_cast<double>(ndata) << std::setprecision(5) <<  "\titime: "  <<  time_initialising/1000  <<  "\tctime: "  <<  time_in_update_centers/1000  <<  "\tutime: "  <<  time_in_update_sample_info/1000  <<  "\trtime: "  <<  time_in_redistribute/1000  <<  "\tttime: "  <<  time_total/1000  <<   "\tlg2 nc(c): "  << std::log2(ncalcs_in_update_centers) <<  "\tlg2 nc: "  <<  std::log2(ncalcs_total)  <<  "\t pc: " << static_cast<double> (metric.get_rel_calccosts())  ;
      
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
    
    
    void go(){
      
      bool with_tests = false;
      if (with_tests == true){
        mowri << "\n\nCOMPILED WITH TESTS ENABLED : WILL BE SLOW" <<zentas::Endl;
      }
      
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

      //mowri << "round :" << -1 << "\t tenergy : " <<  E_total << " \t itime : " << time_initialising << " \t ctime : " << time_in_update_centers << " \t utime : " << time_in_update_sample_info << " \t rtime : " << time_in_redistribute << " \t ttime : " << time_total << zentas::Endl;        
  
      t1 = std::chrono::high_resolution_clock::now();
      time_initialising = std::chrono::duration_cast<std::chrono::microseconds>(t1 - tstart).count();
      ncalcs_initialising = metric.get_ncalcs();
      time_total = time_initialising;
      
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
      
      
      //pll_put_samples_in_cluster(0, ndata - K, non_center_IDs);
      
      std::vector<std::thread> threads;      
      for (size_t ti = 0; ti < nthreads; ++ti){
        threads.emplace_back([this,ti,  &non_center_IDs] () { pll_put_samples_in_cluster(get_start(ti, get_nthreads(), 0, ndata - K), get_end(ti, get_nthreads(), 0, ndata - K), non_center_IDs); } );
      }
      for (auto & t : threads){
        t.join();
      }

    }
      
    void pll_put_samples_in_cluster(size_t i_a, size_t i_z, std::vector<size_t> & non_center_IDs){

      for (size_t i = i_a; i < i_z; ++i){
        put_sample_in_cluster(non_center_IDs[i]);
      }
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
  
    virtual void custom_cluster_statistics_test() = 0;//{}

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
            mowri << "\n" << k_first_nearest << "  " << nearest_1_infos[k][j].a_x << zentas::Endl;
            mowri << std::setprecision(20) <<  d_first_nearest << "  " << nearest_1_infos[k][j].d_x << zentas::Endl;
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
};




} //namespace nszen
 

#endif
