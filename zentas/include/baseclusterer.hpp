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

#include "skeletonclusterer.hpp"


#include "tdata.hpp"
#include "tdatain.hpp"
#include "tmetric.hpp"




#include "zentasinfo.hpp"







namespace nszen{

/* Another k-means++ helper class */
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






template <class TDataIn, class TMetric>
struct BaseClustererInitBundle{

  typedef typename TMetric::Initializer TMetricInitializer;

  size_t K;
  const TDataIn * const ptr_datain;
  const size_t * const center_indices_init;
  std::string initialisation_method;
  size_t seed;
  double max_time;
  double min_mE;
  size_t * const indices_final;
  size_t * const labels;
  size_t nthreads;
  size_t max_rounds;
  std::string energy;
  bool with_tests;
  const TMetricInitializer * ptr_metric_initializer;
  const EnergyInitialiser * ptr_energy_initialiser;
  std::chrono::time_point<std::chrono::high_resolution_clock> bigbang;

  
  BaseClustererInitBundle(size_t K, const TDataIn * const ptr_datain, const size_t * const center_indices_init, std::string initialisation_method_, size_t seed, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t max_rounds, std::string energy, bool with_tests,  const TMetricInitializer & metric_initializer, const EnergyInitialiser & energy_initialiser, const std::chrono::time_point<std::chrono::high_resolution_clock> & bigbang_): 
  
  K(K), ptr_datain(ptr_datain), center_indices_init(center_indices_init), initialisation_method(initialisation_method_), seed(seed), max_time(max_time), min_mE(min_mE), indices_final(indices_final), labels(labels), nthreads(nthreads), max_rounds(max_rounds), energy(energy), with_tests(with_tests), ptr_metric_initializer(&metric_initializer), ptr_energy_initialiser(&energy_initialiser), bigbang(bigbang_) {
    
  }

};





template <class TMetric, class TData>
class BaseClusterer : public SkeletonClusterer{


  
  /* A comment to put somewhere : when two clusters are equally near to a sample, the cluster to which it is assigned is unpredictable. One might have a statement that it goes to the one with lower index for example, but this would mean that if a center is at equal distance to the nearest it cannot be eliminated. This is less efficient, and well, I don't want to go back and change vast swathes of code. 
   * 
   * As a result of this unpredictability, l0, l1, l2 can differ: the important point is that random proposals depend on cluster sizes, and so while the energy and proposal evaluation do not care whether a sample is assigned to a cluster at equal distance to another, the subsequent proposals may differ. 
   * 
   * */
  
  public:

    typedef typename TMetric::Initializer TMetricInitializer;
    typedef typename TData::DataIn DataIn;
    typedef typename TData::RefinementCenterData RefinementCenterData;

  private:

    //CommonBase cb;
    TData centers_data;
    RefinementCenterData rf_data;
    std::vector<TData> cluster_datas; 
    const DataIn * const ptr_datain;
    TMetric metric;
    std::unique_ptr<KmooBundle<TData>> up_kmoo_bundle;
  
  public: 


/////////////////////////////////////
/////////// public //////////////////
/////////////////////////////////////
  public:
    BaseClusterer(size_t K, /* number of clusters */
     const DataIn & datain, /* initialisating data */
     const size_t * const center_indices_init_, /* The K sample indices to initialise with (if initialisation_method is from_indices_init)*/
     std::string init_method, /* */
     size_t random_sd, /* starting seed for generating random numbers, uses a custom "std::default_random_engine"  */
     double max_time, /* will stop in go ( ) at first opportunity after max_time */
     double min_mE, /* will stop in go ( ) at first opportunity if mE is less than min_mE */
     size_t * const indices_final, /* the K final indices (to populate) */
     size_t * const labels, /* the assigned cluster of the datain.ndata samples */
     size_t nthreads,
     size_t max_rounds, 
     std::string energy, 
     bool with_tests,
     const TMetricInitializer & metric_initializer, 
     const EnergyInitialiser & energy_initialiser,
     std::chrono::time_point<std::chrono::high_resolution_clock> bigbang): 
    
    
    SkeletonClusterer(K, bigbang, center_indices_init_, datain.get_ndata(), init_method, max_time, min_mE, indices_final, labels, nthreads, max_rounds, energy, with_tests, random_sd, energy_initialiser),
    
    centers_data(datain, true), rf_data(datain, true), ptr_datain(& datain), metric(datain, nthreads, metric_initializer), up_kmoo_bundle(new KmooBundle<TData>(ptr_datain)) {}

    BaseClusterer(const BaseClustererInitBundle<DataIn, TMetric> & ib): BaseClusterer(ib.K, *(ib.ptr_datain), ib.center_indices_init, ib.initialisation_method, ib.seed, ib.max_time, ib.min_mE, ib.indices_final, ib.labels, ib.nthreads, ib.max_rounds, ib.energy, ib.with_tests, *(ib.ptr_metric_initializer), *(ib.ptr_energy_initialiser), ib.bigbang) {}

    
    std::string string_for_sample(size_t k, size_t j) {
      return cluster_datas[k].string_for_sample(j);
    }
    
    std::string string_for_center(size_t k){
      return centers_data.string_for_sample(k);
    }


     void final_push_into_cluster_basic(size_t i, size_t nearest_center, double min_distance){
      cluster_datas[nearest_center].append(ptr_datain->at_for_move(i));
      nearest_1_infos[nearest_center].emplace_back(nearest_center, min_distance, f_energy(min_distance));
      sample_IDs[nearest_center].push_back(i);      
    }

    

    void reset_multiple_sample_infos(size_t k_to, size_t j_a, size_t j_z){
      for (size_t j = j_a; j < j_z; ++j){
        reset_sample_infos(k_to, j);
      }
    }
    
    



    
    void base_put_sample_in_cluster(size_t i) {
      
      
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
    
    void triangular_put_sample_in_cluster(size_t i, const double * const cc) {
      std::unique_ptr<double []> up_distances (new double [K]);
      double min_distance = std::numeric_limits<double>::max();
      double second_min_distance = std::numeric_limits<double>::max();

      size_t nearest_center = 0;
      
      //double xi = (1 - 1e-6);
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
        throw zentas::zentas_error(errm.str());
      }
 
      if (ndata_1 != ndata_2){
        throw zentas::zentas_error("ndata_1 != ndata_2 in test_parameters_to_tkai");        
      }
    }




     void kmpp_inner(size_t i, size_t k, double a_distance, P2Bundle & p2bun){
    
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
     * is_exhausted reports whether the every sample has a center a distance 0 (when samples duplicated). 
     *  */

    void triangular_kmeanspp_after_initial(TData & c_dt, size_t * const c_ind, double * const cc, P2Bundle & p2bun, const TData * const ptr_p2bun_dt, size_t k0, size_t k1, bool from_full_data, bool & is_exhausted){
    
      if (from_full_data != (ptr_p2bun_dt == nullptr)){
        throw zentas::zentas_error("logic error in triangular_kmeanspp_after_initial : from_full_data != (ptr_p2bun_dt == nullptr)");
      }
      
      size_t kmeanspp_ndata = from_full_data ? ndata : ptr_p2bun_dt->get_ndata();
      
      if (kmeanspp_ndata -1 < k1 - k0){
        std::stringstream ss;
        ss << "triangular_kmeanspp_after_initial, attempting to find more centers than data - 1 : ";
        ss << "kmeanspp_ndata = " << kmeanspp_ndata << "   and   k1 - k0 = " << k1 - k0; 
        throw zentas::zentas_error(ss.str());        
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
      
      /* whether to revert to uniform sampling when all samples already have a center at distance 0 (obviously duplicated data)*/
      bool try_uniform_when_exhausted;

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
        
        /* we can handle the exhausted case here */
        try_uniform_when_exhausted = true;
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

        /* handling the exhausted case here is too complex (cross batch indices needed, argh)
         * , rather just bail and try again from full data */        
        try_uniform_when_exhausted = false;
      
      }
      
      
      
      /* k-means++, at last */
      for (size_t k = k0; k < k1; ++k){
        
        if (v_cum_nearest_energies.back() != 0){
          is_exhausted = false;
          update_c_ind_c_dt(k);
        }
        
        else{
          is_exhausted = true;
          if (try_uniform_when_exhausted == false){
            /* bailing, because 
             * (1) every sample has a center a distance 0 : k-means++ fails
             * (2) must not switch to uniform sampling.
             * */
            return;
          }
          
          else{
            if (from_full_data == false){
              std::stringstream ss;
              ss << "try uniform when exhausted is true, but from full data is false, this looks like a logic error."; 
              throw zentas::zentas_error(ss.str());
            }

            size_t new_index = dis(gen)%(K);
            while (std::find(c_ind, c_ind+k, new_index) != c_ind+k){
              new_index = dis(gen)%(K);
            }
            c_ind[k] = new_index;
            c_dt.append(ptr_datain->at_for_move(c_ind[k]));
          }
        }
 
        
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


    

    void triangular_kmeanspp(size_t * const centers, double * const cc, P2Bundle & p2bun, TData & c_dt, bool & is_exhausted){


      for (size_t i = 0; i < ndata; ++i){
        p2bun.ori(i) = i;
      }
      
      size_t k0 = 0;
      size_t k1 = K;
    
      triangular_kmeanspp_after_initial(c_dt, centers, cc, p2bun, nullptr, k0, k1, true, is_exhausted);
    }




    
        
    void triangular_kmeanspp_aq2(size_t * const centers, double * const cc, P2Bundle & p2bun, TData & c_dt, size_t n_bins, bool & is_exhausted){

      
      /* experiments so far show that multithreading does not help here, can hurt. what's weird is that even if nthreads = 1 in 
       * the pll version, it's sig slower than the serial version.  */
      bool multithread_kmpp = true;
      
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
        if (nthreads == 1 || multithread_kmpp == false){
          
          //auto ndat = p2buns[bin].get_ndata();
          //std::thread t([bin, k0, ndat, &update_nearest_info](){
            //update_nearest_info(bin, 0, k0, 0, ndat);
          //});
          //t.join();
          
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
        triangular_kmeanspp_after_initial(c_dt, centers, cc, p2buns[bin], &p2buns_dt[bin], k0, k1, false, is_exhausted);
        /* is_exhausted is true, we will bail as handling this here (with n_bins > 1) is too complex. Will try again with n_bins = 1 */
        if (is_exhausted){
          return;
        }
        
      }
      
      /* update all until the tail k's (*&^*) */
      for (size_t bin = 0; bin < n_bins; ++bin){
        size_t k1 = ((bin + 1)*non_tail_k)/n_bins;

        
        if (nthreads == 1 || multithread_kmpp == false){
          //auto ndat = p2buns[bin].get_ndata();
          //std::thread t([bin, k1, non_tail_k, ndat, &update_nearest_info](){
            //update_nearest_info(bin, k1, non_tail_k, 0, ndat);
          //});
          //t.join();
          
          
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
      triangular_kmeanspp_after_initial(c_dt, centers, cc, p2bun, nullptr, K - tail_k, K, true, is_exhausted);
      /* is_exhausted is true, we will bail as handling this here (with n_bins > 1) is too complex. Will try again with n_bins = 1 */
      if (is_exhausted){
        return;
      }

    }
    

    virtual double get_rel_calccosts() override final{
      return metric.get_rel_calccosts();
    }
    
    virtual size_t get_ncalcs() override final{
      return metric.get_ncalcs();
    }
    
      
        
     size_t get_a1(size_t k, size_t j) const{
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
    
     size_t get_ndata(size_t k){
      return cluster_datas[k].get_ndata();
    }
    
     size_t get_nthreads(){
      return nthreads;
    }
    
     double get_nthreads_fl(){
      return nthreads_fl;
    }
    
    
    
    
    
    
     void swap_center_with_sample(size_t k, size_t j){
      cluster_has_changed[k] = true;
      
      size_t c_id_k = center_IDs[k];
      center_IDs[k] = sample_IDs[k][j];
      /* jn : bug fix, voronoi 11 june 2017. the following was missing: */
      sample_IDs[k][j] = c_id_k; 
      
      nszen::swap<TData>(centers_data, k, cluster_datas[k], j);
      
      
      reset_sample_infos_basic(k, j);
    }
          
    
     double get_E_total(){
      return E_total;
    }
    
     double get_cluster_energy(size_t k){
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
     void remove_with_tail_pull(size_t k, size_t j){
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
        throw zentas::zentas_error("(see above)");
      }
      size_t i = dis(gen)%(ndata - K);
      size_t k_below = 0;
      while (i >= up_cum_ndatas[k_below]){
        ++k_below;
      }
 
      return k_below;
    }
    
     void signal_cluster_change(size_t k){
      cluster_has_changed[k] = true;
    }
    
    size_t draw_j_uniform(size_t k){
      return dis(gen)%get_ndata(k);
    }
    
    //set the distance from center k to the j1'th element of cluster k1.
     void set_center_sample_distance(size_t k, size_t k1, size_t j1, double & distance) {
      metric.set_distance(centers_data.at_for_metric(k), cluster_datas[k1].at_for_metric(j1), distance);
    }

     void set_center_sample_distance(size_t k, size_t k1, size_t j1, double threshold, double & distance) {
      metric.set_distance(centers_data.at_for_metric(k), cluster_datas[k1].at_for_metric(j1), threshold, distance);
    }
        
     void set_center_sampleID_distance(size_t k, size_t i, double threshold, double & distance) {
      metric.set_distance(centers_data.at_for_metric(k), ptr_datain->at_for_metric(i), threshold, distance);
    }

     void set_sampleID_sampleID_distance(size_t i1, size_t i2, double threshold, double & distance) {
      metric.set_distance(ptr_datain->at_for_metric(i1), ptr_datain->at_for_metric(i2), threshold, distance);
    }

     void set_sampleID_sampleID_distance(size_t i1, size_t i2, double & distance) {
      metric.set_distance(ptr_datain->at_for_metric(i1), ptr_datain->at_for_metric(i2), distance);
    }
        
        
     double get_center_sample_distance(size_t k, size_t k1, size_t j1) {
      double distance;
      set_center_sample_distance(k, k1, j1, distance);
      return distance;
    }
    
     double get_sample_sample_distance(size_t k , size_t j1, size_t j2) {
      double adistance;
      metric.set_distance(cluster_datas[k].at_for_metric(j1), cluster_datas[k].at_for_metric(j2), adistance);
      return adistance;
    }
    
     void set_sample_sample_distance(size_t k1, size_t j1, size_t k2, size_t j2, double threshold, double & adistance) {
      metric.set_distance(cluster_datas[k1].at_for_metric(j1), cluster_datas[k2].at_for_metric(j2), threshold, adistance);
    }
    
    
     void set_center_center_distance(size_t k1, size_t k2, double threshold, double & adistance) {
      metric.set_distance(centers_data.at_for_metric(k1), centers_data.at_for_metric(k2), threshold, adistance);
    }
     void set_center_center_distance(size_t k1, size_t k2, double & adistance) {
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
          throw zentas::zentas_error("n_bins passed to kmeans++ should be positive.");
        }
      }

      
      up_kmoo_bundle->reset(K, ndata);
      
      bool is_exhausted = false;

      if (n_bins != 1){
        triangular_kmeanspp_aq2(center_indices_init, up_kmoo_bundle->cc.data(), up_kmoo_bundle->p2bun, up_kmoo_bundle->c_dt, n_bins, is_exhausted);
        if (is_exhausted){
          mowri << "exhausted in " << initialisation_method << ", will revert to kmeans++-1" << zentas::Endl;
          up_kmoo_bundle->reset(K, ndata);
        }
      }

      
      if (is_exhausted == true || n_bins == 1){
        is_exhausted = false;
        triangular_kmeanspp(center_indices_init, up_kmoo_bundle->cc.data(), up_kmoo_bundle->p2bun, up_kmoo_bundle->c_dt, is_exhausted);
        if (is_exhausted){
          mowri << "exhausted in kmeans++-1, used uniform sampling to complete initialisation" << zentas::Endl;
        }
      }
      

      
      
      /* rearrange so that center_indices_init are in order */
      std::vector<std::array<size_t, 2>> vi (K);
      for (size_t k = 0; k < K; ++k){
        std::get<0>(vi[k]) = center_indices_init[k];
        std::get<1>(vi[k]) = k;
      }
      
      auto fg = std::less<size_t>();
      std::sort(vi.begin(), vi.end(), [&fg](std::array<size_t, 2> & lhs, std::array<size_t, 2> & rhs){ return fg(std::get<0>(lhs), std::get<0>(rhs)); });
      
      std::vector<size_t> old_to_new (K);
      for (unsigned k = 0; k < K; ++k){
        center_indices_init[k] = std::get<0>(vi[k]);
        old_to_new[std::get<1>(vi[k])] = k;
      }
      
      for (unsigned i = 0; i < ndata; ++ i){
        up_kmoo_bundle->p2bun.k_1(i) = old_to_new[up_kmoo_bundle->p2bun.k_1(i)];
        up_kmoo_bundle->p2bun.k_2(i) = old_to_new[up_kmoo_bundle->p2bun.k_2(i)];
      }
      
       
    }

    void initialise_center_indices(){
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
      
      /* initialisation with kmeans++ */
      else if (initialisation_method == "kmeans++" || initialisation_method.substr(0, 9) == "kmeans++-"){
        initialise_with_kmeanspp();
      }
      
      else {
        std::stringstream vims_ss;
        vims_ss << "The valid strings for initialisation_method are [from_indices_init, uniform, kmeans++-INT, afk-mc2-INT (for some positive INT)]";
        throw zentas::zentas_error(vims_ss.str());
      }
  
      std::sort(v_center_indices_init.begin(), v_center_indices_init.end());
      std::copy(center_indices_init, center_indices_init + K, center_IDs);

    }


    void output_halt_kmedoids_reason(){
      
      unsigned n_reasons = 0;
      mowri << "halted because: " << zentas::Endl;
      if (time_total >= max_time_micros){
        mowri << "  (" << n_reasons+1 << ") exceeded max_time (" << max_time_micros/1000. << ")" << zentas::Endl;
        ++n_reasons;
      }
      if (round >= max_rounds){
        mowri << "  (" << n_reasons+1 << ") exceeded max_rounds (" << max_rounds << ")" << zentas::Endl;
        ++n_reasons;
      }
      if (E_total / static_cast<double>(ndata) < min_mE){
        mowri << "  (" << n_reasons+1 << ") mE below termination min_mE (" << min_mE << ")" << zentas::Endl;        
        ++n_reasons;
      }
      if (n_reasons == 0){
        mowri << "   round without any center update" << zentas::Endl;
      }
      
    }

    bool halt_kmedoids(){
      bool do_not_halt = (time_total < max_time_micros) && (round < max_rounds) && (E_total / static_cast<double>(ndata)) >= min_mE;
      return (do_not_halt == false);
    };

    bool halt_refinement(){
      /* TODO. */
      return false;
    }

    void set_t_ncalcs_update_centers_start(){
      t_update_centers_start = std::chrono::high_resolution_clock::now();
      ncalcs_update_centers_start = metric.get_ncalcs();
    }


    void update_t_ncalcs_center_end(){
      t_update_centers_end = std::chrono::high_resolution_clock::now();
      time_in_update_centers += std::chrono::duration_cast<std::chrono::microseconds>(t_update_centers_end - t_update_centers_start).count();
      ncalcs_update_centers_end = metric.get_ncalcs();
      ncalcs_in_update_centers += ncalcs_update_centers_end - ncalcs_update_centers_start;
    }


    void update_t_ncalcs_sample_update_end(){
      t_update_sample_info_end = std::chrono::high_resolution_clock::now();
      time_in_update_sample_info += std::chrono::duration_cast<std::chrono::microseconds>(t_update_sample_info_end - t_update_centers_end).count();
      ncalcs_update_sample_info_end = metric.get_ncalcs();
      ncalcs_in_update_sample_info += ncalcs_update_sample_info_end - ncalcs_update_centers_end;
    }
 

    void update_t_redistibute_end(){ 
      t_redistribute_end = std::chrono::high_resolution_clock::now();
      time_in_redistribute += std::chrono::duration_cast<std::chrono::microseconds>(t_redistribute_end - t_update_sample_info_end).count();
    }

    void update_t_update_all_cluster_stats_end(){
      t_update_all_cluster_statistics_end = std::chrono::high_resolution_clock::now();
      time_in_update_all_cluster_statistics += std::chrono::duration_cast<std::chrono::microseconds>(t_update_all_cluster_statistics_end - t_redistribute_end).count();  
    }
  

    virtual bool refine_centers(){
      throw zentas::zentas_error("refine_centers not implemented");
    }
    
    void refine_center_center_info(){
      throw zentas::zentas_error("refine_center_center_info not implemented");
    }
    void refine_sample_info(){
      throw zentas::zentas_error("refine_sample_info not implemented");
    }
    

    
    void run_kmedoids(){
      if (in_refinement == true){
        throw zentas::zentas_error("in run_kmedoids, but in_refinement is true. This is not correct");
      }
      clustering_loop_scaffolding();
      mowri << get_equals_line(get_round_summary().size());
      output_halt_kmedoids_reason();      
    }
    
    void run_refinement(){
      if (in_refinement == false){
        throw zentas::zentas_error("in run_refinement, but in_refinement is false. This is not correct");
      }

      clustering_loop_scaffolding();
      mowri << get_equals_line(get_round_summary().size());
      mowri << "refinement complete, say something intelligent.";
    }
    
    void clustering_loop_scaffolding(){


      /* TODO : in refinement should be outside, if possible */
      auto bc_centers = [this](){if (in_refinement){ return refine_centers();}else{ return update_centers();}};
      auto bc_center_center_info = [this](){if (in_refinement){refine_center_center_info();}else{ update_center_center_info();}};
      auto bc_sample_info = [this](){if (in_refinement){refine_sample_info();}else{ update_sample_info();}};
      auto bc_halt =        [this](){if (in_refinement){return halt_refinement();}else{ return halt_kmedoids();}};

      while (bc_halt() == false){

        /* ************** *
        * UPDATE CENTERS *
        * ****************/
        set_t_ncalcs_update_centers_start();
        bool modified_centers = bc_centers();

        if (modified_centers == false){
          update_t_ncalcs_center_end();
          mowri << get_round_summary() << zentas::Endl;
          break;
        }

        /* ************************* *
        * UPDATE CENTER CENTER INFO *
        * ***************************/
        bc_center_center_info();
        if (with_tests == true) post_center_update_test();
        update_t_ncalcs_center_end();

         
        for (auto & x : to_leave_cluster){
          x.resize(0);
        }
        
        /* ****************** *
        * UPDATE SAMPLE INFO *
        * ********************/
        bc_sample_info();
        if (with_tests == true) post_sample_update_test();
        update_t_ncalcs_sample_update_end();

        
        /* ************ *
        * REDISTRIBUTE *
        * **************/
        redistribute();
        if (with_tests == true) post_redistribute_test();        
        update_t_redistibute_end();

        
        /* ************************* *
        * UPDATE CLUSTER STATISTICS *
        * ***************************/
        update_all_cluster_statistics();
        if (with_tests == true) post_update_statistics_test();
        update_t_update_all_cluster_stats_end();


        if (with_tests == true) mowri << zentas::Endl;
        mowri << get_round_summary() << zentas::Endl;
        ++round;
      }
    }
    
      
    void populate_labels(){
      for (size_t k = 0; k < K; ++k){
        labels[center_IDs[k]] = k;
        for (size_t j = 0; j < get_ndata(k); ++j){
          labels[sample_IDs[k][j]] = k;
        }
      }
    }
    
    void prepare_for_refinement();
        
    void go(){
      initialise_all();
      run_kmedoids();
      prepare_for_refinement();
      run_refinement();
      populate_labels();
    }

      
    void initialise_all() {
      if (with_tests == true){  
        mowri << "\n\nRUNNING WITH TESTS ENABLED : WILL BE SLOW" <<zentas::Endl;
      }

      mowri << 
R"((The prevent output to terminal, set capture_output to false)
(For a description of column statistics, consider function get_output_verbose_string())
)";      
      mowri << get_equals_line(get_round_summary().size());
      
      if (with_tests == true){
        fundamental_triangle_inequality_test();
      }



      tstart_initialise_centers = std::chrono::high_resolution_clock::now();       
      time_prehistory = std::chrono::duration_cast<std::chrono::microseconds>(tstart_initialise_centers - bigbang).count();
      
      initialise_center_indices();
      
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
      tstart = std::chrono::high_resolution_clock::now();

      //initialisation
      set_center_center_info();
      put_samples_in_clusters();                  
      set_all_cluster_statistics();

      /* checkpoint: all assignments and cluster statistics must be correct. Tests to pass : all */
      if (with_tests == true){
        post_initialisation_test();
      }

      auto t1 = std::chrono::high_resolution_clock::now();
      time_initialising = std::chrono::duration_cast<std::chrono::microseconds>(t1 - tstart).count();
      ncalcs_initialising = metric.get_ncalcs();
      //time_total = time_prehistory + time_initialising + time_to_initialise_centers;
      mowri << get_round_summary() << zentas::Endl;
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
        std::lock_guard<std::mutex> lock(mutex0);
        to_leave_cluster[k].push_back(j);
      }
    }
    


/////////////////////////////////////
/////////// private /////////////////
/////////////////////////////////////    
  private:

    virtual void put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) = 0;

    virtual void initialise_with_kmeanspp() = 0;
        
    /* called from put_sample_in_clusters (virtual functions are not nonsensical) */    
    virtual  void put_sample_custom_in_cluster(size_t, size_t, const double * const) {}

    
    /* when implementing this, one must guarantee that cluster_has_changed is correctly set in this function */
    virtual  void reset_sample_custom(size_t, size_t, size_t, const double * const) {}
    
    virtual bool update_centers() = 0;
         
    virtual void update_sample_info() = 0;
      
    /* functions called from within redistribute */
    virtual  void custom_append(size_t, size_t, size_t) {}
    
    virtual  void custom_replace_with_last(size_t, size_t) {}
    
    virtual  void custom_replace_with(size_t, size_t, size_t, size_t) {};
    
    virtual  void custom_remove_last(size_t) {};

    /* functions called from within set_cluster_statistics */
    virtual  void increment_custom_cluster_statistics(size_t, size_t) {};

    virtual void set_normalised_custom_cluster_statistics(size_t k) = 0;
    virtual void set_to_zero_custom_cluster_statistics(size_t k) = 0;

    virtual std::string get_round_summary() = 0;    

     void final_push_into_cluster(size_t i, size_t nearest_center, double min_distance, const double * const distances){
      //get your lock on, time for polyphonics. 
      std::lock_guard<std::mutex> lockraii(mutex0);
      
      /* the common part of pushing into a cluster */
      final_push_into_cluster_basic(i, nearest_center, min_distance);
      
      /* as usual: the final parameter up_distances.get() guarantees correctness
       * only of lowest and second lowest, other values may exceed */
      /* TODO : there may be computation in here which should not be under the lock */
      put_sample_custom_in_cluster(i, nearest_center, distances);

    }
    
    
     void final_push_into_cluster_post_kmeanspp(size_t i, size_t k1, size_t k2, double d1, double d2){
      std::lock_guard<std::mutex> lockraii(mutex0);
      final_push_into_cluster_basic(i, k1, d1);
      put_nearest_2_infos_margin_in_cluster_post_kmeanspp(k1, k2, d2, f_energy(d2));
    }


        
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
      // shuffle non center IDs      
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
      
      up_kmoo_bundle.reset();

    }
      
    void pll_put_samples_in_cluster(size_t i_a, size_t i_z, std::vector<size_t> & non_center_IDs){

      std::string prefix = "kmeans++";
      if (initialisation_method.substr(0, prefix.size()) == prefix){
        for (size_t i = i_a; i < i_z; ++i){

          size_t nc_ID = non_center_IDs[i];
          size_t k1 = up_kmoo_bundle->p2bun.k_1(nc_ID); 
          size_t k2 = up_kmoo_bundle->p2bun.k_2(nc_ID);
          double d1 = up_kmoo_bundle->p2bun.d_1(nc_ID);
          double d2 = up_kmoo_bundle->p2bun.d_2(nc_ID);
          final_push_into_cluster_post_kmeanspp(nc_ID, k1, k2, d1, d2);

        }
      }      

      else{        
        for (size_t i = i_a; i < i_z; ++i){
          put_sample_in_cluster(non_center_IDs[i]);
        }
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
        E_total += cluster_energies[k];
      }
    }
    

    
    virtual void put_sample_in_cluster(size_t i) = 0;
    
    
    /* redistribute_order is the order in which clusters `send away samples' */
    /* this function : set the vector to be indices {0,...K-1} in a specific order */
    /* clarans is best when k_to is first (done in virtual of baseclarans) */
    virtual void set_redistribute_order(std::vector<size_t> & redistribute_order) = 0;
    

    
    
    // single threaded redistribution. mulithreading here will be very tricky.
    // requires that to_leave_cluster be reliably set. 
    // centers are unchanged by this function. this function simply moves samples between clusters.
    // additional complexity required to maintain randomness (immigrant samples inserted into random indices)
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
      for (auto & k : redistribute_order){
        if (nthreads > 1){
          std::sort(to_leave_cluster[k].begin(), to_leave_cluster[k].end());
        }
        // thought : it would be nice if one could auto over a vector in reverse
        for (size_t ji = to_leave_cluster[k].size(); ji-- > 0;){
          
          j = to_leave_cluster[k][ji];
          
          /* size_t max indicates that the mover has moved to the tail, because of an earlier immigration */
          if (j != std::numeric_limits<size_t>::max()){
            k_new = nearest_1_infos[k][j].a_x;
            
            cluster_has_changed[k] = true;
            cluster_has_changed[k_new] = true;
            
            // now for the redistribution 
            // (1) make space somewhere in k_new (k,j) will be inserted at k_new, j_new:
            j_new = dis(gen)%(get_ndata(k_new) + 1);
            
            
            /* try to swap with a non-mover, as it makes swapping quicker. but if after 10 attempts no luck, 
             * can still swap with movers */
            size_t insert_attempts = 0;
            while (j_new != get_ndata(k_new) && nearest_1_infos[k_new][j_new].a_x != k_new && insert_attempts < 10){
              j_new = dis(gen)%(get_ndata(k_new) + 1);
              ++insert_attempts;
            }
                        
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
              
              /* putting k_new, j_new on the tail to make space for k, j */
              nearest_1_infos[k_new].push_back(nearest_1_infos[k_new][j_new]);
              sample_IDs[k_new].push_back(sample_IDs[k_new][j_new]);
              cluster_datas[k_new].append(cluster_datas[k_new].at_for_move(j_new));
              custom_append(k_new, k_new, j_new);
              
              /* putting k, j in where space has been made for it */
              nearest_1_infos[k_new][j_new] = nearest_1_infos[k][j];
              sample_IDs[k_new][j_new] = sample_IDs[k][j];
              cluster_datas[k_new].replace_with(j_new, cluster_datas[k].at_for_move(j));
              custom_replace_with(k_new, j_new, k,j);
              
              
              /* k_new, j_new is still going to move when k_new immigration starts. 
               * Record that k_new, j_new has moved to tail, and leave a notive (max size_t) 
               * so that when k_new immigration starts we know */
              if (target_is_a_mover && redistribute_rank[k_new] > redistribute_rank[k]){
                /* if to_leave_cluster is sorted we can improve this search. */
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
      mowri << "((post_initialisation_test):" << zentas::Flush;
      mowri << "(injective)" << zentas::Flush;
      injective_ID_test();
      mowri << "(ndata)" << zentas::Flush;
      ndata_tests();
      mowri << "(info)" << zentas::Flush;
      info_tests();
      mowri << "(statistics)" << zentas::Flush;
      mowri << ")" << zentas::Endl;
      cluster_statistics_test();
    }
    
    virtual void center_center_info_test(){}
    
    
    void post_center_update_test(){
      mowri << "((post_center_update_test):" << zentas::Flush;
      mowri << "(injective)" << zentas::Flush;
      injective_ID_test();
      mowri << "(ndata)" << zentas::Flush;
      ndata_tests();
      mowri << "(center_center_info)" << zentas::Flush;
      center_center_info_test();
      mowri << ")" << zentas::Endl;
    }
    
    void post_sample_update_test(){
      mowri << "(post_sample_update_test)" << zentas::Flush;
      mowri << "(ndata)" << zentas::Flush;
      ndata_tests();
      mowri << "(info)" << zentas::Flush;
      info_tests();
      mowri << "(to_leave_cluster)" << zentas::Flush;
      to_leave_cluster_test();
      mowri << "(injective)" << zentas::Flush;
      injective_ID_test();
      mowri << ")" << zentas::Endl;
    }
    
    void post_redistribute_test(){
      mowri << "(post_redistribute_test)" << zentas::Flush;
      mowri << "(ndata)" << zentas::Flush;
      ndata_tests();
      mowri << "(injective)" << zentas::Flush;
      injective_ID_test();
      mowri << "(as_assigned_test)" << zentas::Flush;
      as_assigned_test();  
      mowri << ")" << zentas::Endl;
    }
    
    void post_update_statistics_test(){
      mowri << "(injective)" << zentas::Flush;
      injective_ID_test();
      mowri << "(cluster_statistics)" << zentas::Flush;
      cluster_statistics_test();
      mowri << ")" << zentas::Endl;
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
              throw zentas::zentas_error(errm);
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
            throw zentas::zentas_error(errm);
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
          throw zentas::zentas_error(errs);
        }
        E__0 += E__0k;
      }
      
      if (E__0 != get_E_total()){
        throw zentas::zentas_error(errs);
      }
      
      custom_cluster_statistics_test();
      
    }
  
    void fundamental_triangle_inequality_test(){
      
      size_t n_tests = ndata*ndata;
      
      for (unsigned ti = 0; ti < n_tests; ++ti){
        size_t index0 = dis(gen)%ndata;
        size_t index1 = dis(gen)%ndata;
        size_t index2 = dis(gen)%ndata;
        double d01, d02, d12;
        
        set_sampleID_sampleID_distance(index0, index1, d01);
        set_sampleID_sampleID_distance(index0, index2, d02);
        set_sampleID_sampleID_distance(index1, index2, d12);
        
        if ((1 - 1e-5)*d02 > d01 + d12){
          std::stringstream ss;
          ss << "Triangle inequality failed :\n";
          ss << "Sample 0: " << ptr_datain->string_for_sample(index0) << "\n";
          ss << "Sample 1: " << ptr_datain->string_for_sample(index1) << "\n";
          ss << "Sample 2: " << ptr_datain->string_for_sample(index2) << "\n";
          ss << "d02 = " << d02 << ", d01 = " << d01 << ", d12 = " << d12 << ", d01 + d12 = " << (d01 + d12);
          throw zentas::zentas_error(ss.str());
        }
      }
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
          
  
          if (d_first_nearest != nearest_1_infos[k][j].d_x){

            std::stringstream errm;
            
            errm << "error detected : d_first_nearest != d_first_nearest\n";
            errm << "k=" << k << "\n";
            errm << "j=" << j << "\n";
            errm << "the " << j << "'th sample in cluster " << k << " is " << string_for_sample(k, j) << "\n";
            errm << "cluster size is " << nearest_1_infos[k].size() << "\n";
            errm << std::setprecision(20);
            errm << "get_a1(k,j)=" << get_a1(k,j) << "\n";
            errm << "just computed first nearest center index: " << k_first_nearest << "\n";
            errm << "the " << k_first_nearest << "'th center is: " << string_for_center(k_first_nearest) << "\n";
            errm <<  "just computed distance to this center is: " << d_first_nearest << "\n";
            errm << "the recorded first nearest center index: " << nearest_1_infos[k][j].a_x << "\n";
            errm << "the " << nearest_1_infos[k][j].a_x << "'th center is " << string_for_center(nearest_1_infos[k][j].a_x) << "\n";
            errm << "the recorded distance to this center is " << nearest_1_infos[k][j].d_x << "\n";

            throw zentas::zentas_error(errm.str());


            
          }
  
          if (e_first_nearest != nearest_1_infos[k][j].e_x){
            throw zentas::zentas_error(errm + "e_first_nearest != nearest_1_infos[k][j].e_x");
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
          throw zentas::zentas_error(errm);
        }
        ndata_internal += get_ndata(k);
      }
      
      custom_ndata_test();
      
      size_t ndata_target;

      if (in_refinement == true){
        ndata_target = ndata;
      }
      else{
        ndata_target = ndata - K;
      }
      
      
      if (ndata_internal != ndata_target){
        throw zentas::zentas_error(errm);
      }
    }

    void as_assigned_test(){
      for (size_t k = 0; k < K; ++k){
        for (size_t j = 0; j < get_ndata(k); ++j){
          if (nearest_1_infos[k][j].a_x != k){
            std::string errstring = "A sample in cluster " + std::to_string(k) + " has a_x " + std::to_string(nearest_1_infos[k][j].a_x);
            throw zentas::zentas_error(errstring);
          }
        }
      }
    }
    

    void injective_ID_test(){
      //mowri << "(injective_ID_test)" << zentas::Flush;
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
          std::stringstream errm;
          errm <<  "Error caught. ID "  << i << " appears not exactly once. ";
          errm <<  "It appears as centers of [";
          for (auto & x : as_center_ID){
            errm <<  " `" << x << "' ";
          }
          errm <<  " ], as well as member (k,j) {size of k}) in [";
          for (size_t q = 0; q < in_cluster.size(); ++ q){
            errm <<   " ("  << in_cluster[q] << "," << at_index[q] <<  " {"  << get_ndata(in_cluster[q]) <<  "}) ";
          }
          
          errm << "].";
          throw zentas::zentas_error(errm.str());
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
          throw zentas::zentas_error(errm_ss.str());
        }
      }
    }



};


template <class TMetric, class TData>
void BaseClusterer<TMetric, TData>::prepare_for_refinement(){
  //add center to data. /* TODO : can be faster at claransl123, as intercenter distances already have */ 
  for (size_t k = 0; k < K; ++k){
    put_sample_in_cluster(center_IDs[k]);
    /* hint for faster : final_push_into_cluster_basic(center_IDs[k], k, d1); */
  }
  
  in_refinement = true;
}





} //namespace nszen
 

#endif
