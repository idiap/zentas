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
#ifndef ZENTAS_VORONOI1_HPP
#define ZENTAS_VORONOI1_HPP

#include "baseclusterer.hpp"

namespace nszen{
 
template <class TMetric, class TData>
class VoronoiL0 : public BaseClusterer<TMetric, TData> {
  public:
    
    typedef typename TData::DataIn DataIn;
    VoronoiL0(const BaseClustererInitBundle<DataIn, TMetric> & ib): BaseClusterer<TMetric, TData> (ib) {}
  
  public:
    using BaseClusterer<TMetric, TData>::mowri;
    using BaseClusterer<TMetric, TData>::K;
    using BaseClusterer<TMetric, TData>::ndata;
    using BaseClusterer<TMetric, TData>::f_energy;
    using BaseClusterer<TMetric, TData>::set_center_sample_distance;
    using BaseClusterer<TMetric, TData>::get_center_sample_distance;
    using BaseClusterer<TMetric, TData>::get_sample_sample_distance;
    using BaseClusterer<TMetric, TData>::get_E_total;
    using BaseClusterer<TMetric, TData>::get_cluster_energy;
    using BaseClusterer<TMetric, TData>::swap_center_with_sample;
    using BaseClusterer<TMetric, TData>::reset_nearest_info;
    using BaseClusterer<TMetric, TData>::get_e1;
    using BaseClusterer<TMetric, TData>::get_d1;
    using BaseClusterer<TMetric, TData>::get_ndata;
    using BaseClusterer<TMetric, TData>::get_base_summary_string;
    using BaseClusterer<TMetric, TData>::base_put_sample_in_cluster;
    using BaseClusterer<TMetric, TData>::default_initialise_with_kmeanspp;    

   
   private:     
    virtual inline void put_sample_custom_in_cluster(size_t, size_t, const double * const) final override{}
    virtual inline void reset_sample_custom(size_t, size_t, size_t, const double * const) final override{}
    virtual inline void custom_append(size_t, size_t, size_t) final override{}
    virtual inline void custom_replace_with_last(size_t, size_t) final override{}
    virtual inline void custom_replace_with(size_t, size_t, size_t, size_t) final override{}
    virtual inline void custom_remove_last(size_t) final override{}
    virtual inline void increment_custom_cluster_statistics(size_t, size_t) final override{}
    virtual void set_normalised_custom_cluster_statistics(size_t) final override{}
    virtual void set_to_zero_custom_cluster_statistics(size_t) final override{}
    virtual void set_center_center_info() final override{}
    virtual void update_center_center_info() final override{}
    virtual void custom_cluster_statistics_test() final override{}

    virtual void set_redistribute_order(std::vector<size_t> & redistribute_order) override final {
      std::iota(redistribute_order.begin(), redistribute_order.end(), 0);      
    }

    virtual void put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) final override {
    
      (void)k1;
      (void)k2;
      (void)d2;
      (void)e2;
    }

    virtual void initialise_with_kmeanspp() override final {
      default_initialise_with_kmeanspp();
    }


    virtual void put_sample_in_cluster(size_t i) override final {
      base_put_sample_in_cluster(i);
    }
        
    virtual std::string get_round_summary() final override {
      std::stringstream ss;
      ss << get_base_summary_string();
      return ss.str();
    }
    
    virtual void update_sample_info() override final{
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

    
    virtual bool update_centers() override final{
      
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
            E_prop += f_energy(get_sample_sample_distance(k, j_prop, j));
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
};

} //namespace nszen
 

#endif
