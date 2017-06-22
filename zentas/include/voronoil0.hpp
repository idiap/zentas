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

#include "skeletonclusterer.hpp"
#include "extrasbundle.hpp"

namespace nszen{


 
class VoronoiL0 : public SkeletonClusterer {
 
  public:
    
  VoronoiL0(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb); 
  private:     
  
  virtual void put_sample_custom_in_cluster(size_t, size_t, const double * const) final override{}
  virtual void reset_sample_custom(size_t, size_t, size_t, const double * const) final override{}
  virtual void custom_append(size_t, size_t, size_t) final override{}
  virtual void custom_replace_with_last(size_t, size_t) final override{}
  virtual void custom_replace_with(size_t, size_t, size_t, size_t) final override{}
  virtual void custom_remove_last(size_t) final override{}
  virtual void increment_custom_cluster_statistics(size_t, size_t) final override{}
  virtual void set_normalised_custom_cluster_statistics(size_t) final override{}
  virtual void set_to_zero_custom_cluster_statistics(size_t) final override{}
  virtual void set_center_center_info() final override{}
  virtual void update_center_center_info() final override{}
  virtual void custom_cluster_statistics_test() final override{}
  virtual void set_redistribute_order(std::vector<size_t> & redistribute_order) override final;
  virtual void put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) final override;

  virtual void initialise_with_kmeanspp() override final {
    default_initialise_with_kmeanspp();
  }

  virtual void put_sample_in_cluster(size_t i) override final {
    base_put_sample_in_cluster(i);
  }
      
  virtual std::string get_round_summary() final override;
  virtual void update_sample_info() override final;    
  virtual bool update_centers() override final;
  
  
  virtual void refine_sample_info() override final {default_refine_sample_info();}
  virtual void custom_refine_center_center_info() override final {} 
  virtual void custom_initialise_refinement_variables() override final {optimised_refinement = false;}


};







} //namespace nszen




#endif
