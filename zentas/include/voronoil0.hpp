// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_VORONOI1_HPP
#define ZENTAS_VORONOI1_HPP

#include "skeletonclusterer.hpp"
#include "extrasbundle.hpp"

namespace nszen{


 
class VoronoiL0 : public SkeletonClusterer {
 
  public:

  virtual std::string get_kmedoids_method_string() override final {return "Voronoi-0";}
        
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
  

};







} //namespace nszen




#endif
