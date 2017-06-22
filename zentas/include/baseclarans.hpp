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

#ifndef ZENTAS_BASECLARANS_HPP
#define ZENTAS_BASECLARANS_HPP

#include "skeletonclusterer.hpp"
#include "extrasbundle.hpp"

namespace nszen{
//TODO : put this in a namespace. 
struct ClaransStatistics{
  public:

  // max_{j in k} (energy_margins)
  double M;
  // sum_{j in k} (energy_margins)
  double M_star;
  // max_{j in k} (d_first_nearest)
  double R1;
  // max_{j in k} (d_second_nearest)
  double R2;    
  // M / ndata in cluster k
  double m;
  // M_star / ndata in cluster k
  double m_star;
  
  ClaransStatistics():M(0.0), M_star(0.0), R1(0.0), R2(0.0), m(0.0), m_star(0.0){
  }
  
  void set_to_zero();
  void increment(double d_first_nearest, double e_first_nearest, double d_second_nearest, double e_second_nearest);
  void set_normalised_statistics(size_t ndata);
      
};


struct BaseClaransInitBundle{
  size_t max_proposals;
  bool patient;
  
  BaseClaransInitBundle(size_t max_proposals, bool patient): max_proposals(max_proposals), patient(patient) {}

};

 
class BaseClarans : public SkeletonClusterer { 

  private:
  
  std::vector<std::vector<XNearestInfo>> nearest_2_infos;
  
  public:
  
  std::vector<std::vector<double>> energy_margins;
  std::vector<ClaransStatistics> cluster_statistics;
    
  protected:
  
  size_t k_to;
  size_t k_from;
  size_t j_from;
  
  private: 
  
  size_t n_proposals;
  const size_t max_proposals;
  bool patient;
  
  
  public:
  
  BaseClarans(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb):
  SkeletonClusterer(sb),nearest_2_infos(sb.K), energy_margins(sb.K), cluster_statistics(sb.K), n_proposals(0), max_proposals(eb.clarans.max_proposals), patient(eb.clarans.patient) {      
  }
  
  size_t get_max_proposals(){
    return max_proposals;
  }
  
  size_t get_a2(size_t k, size_t j){
    return nearest_2_infos[k][j].a_x;
  }
  
  double get_d2(size_t k, size_t j){
    return nearest_2_infos[k][j].d_x;
  }    
  
  
  double get_e2(size_t k, size_t j){
    return nearest_2_infos[k][j].e_x;
  }
  
  bool update_centers_greedy();
  bool update_centers_patient();
  
  private:
  
  virtual void reset_sample_custom(size_t k, size_t j, size_t nearest_center, const double * const distances) final override;
  virtual void initialise_with_kmeanspp() override final;
  virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial) = 0;
  virtual void set_redistribute_order(std::vector<size_t> & redistribute_order) override final;
  virtual bool update_centers() override final;
  
  /* seems overly conservative using hoeffdinger algo : this function is not used (but implementation retained, at end) */ 
  double get_delta_E_hoeffding_l3(size_t k1, size_t k2, size_t j2, double d_nearest_k1, const double * const cc); 
  void reset_1_nearest_excluding(XNearestInfo & nearest_info, size_t k, const double * const distances);
  void set_center_center_distances(size_t k, double * const distances);
  double get_delta_E_not_k1_l2(size_t k1, size_t k2, size_t j2, const double * const cc, bool serial);
  void update_k_to_sample_info_l2(size_t j_a, size_t j_z, const double * const dists_centers_old_k_to, const double * const cc);
  void pll_update_sample_info_l23(size_t k, size_t j_a, size_t j_z, const double * const dists_centers_min_pr, const double * const cc);
  
  /* Determine the nearest and second nearest of sample j1 of cluster k1, given that distance to center a1 is d1 and distance to a2 is d2. It is not necessary that d1 < d2 at entry. */
  void set_nearest_12_warmstart(size_t k1, size_t j1, size_t & a1, size_t & a2, double & d1, double & d2, const double * const cc);
  void reset_sample_info_direct(size_t k, size_t j, size_t a1, size_t a2, double d1, double d2);
  double get_delta_E_k1_l12(size_t k1, size_t k2, size_t j2, double dist_k1_j2, bool serial);
  double get_delta_E_not_k1_l1(size_t k1, size_t k2, size_t j2, const double * const dists_centers_j2, bool serial);
  void acceptance_call(size_t k1, size_t k2, size_t j2);
  virtual std::string get_round_summary() override final;
  void reset_second_nearest_info(size_t k, size_t j, size_t k_second_nearest, double d_second_nearest, double e_second_nearest);
  virtual void custom_acceptance_call() = 0;
  void set_proposal(size_t & k1, size_t & k2, size_t & j2);
  void update_k_to_k_from_j_from(size_t k_to_in, size_t k_from_in, size_t j_from_in);
  void pll_update_sample_info_l1(size_t k, size_t j_a, size_t j_z, double dists_centers_min_pr_k);
  
  /* ************************************
  * *********** tests ******************
  * ********************************** */     
  virtual void custom_ndata_test() override final;
  void clarans_statistics_test();
  virtual void custom_cluster_statistics_test() override final;
  virtual void increment_custom_cluster_statistics(size_t k, size_t j) final override;
  virtual void set_normalised_custom_cluster_statistics(size_t k) final override;    
  virtual void set_to_zero_custom_cluster_statistics(size_t k) final override;
  void set_second_nearest(size_t k_first_nearest, const double * const distances, size_t & k_second_nearest, double & d_second_nearest);
  

  
  protected:
  
  void unset_clarans_variable_for_optimised_refinement();
  void refresh_energy_margins(size_t k, size_t j);
  double get_delta_hat_l3(size_t k1, size_t k2, size_t j2, double d_nearest_k1, const double * const cc);
  void center_center_info_test_l1(const std::vector<XNearestInfo> & center_nearest_center);
  void acceptance_call_l2(double * const cc, double * const dists_centers_old_k_to, double * const d_min_cc, size_t * const a_min_cc);
  void set_center_center_info_l2(double * const cc, double * const d_min_cc, size_t * const a_min_cc);
  void put_nearest_2_infos_margin_in_cluster_final(size_t k_first_nearest, size_t k_second_nearest, double d_second_nearest, double e_second_nearest);
  virtual void put_nearest_2_infos_margin_in_cluster_post_kmeanspp(size_t k1, size_t k2, double d2, double e2) final override;
  void put_nearest_2_infos_margin_in_cluster(size_t i, size_t k_first_nearest, const double * const distances);
  void reset_sample_nearest_2_infos_margin(size_t k, size_t j, size_t nearest_center, const double * const distances);
  void nearest_2_infos_margin_append(size_t k_new, size_t k, size_t j);    
  void nearest_2_infos_margin_replace_with_last(size_t k, size_t j);
  void nearest_2_infos_margin_replace_with(size_t k1, size_t j1, size_t k2, size_t j2);
  void nearest_2_infos_margin_remove_last(size_t k);    
  void nearest_2_infos_margin_test();    
  void basic_clarans_update_sample_info();
  double get_delta_E_l0(size_t k1, size_t k2, size_t j2);
  void update_sample_info_l1(const double * const dists_centers_old_k_to, const double * const dists_centers_new_k_to);
  void update_sample_info_l23(const double * const dists_centers_old_k_to, const double * const cc);
  double get_delta_E_l1(size_t k1, size_t k2, size_t j2, double d_nearest_k1, bool serial);
  double get_delta_E_l2(size_t k1, size_t k2, size_t j2, double d_nearest_k1, const double * const cc, bool serial);    
  void update_center_center_info_l1(std::vector<XNearestInfo> & center_nearest_center, double * const dists_centers_old_k_to, double * const dists_centers_new_k_to);
  void set_center_center_info_l1(std::vector<XNearestInfo> & center_nearest_center);

};


}




#endif
