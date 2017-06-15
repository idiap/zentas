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

namespace nszen{


template <class TDataIn, class TMetric>
struct BaseClustererInitBundle{

     
  typedef typename TMetric::Initializer TMetricInitializer;

  const SkeletonClustererInitBundle & sc;
  const TDataIn & datain;
  const TMetricInitializer & metric_initializer;
  
  BaseClustererInitBundle(

  const SkeletonClustererInitBundle & sc, 
  const TDataIn & datain,
  const TMetricInitializer & metric_initializer): sc(sc), datain(datain), metric_initializer(metric_initializer) {}
  
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

  TData centers_data;
  RefinementCenterData rf_data;
  std::vector<TData> cluster_datas; 
  const DataIn * const ptr_datain;
  TMetric metric;
  std::unique_ptr<TData> ptr_kmoo_c_dt; 
  std::vector<TData> p2buns_dt;

/////////////////////////////////////
/////////// public //////////////////
/////////////////////////////////////
  public:
  
      
  BaseClusterer(const SkeletonClustererInitBundle & sc, const DataIn & datain, const TMetricInitializer & metric_initializer): 
  SkeletonClusterer(sc), 
  centers_data(datain, true), rf_data(datain, true), ptr_datain(& datain), metric(datain, sc.nthreads, metric_initializer), ptr_kmoo_c_dt( new TData(datain, true)) {}

  BaseClusterer(const BaseClustererInitBundle<DataIn, TMetric> & ib): BaseClusterer(ib.sc, ib.datain, ib.metric_initializer) {}

  /* metric only appears in these functions */
  virtual double get_rel_calccosts() override final{
    return metric.get_rel_calccosts();
  }
  
  virtual size_t get_ncalcs() override final{
    return metric.get_ncalcs();
  }

  virtual void set_cc_pp_distance(size_t k1, size_t k2)  override final  {
    metric.set_distance(ptr_kmoo_c_dt->at_for_metric(k1), ptr_kmoo_c_dt->at_for_metric(k2), std::numeric_limits<double>::max(), kmoo_cc[k1*K + k2]);
  }

  virtual void set_center_sample_pp_distance(size_t k, size_t bin, size_t i, double & adistance)  override final  {
    metric.set_distance(ptr_kmoo_c_dt->at_for_metric(k), p2buns_dt[bin].at_for_metric(i), aq2p_p2buns[bin].d_2(i), adistance);
  }
  
  virtual void set_center_sample_distance(size_t k, size_t k1, size_t j1, double threshold, double & distance)  override final   {
    metric.set_distance(centers_data.at_for_metric(k), cluster_datas[k1].at_for_metric(j1), threshold, distance);
  }

  virtual void set_center_sampleID_distance(size_t k, size_t i, double threshold, double & distance)  override final  {
    metric.set_distance(centers_data.at_for_metric(k), ptr_datain->at_for_metric(i), threshold, distance);
  }

  virtual void set_sampleID_sampleID_distance(size_t i1, size_t i2, double threshold, double & distance)  override final  {
    metric.set_distance(ptr_datain->at_for_metric(i1), ptr_datain->at_for_metric(i2), threshold, distance);
  }

  virtual void set_sample_sample_distance(size_t k1, size_t j1, size_t k2, size_t j2, double threshold, double & adistance)  override final   {
    metric.set_distance(cluster_datas[k1].at_for_metric(j1), cluster_datas[k2].at_for_metric(j2), threshold, adistance);
  }
  
  virtual void set_center_center_distance(size_t k1, size_t k2, double threshold, double & adistance)  override final  {
    metric.set_distance(centers_data.at_for_metric(k1), centers_data.at_for_metric(k2), threshold, adistance);
  }
    
  virtual void append_to_centers_from_ID(size_t i) override final{
      centers_data.append(ptr_datain->at_for_move(i));
  }
  
  virtual void add_empty_cluster() override final{
    cluster_datas.emplace_back(*ptr_datain, true);
  }

  protected:
  /* TData only appears in these functions */
  virtual std::string string_for_sample(size_t k, size_t j) override final {
    return cluster_datas[k].string_for_sample(j);
  }
  
  virtual std::string string_for_center(size_t k) override final {
    return centers_data.string_for_sample(k);
  }
  
  virtual void append_from_ID(size_t k, size_t i) override final {
    cluster_datas[k].append(ptr_datain->at_for_move(i));
  }
  
  virtual void append_across(size_t k_to, size_t k_from, size_t j_from) override final{      
    cluster_datas[k_to].append(cluster_datas[k_from].at_for_move(j_from));
  }

  virtual void append_pp_from_ID(size_t i) override final {
    ptr_kmoo_c_dt->append(ptr_datain->at_for_move(i));
  }

  virtual void append_pp_from_bin(size_t bin, size_t j) override final{
    ptr_kmoo_c_dt->append(p2buns_dt[bin].at_for_move(j));
  }
  
  virtual void append_aq2p_p2buns(size_t bin, size_t i) override final{
    p2buns_dt[bin].append(ptr_datain->at_for_move(i));
  }

 virtual size_t get_ndata(size_t k) override final{
    return cluster_datas[k].get_ndata();
  }
  
 
  virtual void swap_center_data(size_t k, size_t j) override final{
    nszen::swap<TData>(centers_data, k, cluster_datas[k], j);
  }
  
  virtual void replace_with_last_element(size_t k, size_t j) override final{
    replace_with_last<TData> (cluster_datas[k], j);      
  }
  
  virtual void remove_last(size_t k) override final{
    cluster_datas[k].remove_last();
  }
      
  virtual void  print_centers() override final{
    mowri << centers_data.get_string() << zentas::Flush;
  }

  virtual void replace_with(size_t k_to, size_t j_to, size_t k_from, size_t j_from) override final{
    cluster_datas[k_to].replace_with(j_to, cluster_datas[k_from].at_for_move(j_from));
  } 

  virtual size_t get_p2buns_dt_ndata(unsigned aq2p_bin) override final{
    return p2buns_dt[aq2p_bin].get_ndata();
  }

  virtual void reset_p2buns_dt(unsigned n_bins) override final{
    p2buns_dt = std::vector<TData> (n_bins,{*ptr_datain, true});
  }

  virtual void centers_replace_with_sample(size_t k1, size_t k2, size_t j2) override final{
    centers_data.replace_with(k1, cluster_datas[k2].at_for_move(j2));
  }

  virtual void kmoo_finish_with() override final {
    kmoo_cc.resize(0);
    kmoo_p2bun = P2Bundle(0);
    ptr_kmoo_c_dt.reset();// = TData(datain, true);        
  }          

  virtual std::string string_from_ID(size_t i) override final{
    return ptr_datain->string_for_sample(i);
  }
  



};

} //namespace nszen
 

#endif
