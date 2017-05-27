#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP

#include <functional>
#include <random>
#include <sstream>

#include "outputwriter.hpp"

namespace nszen{

size_t extract_INT(std::string initialisation_method, size_t prefix_length);
  
template <class TMetric, class TDataIn>
void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const TDataIn & datain, TMetric & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy){


  size_t chain_length = extract_INT(initialisation_method, 8);
  
  mowri << "chain length : " << chain_length << zentas::Endl;
  
  /* number of attempted moves before sampling */
  
  double adistance;
  
  
  std::vector<bool> is_center (ndata, false);
  /* the first center is chosen at random */
  size_t index0 = dis(gen)%ndata;
  center_indices_init[0] = index0;
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
                          datain.at_for_metric(center_indices_init[kp]), 
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
                            datain.at_for_metric(center_indices_init[kp]), 
                            adistance);
        e_proposal = std::min<double>(e_proposal, f_energy(adistance));
      }

        
      if ((e_proposal/e_current)*((e0[sample_index_current]*2*ndata + e0_sum)/(e0[sample_index_proposal]*2*ndata + e0_sum))  >  dis_uni01(gen)){
        e_current = e_proposal;
        sample_index_current = sample_index_proposal;
      }
    }
    is_center[sample_index_current] = true;
    center_indices_init[k] = sample_index_current;          
  }
}


void populate_from_indices_init(const size_t * const center_indices_init_in, size_t * const center_indices_init, size_t K, size_t ndata);

void populate_uniformly(size_t * const center_indices_init, size_t K, size_t ndata, std::uniform_int_distribution<size_t> & dis, std::default_random_engine & gen);

}

#endif
