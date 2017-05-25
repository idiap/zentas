#include "initialisation.hpp"


#include "tdatain.hpp"
#include "tmetric.hpp"


namespace nszen{


///* ahhhh !!!! how to do this cleanly */

//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const SparseVectorDataRootedIn <double> & datain, LpMetric<SparseVectorDataRootedIn <double> >  & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const SparseVectorDataRootedIn <float> & datain, LpMetric<SparseVectorDataRootedIn <float> >  & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);



//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const SparseVectorDataUnrootedIn <double> & datain, LpMetric<SparseVectorDataUnrootedIn <double> >  & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const SparseVectorDataUnrootedIn <float> & datain, LpMetric<SparseVectorDataUnrootedIn <float> >  & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const DenseVectorDataRootedIn <double> & datain, LpMetric<DenseVectorDataRootedIn <double> >  & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, const DenseVectorDataRootedIn <float> & datain, LpMetric<DenseVectorDataRootedIn <float> >  & metric, size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, 
//const DenseVectorDataUnrootedIn <double> & datain, 
//LpMetric<DenseVectorDataUnrootedIn <double> >  & metric, 
//size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, 
//const DenseVectorDataUnrootedIn <float> & datain, 
//LpMetric<DenseVectorDataUnrootedIn <float> >  & metric, 
//size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, 
//const StringDataUnrootedIn <int>  & datain, 
//LevenshteinMetric<StringDataUnrootedIn <int> >  & metric, 
//size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, 
//const StringDataUnrootedIn <char>  & datain, 
//LevenshteinMetric<StringDataUnrootedIn <char> >  & metric, 
//size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, 
//const StringDataRootedIn <int>  & datain, 
//LevenshteinMetric<StringDataRootedIn <int> >  & metric, 
//size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);


//template void populate_afk_mc2(std::string initialisation_method, size_t * const center_indices_init, 
//const StringDataRootedIn <char>  & datain, 
//LevenshteinMetric<StringDataRootedIn <char> >  & metric, 
//size_t K, size_t ndata, zentas::outputwriting::OutputWriter & mowri, std::default_random_engine & gen,  std::uniform_int_distribution<size_t> & dis, std::function<double(double)> & f_energy);






/* TODO : move to a cpp  */
void populate_from_indices_init(const size_t * const center_indices_init_in, size_t * const center_indices_init, size_t K, size_t ndata){    
  for (size_t k = 0; k < K; ++k){
    /* confirm uniqueness and range. TODO : do I want to do this? performance? */
    if (center_indices_init_in[k] >= ndata){
      throw std::runtime_error("initialising center index out of bounds in BaseClusterer constructor");
    }
    for (unsigned j = 0; j < k; ++j){
      if (center_indices_init_in[k] == center_indices_init_in[j]){
        throw std::runtime_error("initialising center index at j and k are the same in BaseClusterer constructor");
      }
    }
    
    center_indices_init[k] = center_indices_init_in[k];
  }
}


void populate_uniformly(size_t * const center_indices_init, size_t K, size_t ndata, std::uniform_int_distribution<size_t> & dis, std::default_random_engine & gen){
  bool accepted;
  size_t proposed_i;
  for (size_t k = 0; k < K; ++k){
    accepted = false;
    while (accepted == false){
      accepted = true;
      proposed_i = dis(gen)%ndata;
      
      for (size_t k_m = 0; k_m < k; ++k_m){
        if (center_indices_init[k_m] == proposed_i){
          accepted = false;
        }
      }
    }
    if (accepted == true){
      center_indices_init[k] = proposed_i;
    }
  }
}

}
