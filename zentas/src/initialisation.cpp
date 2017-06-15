#include "initialisation.hpp"


#include "tdatain.hpp"
#include "tmetric.hpp"

#include "zentaserror.hpp"

namespace nszen{

namespace init{

size_t extract_INT(std::string initialisation_method, size_t prefix_length){
  
  if (initialisation_method.size() < prefix_length + 1){
    std::stringstream errm_ss;
    errm_ss << "invalid initialisation_method " << initialisation_method << ". It is not of the form `prefix'-INT, where INT is positive."; 
    throw zentas::zentas_error(errm_ss.str());
  }
  
  std::string digit_substring = initialisation_method.substr(prefix_length, initialisation_method.size() - prefix_length);
  auto striter = digit_substring.begin(); 
  while (striter != digit_substring.end()){
    char x = *striter;
    if (std::isdigit(x) == false){
      std::stringstream errm_ss;
      errm_ss << "Unexpected character while attempting to extract integer from " << initialisation_method << " (" << digit_substring << ")" << ", `" << x << "'";
      throw zentas::zentas_error(errm_ss.str());
    }
    ++striter;
  }
  
  return std::stoi(digit_substring);
}




void populate_from_indices_init(const size_t * const center_indices_init_in, size_t * const center_indices_init, size_t K, size_t ndata){    
  for (size_t k = 0; k < K; ++k){
    /* confirm uniqueness and range */
    if (center_indices_init_in[k] >= ndata){
      throw zentas::zentas_error("initialising center index out of bounds in BaseClusterer constructor");
    }
    for (unsigned j = 0; j < k; ++j){
      if (center_indices_init_in[k] == center_indices_init_in[j]){
        throw zentas::zentas_error("initialising center index at j and k are the same in BaseClusterer constructor");
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
}
