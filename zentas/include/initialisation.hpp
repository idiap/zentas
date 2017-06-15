#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP

#include <functional>
#include <random>
#include <sstream>

#include "outputwriter.hpp"

namespace nszen{

namespace init{
size_t extract_INT(std::string initialisation_method, size_t prefix_length);
  
void populate_from_indices_init(const size_t * const center_indices_init_in, size_t * const center_indices_init, size_t K, size_t ndata);

void populate_uniformly(size_t * const center_indices_init, size_t K, size_t ndata, std::uniform_int_distribution<size_t> & dis, std::default_random_engine & gen);

}
}

#endif
