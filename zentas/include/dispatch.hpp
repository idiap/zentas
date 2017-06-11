#ifndef ZEN_DISPATCH_HPP
#define ZEN_DISPATCH_HPP

#include "energyinit.hpp"

namespace nszen{

void scrutinize_input_1(const EnergyInitialiser & energy_initialiser, std::string energy, size_t K, std::string algorithm, size_t level, size_t ndata);

template <typename TData, typename TMetric>
void dispatch(std::string algorithm, size_t level, const typename TData::InitBundle & datain_ib, size_t K, const size_t * const indices_init, std::string initialisation_method, size_t max_proposals, size_t seed, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t max_rounds, bool patient, std::string energy, bool with_tests, const typename TMetric::Initializer & metric_initializer, const EnergyInitialiser & energy_initialiser, std::chrono::time_point<std::chrono::high_resolution_clock> bigbang){
  
  
  
  
  typedef typename TData::DataIn DataIn;  
  
  DataIn datain(datain_ib);

  BaseClustererInitBundle<DataIn, TMetric> ib(K, &datain, indices_init, initialisation_method, seed, max_time, min_mE, indices_final, labels, nthreads, max_rounds, energy, with_tests, metric_initializer, energy_initialiser, bigbang);
  
  BaseClaransInitBundle clib(max_proposals, patient);
  
  if (algorithm.compare("clarans") == 0){
    if (level == 0){
      //throw zentas::zentas_error("clarans l0 not enabled, grep seow340jkosdm4 and uncomment here to enable");      
      nszen::ClaransL0<TMetric , TData > cc (ib, clib);
      cc.go();
      
    }  
    if (level == 1){
      //throw zentas::zentas_error("clarans l1 not enabled, grep sdrfoweinsdima and uncomment here to enable");
      nszen::ClaransL1<TMetric, TData > cc (ib, clib);
      cc.go();
    }
    if (level == 2){
      //throw zentas::zentas_error("clarans l2 not enabled, grep sdmi4sdfsdollll and uncomment here to enable");
      nszen::ClaransL2<TMetric, TData > cc (ib, clib);
      cc.go();
    }
    if (level == 3){
      //throw zentas::zentas_error("clarans l3 not enabled, grep sdmi4sdfsdfollll and uncomment here to enable");
      nszen::ClaransL3<TMetric, TData > cc (ib, clib);
      cc.go();
    }
  }
  
  else if (algorithm.compare("voronoi") == 0){
    //throw zentas::zentas_error("voronoi not enabled, grep dfseimmgrfiddiddidiid and uncomment here to enable");
    if (level == 0){
      nszen::VoronoiL0<TMetric,  TData > cc (ib);
      cc.go();
    }
  }
}


/* This is the place to do all kinds of tests on the input: all user calls (R/Python/Terminal) will pass through this function */
template <typename TData, typename TMetric>
void zentas_base(
const typename TData::InitBundle & datain_ib, size_t K, const size_t * const indices_init, std::string initialisation_method, std::string algorithm, size_t level, size_t max_proposals, bool capture_output, std::string & text, size_t seed, double max_time, double min_mE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t max_rounds, bool patient, std::string energy, bool with_tests, const typename TMetric::Initializer & metric_initializer, const EnergyInitialiser & energy_initialiser, const std::chrono::time_point<std::chrono::high_resolution_clock> & bigbang){

  /* used during experiments to see if openblas worth the effort. Decided not. 
  //openblas_set_num_threads(1);
  */

  #ifndef COMPILE_FOR_R
  std::stringstream buffer;
  auto cout_buff = std::cout.rdbuf();
  if (capture_output == true){
    std::cout.rdbuf(buffer.rdbuf());
  }
  std::ofstream nowhere;
  #endif
  
  
  scrutinize_input_1(energy_initialiser, energy, K, algorithm, level, datain_ib.ndata);

  
  if (initialisation_method == "from_init_indices" && indices_init == nullptr){
    throw zentas::zentas_error(R"(initialisation_method == "from_init_indices" && indices_init == nullptr) is true)");
  }
  
  //if (initialisation_method != "from_init_indices"){  
    //indices_init = nullptr;
  //}
  

  dispatch <TData, TMetric> (algorithm, level, datain_ib, K, indices_init, initialisation_method, max_proposals, seed, max_time, min_mE, indices_final, labels, nthreads, max_rounds, patient, energy, with_tests, metric_initializer, energy_initialiser, bigbang);
  
  
  #ifndef COMPILE_FOR_R
  if (capture_output == true){
    text = buffer.str();
    std::cout.rdbuf(cout_buff);
  }
  
  else{
    text = "capture_output was false, so nothing here";
  }
  #endif
}

}

#endif
