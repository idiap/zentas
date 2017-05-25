#ifndef ZEN_DISPATCH_HPP
#define ZEN_DISPATCH_HPP

namespace nszen{

template <typename TData, typename TMetric>
void dispatch(std::string algorithm, size_t level, const typename TData::InitBundle & datain_ib, size_t K, const size_t * const indices_init, std::string initialisation_method, size_t max_proposals, size_t seed, double maxtime, double minmE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t maxrounds, bool patient, std::string energy, const typename TMetric::Initializer & metric_initializer, const EnergyInitialiser & energy_initialiser){
  
  typedef typename TData::DataIn DataIn;  
  
  DataIn datain(datain_ib);

  BaseClustererInitBundle<DataIn, TMetric> ib(K, &datain, indices_init, initialisation_method, seed, maxtime, minmE, indices_final, labels, nthreads, maxrounds, energy, metric_initializer, energy_initialiser);
  
  BaseClaransInitBundle clib(max_proposals, patient);
  
  if (algorithm.compare("clarans") == 0){
    if (level == 0){
      //throw std::runtime_error("clarans l0 not enabled, grep seow340jkosdm4 and uncomment here to enable");      
      nszen::ClaransL0<TMetric , TData > cc (ib, clib);
      cc.go();
      
    }  
    if (level == 1){
      //throw std::runtime_error("clarans l1 not enabled, grep sdrfoweinsdima and uncomment here to enable");
      nszen::ClaransL1<TMetric, TData > cc (ib, clib);
      cc.go();
    }
    if (level == 2){
      //throw std::runtime_error("clarans l2 not enabled, grep sdmi4sdfsdollll and uncomment here to enable");
      nszen::ClaransL2<TMetric, TData > cc (ib, clib);
      cc.go();
    }
    if (level == 3){
      //throw std::runtime_error("clarans l3 not enabled, grep sdmi4sdfsdfollll and uncomment here to enable");
      nszen::ClaransL3<TMetric, TData > cc (ib, clib);
      cc.go();
    }
  }
  
  else if (algorithm.compare("voronoi") == 0){
    //throw std::runtime_error("voronoi not enabled, grep dfseimmgrfiddiddidiid and uncomment here to enable");
    if (level == 0){
      nszen::VoronoiL0<TMetric,  TData > cc (ib);
      cc.go();
    }
  }
}


/* This is the place to do all kinds of tests on the input: all user calls (R/Python/Terminal) will pass through this function */
template <typename TData, typename TMetric>
void zentas_base(
const typename TData::InitBundle & datain_ib, size_t K, const size_t * const indices_init, std::string initialisation_method, std::string algorithm, size_t level, size_t max_proposals, bool capture_output, std::string & text, size_t seed, double maxtime, double minmE, size_t * const indices_final, size_t * const labels, size_t nthreads, size_t maxrounds, bool patient, std::string energy, const typename TMetric::Initializer & metric_initializer, const EnergyInitialiser & energy_initialiser){

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

  if (energy_initialiser.get_critical_radius() <= 0 && energy.compare("squarepotential") == 0){
    throw std::runtime_error("critical radius <= 0 is not allowed for squarepotential energy");
  }
  
  else if (energy_initialiser.get_critical_radius() > 0 && energy.compare("squarepotential") != 0){
    throw std::runtime_error("critical radius > 0 is only allowed of with squarepotential energy");
  }

  if (energy_initialiser.get_exponent_coeff() <= 0 && energy.compare("exp") == 0){
    throw std::runtime_error("exponent_coeff <= 0 is not allowed for exp energy");
  }
  
  else if (energy_initialiser.get_exponent_coeff() > 0 && energy.compare("exp") != 0){
    throw std::runtime_error("exponent_coeff > 0 is only allowed with exp energy");
  }
  
  if (K <= 1){
    throw std::runtime_error("K > 1 is a strict requirement");
  }
  
  if (K >= datain_ib.ndata){
    throw std::runtime_error("K < ndata is a strict requirement");
  }

  /* checking for (algorithm, level) compatibility */
  bool algorithm_level_ok = false;  
  if (algorithm.compare("voronoi") == 0){
    if (level == 0){
      algorithm_level_ok = true;
    }
  }
  
  else if (algorithm.compare("clarans") == 0){
    if (level == 0 || level == 1 || level == 2 || level == 3){
      algorithm_level_ok = true;
    }
  }
    
  if (algorithm_level_ok == false){
    throw std::runtime_error("Something wrong with algorithm : " + algorithm + "  level : " + std::to_string(level));
  }  
  
  

  
  if (initialisation_method == "from_init_indices" && indices_init == nullptr){
    throw std::runtime_error(R"(initialisation_method == "from_init_indices" && indices_init == nullptr) is true)");
  }
  
  //if (initialisation_method != "from_init_indices"){  
    //indices_init = nullptr;
  //}
  

  dispatch <TData, TMetric> (algorithm, level, datain_ib, K, indices_init, initialisation_method, max_proposals, seed, maxtime, minmE, indices_final, labels, nthreads, maxrounds, patient, energy, metric_initializer, energy_initialiser);
  
  
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
