
#include <vector>

#include <algorithm>
#include <map>

#include <sstream>
#include "energyinit.hpp"

#include "zentaserror.hpp"
namespace nszen{
  
void scrutinize_input_1(const EnergyInitialiser & energy_initialiser, std::string energy, size_t K, std::string algorithm, size_t level, size_t ndata){

  if (energy_initialiser.get_critical_radius() <= 0 && energy.compare("squarepotential") == 0){
    throw zentas::zentas_error("critical radius <= 0 is not allowed for squarepotential energy");
  }
  
  else if (energy_initialiser.get_critical_radius() > 0 && energy.compare("squarepotential") != 0){
    throw zentas::zentas_error("critical radius > 0 is only allowed of with squarepotential energy");
  }

  if (energy_initialiser.get_exponent_coeff() <= 0 && energy.compare("exp") == 0){
    throw zentas::zentas_error("exponent_coeff <= 0 is not allowed for exp energy");
  }
  
  else if (energy_initialiser.get_exponent_coeff() > 0 && energy.compare("exp") != 0){
    throw zentas::zentas_error("exponent_coeff > 0 is only allowed with exp energy");
  }
  
  if (K <= 1){
    throw zentas::zentas_error("K > 1 is a strict requirement");
  }
  
  if (K >= ndata){
    throw zentas::zentas_error("K < ndata is a strict requirement");
  }
  
  std::map <std::string, std::vector<size_t> > alg_levs = {
    {{"voronoi", {0}}, {"clarans", {0,1,2,3}}}};

  /* checking for (algorithm, level) compatibility */
  bool algorithm_level_ok = false;  
  
  std::vector<std::string> supported_algs {};
  for (auto & x : alg_levs){
    auto alg = std::get<0>(x);
    auto levs = std::get<1>(x);
    supported_algs.push_back(alg);
    if (algorithm.compare(alg) == 0){
      if (std::find(levs.begin(), levs.end(), level) != levs.end()){
        algorithm_level_ok = true;
      }
    }
  }
  
  //if (algorithm.compare("voronoi") == 0){
    //if (level == 0){
      //algorithm_level_ok = true;
    //}
  //}
  
  //else if (algorithm.compare("clarans") == 0){
    //if (level == 0 || level == 1 || level == 2 || level == 3){
      //algorithm_level_ok = true;
    //}
  //}
    
  if (algorithm_level_ok == false){
    std::stringstream errm;
    errm << "Something wrong with algorithm-level combination:    algorithm=" << algorithm << "  level=" << level << ". Either the algorithm does not exist, or the level of optimisation is not implemented for the algorithm. Currently the supported algorithms {levels} are: \n";
    for (auto & x : supported_algs){
      errm << x << " {";
      for (auto & l : alg_levs.at(x)){
        errm << " " << l << " ";
      }
      errm << "}\n";
    } 

    throw zentas::zentas_error(errm.str());
  }  
  
}

}
