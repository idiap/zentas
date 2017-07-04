#include "lpmetric.hpp"

namespace nszen{



    
LpMetricInitializer::LpMetricInitializer(char p_, bool do_refinement_, std::string rf_alg_, size_t rf_max_rounds_, double rf_max_time_): p(p_), do_refinement(do_refinement_), rf_alg(rf_alg_), rf_max_rounds(rf_max_rounds_), rf_max_time(rf_max_time_) {}
    
LpMetricInitializer::LpMetricInitializer(): p('!') {}
    
void LpMetricInitializer::reset(std::string metric, bool do_refinement_, std::string rf_alg_, size_t rf_max_rounds_, double rf_max_time_){
  
  do_refinement = do_refinement_;
  rf_alg = rf_alg_;
  rf_max_rounds = rf_max_rounds_;
  rf_max_time = rf_max_time_;

  
  if (metric.compare("li") == 0){
    p = 'i';
  }
    
  else if (metric.compare("l2") == 0){
    p = '2';
  }
  
  else if (metric.compare("l1") == 0){
    p = '1';
  }

  else if (metric.compare("l0") == 0){
    p = '0';
  }
  
  else{
    std::stringstream ss;
    ss << "Currently, only li (inf) & l0 & l1 & l2 metrics are implemented for vector data, not `" << metric << "'."; 
    throw zentas::zentas_error(ss.str());
  }
}





}






