#include "lpmetric.hpp"

namespace nszen{



    
LpMetricInitializer::LpMetricInitializer(char p_, bool do_refinement_): p(p_), do_refinement(do_refinement_){}
    
LpMetricInitializer::LpMetricInitializer(): p('!') {}
    
void LpMetricInitializer::reset(std::string metric, bool do_refinement_){
  
  do_refinement = do_refinement_;
  
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






