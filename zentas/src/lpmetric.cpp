#include "lpmetric.hpp"

namespace nszen{



    
LpMetricInitializer::LpMetricInitializer(char p): p(p){}
    
LpMetricInitializer::LpMetricInitializer(): p('!') {}
    
void LpMetricInitializer::reset(std::string metric){
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






