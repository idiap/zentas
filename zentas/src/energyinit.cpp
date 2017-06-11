#include "energyinit.hpp"

namespace nszen{
  
EnergyInitialiser::EnergyInitialiser():critical_radius(0), exponent_coeff(0) {}
EnergyInitialiser::EnergyInitialiser(double critical_radius, double exponent_coeff):critical_radius(critical_radius), exponent_coeff(exponent_coeff) {}
    
double EnergyInitialiser::get_critical_radius () const {
  return critical_radius;
}

double EnergyInitialiser::get_exponent_coeff() const {
  return exponent_coeff;
}
   
}
