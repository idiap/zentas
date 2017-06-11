#ifndef ZENTAS_ENERGYINIT_HPP
#define ZENTAS_ENERGYINIT_HPP

namespace nszen{

struct EnergyInitialiser{

  private:
    double critical_radius;
    double exponent_coeff;
    
  public: 
    EnergyInitialiser();
    EnergyInitialiser(double critical_radius, double exponent_coeff);
    
    double get_critical_radius () const;
    
    double get_exponent_coeff() const;
   
};

}


#endif
