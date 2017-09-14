// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#include <zentas/energyinit.hpp>

namespace nszen
{

EnergyInitialiser::EnergyInitialiser() : critical_radius(0), exponent_coeff(0) {}
EnergyInitialiser::EnergyInitialiser(double cr, double ec)
  : critical_radius(cr), exponent_coeff(ec)
{
}

double EnergyInitialiser::get_critical_radius() const { return critical_radius; }

double EnergyInitialiser::get_exponent_coeff() const { return exponent_coeff; }
}
