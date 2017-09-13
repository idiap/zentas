// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_CLARANSL3_HPP
#define ZENTAS_CLARANSL3_HPP

#include "baseclaransl23.hpp"



namespace nszen{

class ClaransL3 : public BaseClaransL23 {

  public: 

  virtual std::string get_kmedoids_method_string() override final {return "clarans-3";}    
    
  ClaransL3(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb): BaseClaransL23 (sb, eb) {}
    
  private:

  virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial){
    (void)serial;      
    return get_delta_hat_l3(k1, k2, j2, get_d_min_cc(k1), get_cc());
  }
};

}

#endif
