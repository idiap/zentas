/*Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

zentas is a k-medoids library written in C++ and Python. This file is part of zentas.
zentas is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License version 3 as published by the Free Software Foundation.
zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details. You should have received a copy of
the GNU General Public License along with zentas. If not, see
<http://www.gnu.org/licenses/>.
*/
#ifndef ZENTAS_CLARANSL3_HPP
#define ZENTAS_CLARANSL3_HPP

#include "baseclaransl23.hpp"



namespace nszen{

class ClaransL3 : public BaseClaransL23 {

  public: 

  ClaransL3(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb): BaseClaransL23 (sb, eb) {}
    
  private:

  virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial){
    (void)serial;      
    return get_delta_hat_l3(k1, k2, j2, get_d_min_cc(k1), get_cc());
  }
};

}

#endif
