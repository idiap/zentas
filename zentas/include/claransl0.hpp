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

#ifndef ZENTAS_CLARANS_HPP
#define ZENTAS_CLARANS_HPP



#include "baseclarans.hpp"



namespace nszen{

 
class ClaransL0 : public BaseClarans{

  public:

    virtual std::string get_kmedoids_method_string() override final {return "clarans-0";}    
    
    ClaransL0(const SkeletonClustererInitBundle & sb, const ExtrasBundle & eb);
    
   private: 
   
    /* There are no custom quantities to
     * monitor for a Sample other than those 
     * defined in BaseClarans, use functions 
     * defined in BaseClarans */
    virtual void put_sample_custom_in_cluster(size_t i, size_t k_nearest, const double * const distances) final override;
    //virtual void reset_sample_custom(size_t k, size_t j, size_t nearest_center, const double * const distances) final override;
    virtual void custom_append(size_t k_new, size_t k, size_t j) final override;
    virtual void custom_replace_with_last(size_t k, size_t j) final override;
    virtual void custom_replace_with(size_t k1, size_t j1, size_t k2, size_t j2) final override;
    virtual void custom_remove_last(size_t k) final override;
    virtual void update_sample_info() override final;
    /* We do not record any 
     * center-center variables */
    virtual void set_center_center_info() final override;
    virtual void update_center_center_info() final override;
    virtual void custom_acceptance_call() override final;
    virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial) override final;
    virtual void custom_info_test() override final;
    virtual void put_sample_in_cluster(size_t i) override final; 
};

} //namespace nszen
 

#endif
