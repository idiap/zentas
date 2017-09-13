// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

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
