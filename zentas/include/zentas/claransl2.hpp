// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_CLARANSL2_HPP
#define ZENTAS_CLARANSL2_HPP

#include <zentas/baseclaransl23.hpp>

namespace nszen
{

class ClaransL2 : public BaseClaransL23
{

  public:
  virtual std::string get_kmedoids_method_string() override final { return "clarans-2"; }

  ClaransL2(const SkeletonClustererInitBundle& sb, const ExtrasBundle& eb) : BaseClaransL23(sb, eb)
  {
  }

  private:
  virtual double get_delta_E(size_t k1, size_t k2, size_t j2, bool serial) override final
  {
    return get_delta_E_l2(k1, k2, j2, get_d_min_cc(k1), get_cc(), serial);
  }
};

}  // namespace nszen

#endif
