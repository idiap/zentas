#include "baseclusterer.hpp"


namespace nszen{

void P2Bundle::initialise_data(size_t n){
  ndata = n;
  d12_.reset(new std::array<double, 2> [n]);
  k_1_.reset(new size_t[n]);
  k_2_.reset(new size_t[n]);
  ori_.reset(new size_t[n]);
  
  for (unsigned i = 0; i < n; ++i){
    d_1(i) = std::numeric_limits<double>::max();
    d_2(i) = std::numeric_limits<double>::max();
  }
  for (unsigned i = 0; i < n; ++i){
    /* note that in the first round of kmeanspp, the lookup in cc requires this to be zero  */ 
    k_1(i) = 0; 
  }
}

P2Bundle::P2Bundle(size_t n){
  initialise_data(n);
}

P2Bundle::P2Bundle(const std::vector<size_t> & ori){
  initialise_data(ori.size());
  for (size_t i = 0; i < ndata; ++i){
    ori_[i] = ori[i];
  }
}


void XNearestInfo::reset(size_t new_a_x, double new_d_x, double new_e_x){
  a_x = new_a_x;
  d_x = new_d_x;
  e_x = new_e_x;
}

void XNearestInfo::reset(XNearestInfo & nearest_x_infos){
  a_x = nearest_x_infos.a_x;
  d_x = nearest_x_infos.d_x;
  e_x = nearest_x_infos.e_x;
}

XNearestInfo::XNearestInfo(size_t a_x, double d_x, double e_x):a_x(a_x), d_x(d_x), e_x(e_x){}

std::string XNearestInfo::get_string(){
  return  std::to_string(a_x) + "\t " + std::to_string(d_x) + "\t " + std::to_string(e_x);
}

}
