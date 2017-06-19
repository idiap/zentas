#include "lpmetric.hpp"

namespace nszen{

template <typename TNumber>
void L0Distance<TNumber>::set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) {
  ++ncalcs;
  distance = 0;       
  size_t d;
  for (d = 0; d < dimension; ++d){
    distance += (*(a + d) != *(b + d));
    if (distance > threshold){
      break;
    }
  }
  calccosts += (d == dimension ? d : 1 + d);
}
template void L0Distance<float>::set_distance(const float * const & a, const float * const & b, double threshold, double & distance);
template void L0Distance<double>::set_distance(const double * const & a, const double * const & b, double threshold, double & distance);

template <typename TNumber>
void L1Distance<TNumber>::set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance){
  ++ncalcs;
  distance = 0;       
  double diff;
  size_t d;
  
  for (d = 0; d < dimension; ++d){
    diff = *(a + d) - *(b + d);
    distance += std::abs(diff);
    if (distance > threshold){
      break;
    }
  }
  calccosts += (d == dimension ? d : 1 + d);
}
template void L1Distance<float>::set_distance(const float * const & a, const float * const & b, double threshold, double & distance);
template void L1Distance<double>::set_distance(const double * const & a, const double * const & b, double threshold, double & distance);

template <typename TNumber>
void L2Distance<TNumber>::set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) {
  
  /* Experiment with blas : 
   * with d = 1000, the speed-up was only 10s -> 7s. 
   * Not interesting enough to warrant the additional compilation hassle. 
   * it looked like this:
  >> ++ncalcs;
  >> wblas::copy(dimension, a, 1, worker, 1);
  >> wblas::axpy(dimension, static_cast<TNumber>(-1.), b, 1, worker, 1);
  >> distance = wblas::dot(dimension, worker, 1, worker, 1);
  >> calccosts += dimension;
  */
          
  ++ncalcs;
  distance = 0;       
  double diff;
  size_t d;

  /* Comment on squaring the threshold : We square the threshold 
   * as we'll be working with the squared distance. This might look
   *  dodgey: what if the threshold squared oveflows? This is not
   *  serious as it overflows to std::infinity, which is larger than
   *  any double anyway. Besides, it would be silly to complain 
   * here, as we're assuming that the distance squared is a valid 
   * double (which is a necessary assumption to make the computation).
   *  */
  threshold *= threshold;
  for (d = 0; d < dimension; ++d){
    diff = *(a + d) - *(b + d);
    distance += diff*diff;
    if (distance > threshold){
      break;
    }
  }
  
  calccosts += (d == dimension ? d : 1 + d);
  distance = std::sqrt(distance);

}
template void L2Distance<float>::set_distance(const float * const & a, const float * const & b, double threshold, double & distance);
template void L2Distance<double>::set_distance(const double * const & a, const double * const & b, double threshold, double & distance);


template <typename TNumber>
void L_oo_Distance<TNumber>::set_distance(const TNumber * const & a, const TNumber * const & b, double threshold, double & distance) {
  ++ncalcs;
  distance = 0;       
  size_t d;
  double diff;
  for (d = 0; d < dimension; ++d){
    diff = *(a + d) - *(b + d);
    distance = std::max(distance, std::abs(diff));
    if (distance > threshold){
      break;
    }
  }
  calccosts += (d == dimension ? d : 1 + d);
}
template void L_oo_Distance<float>::set_distance(const float * const & a, const float * const & b, double threshold, double & distance);
template void L_oo_Distance<double>::set_distance(const double * const & a, const double * const & b, double threshold, double & distance);

}
