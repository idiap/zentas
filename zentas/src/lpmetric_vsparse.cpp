#include "lpmetric.hpp"


namespace nszen{




template <typename TNumber>
void L_oo_Distance<TNumber>::set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) {
          
  ++ncalcs;
  distance = 0;
  size_t a_pos = 0;
  size_t b_pos = 0;
  bool indices_remain = (a.size > 0 && b.size > 0);
  double diff;
  while (indices_remain && distance <= threshold){
    if (a.indices[a_pos] == b.indices[b_pos]){
      diff = a.values[a_pos] - b.values[b_pos];
      ++a_pos;
      ++b_pos;
      indices_remain = (a_pos < a.size && b_pos < b.size);
    }
    
    else if (a.indices[a_pos] < b.indices[b_pos]){
      diff = a.values[a_pos];
      ++a_pos;
      indices_remain = a_pos < a.size;
    }
    
    else{
      diff = a.values[b_pos];
      ++b_pos;
      indices_remain = b_pos < b.size;
    }
    distance = std::max(distance, std::abs(diff));
    calccosts += 1;
  }


  
  while (a_pos != a.size && distance < threshold){
    diff = a.values[b_pos];
    distance = std::max(distance, std::abs(diff));
    calccosts += 1;
    ++a_pos;
  }
  
  while (b_pos != b.size and distance < threshold){
    diff = b.values[b_pos];
    distance = std::max(distance, std::abs(diff));
    calccosts += 1;
    ++b_pos;        
  }            
}

template void L_oo_Distance<float>::set_distance(const SparseVectorSample<float> & a, const SparseVectorSample<float> & b, double threshold, double & distance);
template void L_oo_Distance<double>::set_distance(const SparseVectorSample<double> & a, const SparseVectorSample<double> & b, double threshold, double & distance);


template <typename TNumber>
void L0Distance<TNumber>::set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
  ++ncalcs;
  distance = 0;
  size_t a_pos = 0;
  size_t b_pos = 0;
  bool indices_remain = (a.size > 0 && b.size > 0);
  
  //threshold += 0; //why :| :| :| ?
  while (indices_remain && distance <= threshold){
    if (a.indices[a_pos] == b.indices[b_pos]){
      distance += (a.values[a_pos] != b.values[b_pos]);
      ++a_pos;
      ++b_pos;
      indices_remain = (a_pos < a.size && b_pos < b.size);
    }
    
    else if (a.indices[a_pos] < b.indices[b_pos]){
      distance += 1;
      ++a_pos;
      indices_remain = a_pos < a.size;
    }
    
    else{
      distance += 1;
      ++b_pos;
      indices_remain = b_pos < b.size;
    }
    calccosts += 1;
  }
  
  if (distance <= threshold){
    distance += a.size - a_pos;
    distance += b.size - b_pos;
    calccosts += 1;
  }
  

}

template void L0Distance<float>::set_distance(const SparseVectorSample<float> & a, const SparseVectorSample<float> & b, double threshold, double & distance);
template void L0Distance<double>::set_distance(const SparseVectorSample<double> & a, const SparseVectorSample<double> & b, double threshold, double & distance);

template <typename TNumber>
void L1Distance<TNumber>::set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance) {
          
  ++ncalcs;
  distance = 0;      
  double diff;
  size_t a_pos = 0;
  size_t b_pos = 0;
  
  /* if both a and b have unchecked indices remaining */
  bool indices_remain = (a.size > 0 && b.size > 0); 

  while (indices_remain && distance <= threshold){
    if (a.indices[a_pos] == b.indices[b_pos]){
      diff = a.values[a_pos] - b.values[b_pos];
      ++a_pos;
      ++b_pos;
      indices_remain = (a_pos < a.size && b_pos < b.size);
    }
    
    else if (a.indices[a_pos] < b.indices[b_pos]){
      diff = a.values[a_pos];
      ++a_pos;
      indices_remain = a_pos < a.size;
    }
    
    else{
      diff = b.values[b_pos];
      ++b_pos;
      indices_remain = b_pos < b.size;
    }
    
    distance += std::abs(diff);
    calccosts += 1;
  }
  
  
  /* either or a or b or both have exhausted their indices, noe make both exhausted */
  while (a_pos != a.size && distance <= threshold){
    distance += std::abs(a.values[a_pos]);
    calccosts += 1;
    ++a_pos;
  }
  
  while (b_pos != b.size && distance <= threshold){
    distance += std::abs(b.values[b_pos]);
    calccosts += 1;
    ++b_pos;        
  }      
}


template void L1Distance<float>::set_distance(const SparseVectorSample<float> & a, const SparseVectorSample<float> & b, double threshold, double & distance);
template void L1Distance<double>::set_distance(const SparseVectorSample<double> & a, const SparseVectorSample<double> & b, double threshold, double & distance);

    
template <typename TNumber>
void L2Distance<TNumber>::set_distance(const SparseVectorSample<TNumber> & a, const SparseVectorSample<TNumber> & b, double threshold, double & distance){
          
  ++ncalcs;
  distance = 0;      
  double diff;
  size_t a_pos = 0;
  size_t b_pos = 0;
  bool indices_remain = (a.size > 0 && b.size > 0);

  threshold *= threshold;
  
  while (indices_remain && distance <= threshold){
    if (a.indices[a_pos] == b.indices[b_pos]){
      diff = a.values[a_pos] - b.values[b_pos];
      ++a_pos;
      ++b_pos;
      indices_remain = (a_pos < a.size && b_pos < b.size);
    }
    
    else if (a.indices[a_pos] < b.indices[b_pos]){
      diff = a.values[a_pos];
      ++a_pos;
      indices_remain = a_pos < a.size;
    }
    
    else{
      diff = b.values[b_pos];
      ++b_pos;
      indices_remain = b_pos < b.size;
    }
    
    distance += diff*diff;
    calccosts += 1;
  }
  
  while (a_pos != a.size and distance < threshold){
    diff = a.values[a_pos];
    distance += diff*diff;
    calccosts += 1;
    ++a_pos;
  }
  
  while (b_pos != b.size and distance < threshold){
    diff = b.values[b_pos];
    distance += diff*diff;
    calccosts += 1;
    ++b_pos;        
  }
  
  distance = std::sqrt(distance);
}

template void L2Distance<float>::set_distance(const SparseVectorSample<float> & a, const SparseVectorSample<float> & b, double threshold, double & distance);
template void L2Distance<double>::set_distance(const SparseVectorSample<double> & a, const SparseVectorSample<double> & b, double threshold, double & distance);

}

