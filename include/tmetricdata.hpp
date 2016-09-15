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
/* defines classes satisying TMetricData template:
 * public: ndata, set_distance and get_distance.
 * */

#ifndef ZENTAS_METRICDATA_H
#define ZENTAS_METRICDATA_H

namespace nszen{

template <typename TNumber>
class VDataL2{
  public:

    size_t ndata;
    VDataL2(size_t ndata, size_t dimension, const TNumber * const data):ndata(ndata), dimension(dimension), data(data){}

    void set_distance(size_t i, size_t ip, double & distance) const{
      distance = 0;
      double diff = 0;
      for (size_t d = 0; d < dimension; ++d){
       diff = data[i*dimension + d] - data[ip*dimension + d];
       distance += diff*diff;
      }
      distance = std::sqrt(distance);
    }
    
    double get_distance(size_t i, size_t ip) const{
      double adistance = 0;
      set_distance(i, ip, adistance);
      return adistance;
    }

  private:
    size_t dimension;
    const TNumber * data;   
  
};

} //namespace nszen

#endif
