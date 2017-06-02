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
/* defines classes satisying TEnergyFunc template */

#ifndef ZENTAS_ENERGYFUNCS_H
#define ZENTAS_ENERGYFUNCS_H


namespace nszen{
  

  
class Cubic{
  public:
   inline double operator() (double distance) const{
     return distance*distance*distance;
   }
};

class Quadratic{
  public:
   inline double operator() (double distance) const{
     return distance*distance;
   }
};


class Identity{
  public:
   inline double operator() (double distance) const{
     return distance;
   }
};


class Log{
  public:
   inline double operator() (double distance) const{
     return std::log(1. + distance);
   }
};

class SquareRoot{
  public:
   inline double operator() (double distance) const{
     return std::sqrt(distance);
   }
};


class Exponential{
  public:
  
    Exponential(double lambda_in): lambda(lambda_in) {}
    
    inline double operator() (double distance) const{
      return std::exp(lambda * distance) - 1.;
    }
    
  private:
    const double lambda;
  
};


class SquarePotential{
  public:
  
    SquarePotential(double critical_radius_in): critical_radius(critical_radius_in) {}
    
    inline double operator() (double distance) const{
      /* for clarans, return distance > critical_radius would be fine. But for 
       * voronoi, there needs to be some infinitesimally small gradient for
       * centers to move. We thus give the |__| potential a very slight gradient. */
      
      return distance > critical_radius;
      
      //double tiny_scaled_distance = 0.000001*distance/critical_radius;
      
      //if (distance > critical_radius){
        //return 1 + tiny_scaled_distance;
      //}
      
      //else{
        //return tiny_scaled_distance;
      //}
    }
    
  private:
    const double critical_radius;
  
};

} //namespace nszen

#endif

