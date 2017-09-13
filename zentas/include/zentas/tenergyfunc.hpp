// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_ENERGYFUNCS_H
#define ZENTAS_ENERGYFUNCS_H

#include <cmath>

// TODO : how important are the s here? Can we move to cpp?
namespace nszen
{

class Cubic
{
  public:
  double operator()(double distance) const { return distance * distance * distance; }
};

class Quadratic
{
  public:
  double operator()(double distance) const { return distance * distance; }
};

class Identity
{
  public:
  double operator()(double distance) const { return distance; }
};

class Log
{
  public:
  double operator()(double distance) const { return std::log(1. + distance); }
};

class SquareRoot
{
  public:
  double operator()(double distance) const { return std::sqrt(distance); }
};

class Exponential
{
  public:
  Exponential(double lambda_in) : lambda(lambda_in) {}

  double operator()(double distance) const { return std::exp(lambda * distance) - 1.; }

  private:
  const double lambda;
};

class SquarePotential
{
  public:
  SquarePotential(double critical_radius_in) : critical_radius(critical_radius_in) {}

  double operator()(double distance) const
  {
    /* for clarans, return distance > critical_radius would be fine. But for
     * voronoi, there needs to be some infinitesimally small gradient for
     * centers to move. We thus give the |__| potential a very slight gradient. */

    return distance > critical_radius;
  }

  private:
  const double critical_radius;
};

}  // namespace nszen

#endif
