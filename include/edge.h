/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */


/*
 * Description:
 * The edge class is used to solve the
 * Riemann problem on the interface between
 * adjecent elements and return the flux.
 */

#include <Eigen/Dense>
#define EIGEN_MPL2_ONLY

#ifndef EDGE_H
#define EDGE_H

using namespace Eigen;

class edge
{
  private:
  double itsFlux;
  double* itsLValue;
  double* itsRValue;
  
  double itsNormal;

  public:
  edge();

  // Construct a new edge from pointers to left and right value and a edge normal
  edge(double* theLValue, double* theRValue, double theNormal);

  // Solve Riemann problem at the edge
  void evaluateFlux(const double a);

  // Return flux
  double Flux(double theNormal);
};
#endif 