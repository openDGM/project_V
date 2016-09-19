/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * A Class to represent a nodal element.
 * Class encapsulates local grid and the 
 * solution of the element as well as 
 * advection procedures.
 */

#include "edge.h"
#include <Eigen/Dense>
#define EIGEN_MPL2_ONLY

#ifndef ELEMENT_H
#define ELEMENT_H

using namespace Eigen;
using namespace std;

class element
{
  private:
  int       itsN;                  // polynomial order of the nodal element
  edge*     itsLEdge;              // left element boundary edge
  edge*     itsREdge;              // right element boundary edge
  double    itsLFlux,itsRFlux;     // element boundary fluxes
  double    itsXLeft,itsXRight;    // element boundary coordinates
  double    itsLNormal,itsRNormal; // boundary normal vectors
  
  VectorXd itsU;                   // solution of the nodal element
  VectorXd itsX;                   // global coordinates of the nodal element
  VectorXd itsRHSU;                // right hand side of PDE
  VectorXd itsRESU;                // Runge-Kutta residual
  VectorXd itsJ;                   // Jacobian of the nodal element

  public:
  element ();
  // Construct a new element from left and right edge coordinates, normals and polynomial order 
  element (double theXLeft, double theXRight, double theLNormal, double theRNormal, int theN);

  // Update flux values at boundaries from left and right edge
  void updateFluxes();
  // Compute right hand side of the PDE
  void advecRHS(const double a);
  // Integrate solution in time by one RK-Step
  void advecRK1D(const int INTRK, const double dt);

  // Return a pointer to leftmost value of solution vector itsU
  double* connectToLeftBound();
  // Return a pointer to rightmost value of solution vector itsU
  double* connectToRightBound();

  // set pointer to object of type edge to the left of element
  void setLEdge(edge* theLEdge);
  // set pointer to object of type edge to the right of element
  void setREdge(edge* theREdge);

  // set initial condition
  void setU(const VectorXd &theU);

  // Return solution
  VectorXd getU();

  // Return grid
  VectorXd getX();
};
#endif