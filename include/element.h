/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * Classes to represent a nodal elements.
 * Classes encapsulate local grid and the 
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

typedef Map<VectorXd,0,InnerStride<>> EdgeMap;

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

class element2D
{
  private:
  int       itsN;
  edge2D*   itsLeftEdge;
  edge2D*   itsRightEdge;
  edge2D*   itsTopEdge;
  edge2D*   itsBottomEdge;

  VectorXd  itsLeftFlux,itsRightFlux,itsTopFlux,itsBottomFlux;
  pair<double,double> itsBottomLeft, itsBottomRight, itsTopLeft, itsTopRight;
  Vector2d  itsLeftNormal, itsRightNormal, itsTopNormal, itsBottomNormal;

  VectorXd itsU;                   // solution of the nodal element
  VectorXd itsX;                   // global coordinates of the nodal element
  VectorXd itsY;
  VectorXd itsRHSU;                // right hand side of PDE
  VectorXd itsRESU;                // Runge-Kutta residual
  VectorXd itsDXDR;                // dxdr Jacobian of the nodal element
  VectorXd itsDXDS;                // dxds Jacobian of the nodal element
  VectorXd itsDYDR;                // dydr Jacobian of the nodal element
  VectorXd itsDYDS;                // dyds Jacobian of the nodal element

  public:
  element2D ();
  element2D (pair<double,double> theBottomLeft, pair<double,double> theBottomRight, 
             pair<double,double> theTopLeft, pair<double,double> theTopRight, 
             Vector2d theLeftNormal, Vector2d theRightNormal, 
             Vector2d theTopNormal, Vector2d theBottomNormal, 
             int theN);

  void updateFluxes();
  void advecRHS(const Vector2d a);
  void advecRK2D(const int INTRK, const double dt);

  unique_ptr<EdgeMap> connectToLeftBound();
  unique_ptr<EdgeMap> connectToRightBound();
  unique_ptr<EdgeMap> connectToTopBound();
  unique_ptr<EdgeMap> connectToBottomBound();

  void setLeftEdge(edge2D* theLeftEdge);
  void setRightEdge(edge2D* theRightEdge);
  void setTopEdge(edge2D* theTopEdge);
  void setBottomEdge(edge2D* theBottomEdge);

  void setU(const VectorXd &theU);

  VectorXd getU();
  VectorXd getX();
  VectorXd getY();

};
#endif
