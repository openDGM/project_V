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
#include "operators2D.h"
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
  int       itsN;                  // Polynomial order of the nodal element
  edge*     itsLEdge;              // Left element boundary edge
  edge*     itsREdge;              // Right element boundary edge
  double    itsLFlux,itsRFlux;     // Element boundary fluxes
  double    itsXLeft,itsXRight;    // Element boundary coordinates
  double    itsLNormal,itsRNormal; // Boundary normal vectors
  
  VectorXd itsU;                   // Solution of the nodal element
  VectorXd itsX;                   // Gglobal coordinates of the nodal element
  VectorXd itsRHSU;                // Right hand side of PDE
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

  // Set pointer to object of type edge to the left of element
  void setLEdge(edge* theLEdge);
  // Set pointer to object of type edge to the right of element
  void setREdge(edge* theREdge);

  // Set initial condition
  void setU(const VectorXd &theU);

  // Return solution
  VectorXd getU();

  // Return grid
  VectorXd getX();
};

class element2D
{
  private:
  int       itsN;                  // Polynomial order of the nodal element
  edge2D*   itsLeftEdge;           // Left edge of the element
  edge2D*   itsRightEdge;          // Right edge of the element
  edge2D*   itsTopEdge;            // Top edge of the element
  edge2D*   itsBottomEdge;         // Bottom edge of the element

  // Element boundry fluxes
  MatrixXd  itsLeftFlux,itsRightFlux,itsTopFlux,itsBottomFlux;

  // Element corner coordinates
  pair<double,double> itsBottomLeft, itsBottomRight, itsTopLeft, itsTopRight;

  // Element normals
  Vector2d  itsLeftNormal, itsRightNormal, itsTopNormal, itsBottomNormal;

  VectorXd itsU;                   // Solution of the nodal element
  VectorXd itsX;                   // Global X coordinates of the nodal element
  VectorXd itsY;                   // Global Y coordinates of the nodal element
  VectorXd itsRHSU;                // Right hand side of PDE
  VectorXd itsRESU;                // Runge-Kutta residual
  VectorXd itsDXDR;                // dxdr Jacobian of the nodal element
  VectorXd itsDXDS;                // dxds Jacobian of the nodal element
  VectorXd itsDYDR;                // dydr Jacobian of the nodal element
  VectorXd itsDYDS;                // dyds Jacobian of the nodal element
  VectorXd itsJ;                   // Jacobian dxdr*dyds-dxds*dydr

  public:
  element2D ();
  // Construct a new element from corner coordinates, normals and polynomial order
  element2D (pair<double,double> theBottomLeft, pair<double,double> theBottomRight, 
             pair<double,double> theTopLeft, pair<double,double> theTopRight, 
             Vector2d theLeftNormal, Vector2d theRightNormal, 
             Vector2d theTopNormal, Vector2d theBottomNormal, 
             int theN);

  // Update flux values at the boundaries
  void updateFluxes();

  // Compute right hand side of the PDE
  void advecRHS(const Vector2d a);

  // Integrate solution in time by one RK step
  void advecRK2D(const int INTRK, const double dt);

  // Create pointers to boundary values of the element
  unique_ptr<EdgeMap> connectToLeftBound();
  unique_ptr<EdgeMap> connectToRightBound();
  unique_ptr<EdgeMap> connectToTopBound();
  unique_ptr<EdgeMap> connectToBottomBound();

  // Set pointers to the edges around the element
  void setLeftEdge(edge2D* theLeftEdge);
  void setRightEdge(edge2D* theRightEdge);
  void setTopEdge(edge2D* theTopEdge);
  void setBottomEdge(edge2D* theBottomEdge);

  // Set initial condition for U
  void setU(const VectorXd &theU);

  // Return solution
  VectorXd getU();

  // Return X-grid
  VectorXd getX();

  // Return Y-grid
  VectorXd getY();

};
#endif
