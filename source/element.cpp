/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "element.h"
#include "globals.h"
#include "functions1D.h"

using namespace Eigen;
using namespace globals;
using namespace functions1D;

element::element() {}

element::element(double theXLeft, double theXRight, double theLNormal, double theRNormal, int theN) : 
  itsXLeft(theXLeft), 
  itsXRight(theXRight), 
  itsLNormal(theLNormal),
  itsRNormal(theRNormal),
  itsN(theN) 
  {
    itsX = LocalX1D(jacobiGL(0,0,itsN),itsXLeft,itsXRight);
    itsJ = Jacobian(itsX,Dr);
    itsRESU = VectorXd::Zero(itsN+1);
    itsU = VectorXd::Zero(itsN+1);
  }

void element::updateFluxes()
  {
    itsLFlux = itsLEdge->Flux(itsLNormal);
    itsRFlux = itsREdge->Flux(itsRNormal);
  }

double* element::connectToLeftBound()
  {
    return &itsU[0];
  }

double* element::connectToRightBound()
  {
    return &itsU[itsN];
  }

void element::setLEdge(edge* theLEdge)
  {
    itsLEdge = theLEdge;
  }

void element::setREdge(edge* theREdge)
  {
    itsREdge = theREdge;
  }

void element::setU(const VectorXd &theU)
  { 
    itsU = theU;
  }

void element::advecRHS(const double a)
  {
    Vector2d du(-a*itsU[0]+itsLFlux*itsLNormal,a*itsU[itsN]-itsRFlux*itsRNormal);

    itsRHSU = -a/itsJ.array()*(globals::Dr*itsU).array() + (globals::Lift*du).array()/itsJ.array();
  }

void element::advecRK1D(const int INTRK, const double dt)
  {
    itsRESU = globals::rk4a[INTRK]*itsRESU + dt*itsRHSU;
    itsU = itsU + globals::rk4b[INTRK]*itsRESU;
  }

VectorXd element::getU()
  {
    return itsU;
  }

VectorXd element::getX()
  {
    return itsX;
  }
