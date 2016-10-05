/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "element.h"
#include "globals.h"
#include "functions1D.h"
#include "functions2D.h"

using namespace Eigen;
using namespace globals;
using namespace functions1D;
using namespace functions2D;

element2D::element2D() {}

element2D::element2D(pair<double,double> theBottomLeft, pair<double,double> theBottomRight,
                     pair<double,double> theTopLeft, pair<double,double> theTopRight,
                     Vector2d theLeftNormal, Vector2d theRightNormal,
                     Vector2d theTopNormal, Vector2d theBottomNormal,
                     int theN) :
  itsBottomLeft(theBottomLeft),
  itsBottomRight(theBottomRight),
  itsTopLeft(theTopLeft),
  itsTopRight(theTopRight), 
  itsLeftNormal(theLeftNormal),
  itsRightNormal(theRightNormal),
  itsTopNormal(theTopNormal),
  itsBottomNormal(theBottomNormal),
  itsN(theN) 
  {
    itsX = LocalX2D(jacobiGL(0,0,itsN),itsBottomLeft,itsTopRight);
    itsY = LocalY2D(jacobiGL(0,0,itsN),itsBottomLeft,itsTopRight);
    itsDXDR = Jacobian(itsX,Dr);
    itsDXDS = Jacobian(itsX,Ds);
    itsDYDR = Jacobian(itsY,Dr);
    itsDYDS = Jacobian(itsY,Ds);
    itsRESU = VectorXd::Zero((itsN+1)*(itsN+1));
    itsU = VectorXd::Zero((itsN+1)*(itsN+1));
  }

void element2D::updateFluxes()
  {
    itsLeftFlux = itsLeftEdge->Flux(itsLeftNormal);
    itsRightFlux = itsRightEdge->Flux(itsRightNormal);
    itsBottomFlux = itsBottomEdge->Flux(itsBottomNormal);
    itsTopFlux = itsTopEdge->Flux(itsTopNormal);
  }

Map<VectorXd,0,InnerStride<>> element2D::connectToLeftBound()
  {
    return Map<VectorXd,0,InnerStride<>>(itsU.data(),itsN+1,InnerStride<>(itsN+1));
  }

Map<VectorXd,0,InnerStride<>> element2D::connectToRightBound()
  {
    return Map<VectorXd,0,InnerStride<>>(itsU.data()+itsN,itsN+1,InnerStride<>(itsN+1));
  }

Map<VectorXd,0,InnerStride<>> element2D::connectToTopBound()
  {
    return Map<VectorXd,0,InnerStride<>>(itsU.data(),itsN+1,InnerStride<>(1));
  }

Map<VectorXd,0,InnerStride<>> element2D::connectToBottomBound()
  {
    return Map<VectorXd,0,InnerStride<>>(itsU.data()+ItsN*(itsN+1),itsN+1,InnerStride<>(1));
  }


void element2D::setLeftEdge(edge* theLeftEdge)
  {
    itsLeftEdge = theLeftEdge;
  }

void element2D::setRightEdge(edge* theRightEdge)
  {
    itsRightEdge = theRightEdge;
  }

void element2D::setTopEdge(edge* theTopEdge)
  {
    itsTopEdge = theTopEdge;
  }

void element2D::setBottomEdge(edge* theBottomEdge)
  {
    itsBottomEdge = theBottomEdge;
  }
 
void element2D::setU(const VectorXd &theU)
  { 
    itsU = theU;
  }

void element2D::advecRHS(const Vector2d a)
  {
    Vector2d du(-a*itsU[0]+itsLFlux*itsLNormal,a*itsU[itsN]-itsRFlux*itsRNormal);

    itsRHSU = -a/itsJ.array()*(globals::Dr*itsU).array() + (globals::Lift*du).array()/itsJ.array();
  }

void element2D::advecRK2D(const int INTRK, const double dt)
  {
    itsRESU = rk4a[INTRK]*itsRESU + dt*itsRHSU;
    itsU = itsU + rk4b[INTRK]*itsRESU;
  }

VectorXd element2D::getU()
  {
    return itsU;
  }

VectorXd element2D::getX()
  {
    return itsX;
  }

VectorXd elemendt2D::getY()
  {
    return itsY;
  }