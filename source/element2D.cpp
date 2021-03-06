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
    itsDXDR = operators2D::DDr(itsN,itsX);
    itsDXDS = operators2D::DDs(itsN,itsX);
    itsDYDR = operators2D::DDr(itsN,itsY);
    itsDYDS = operators2D::DDs(itsN,itsY);
    itsJ = itsDXDR.array()*itsDYDS.array() - itsDXDS.array()*itsDYDR.array();
    itsRESU = VectorXd::Zero((itsN+1)*(itsN+1));
    itsU = VectorXd::Zero((itsN+1)*(itsN+1));
  }

void element2D::updateFluxes()
  {
    itsLeftFlux = itsLeftEdge->Flux2D();
    itsRightFlux = itsRightEdge->Flux2D();
    itsBottomFlux = itsBottomEdge->Flux2D();
    itsTopFlux = itsTopEdge->Flux2D();
  }

unique_ptr<EdgeMap> element2D::connectToLeftBound()
  {
    return make_unique<EdgeMap>(itsU.data(),itsN+1,InnerStride<>(1));
  }

unique_ptr<EdgeMap> element2D::connectToRightBound()
  {
    return make_unique<EdgeMap>(itsU.data()+itsN*(itsN+1),itsN+1,InnerStride<>(1));
  }

unique_ptr<EdgeMap> element2D::connectToTopBound()
  {
    return make_unique<EdgeMap>(itsU.data()+itsN,itsN+1,InnerStride<>(itsN+1));
  }

unique_ptr<EdgeMap> element2D::connectToBottomBound()
  {
    return make_unique<EdgeMap>(itsU.data(),itsN+1,InnerStride<>(itsN+1));
  }

void element2D::setLeftEdge(edge2D* theLeftEdge)
  {
    itsLeftEdge = theLeftEdge;
  }

void element2D::setRightEdge(edge2D* theRightEdge)
  {
    itsRightEdge = theRightEdge;
  }

void element2D::setTopEdge(edge2D* theTopEdge)
  {
    itsTopEdge = theTopEdge;
  }

void element2D::setBottomEdge(edge2D* theBottomEdge)
  {
    itsBottomEdge = theBottomEdge;
  }
 
void element2D::setU(const VectorXd &theU)
  { 
    itsU = theU;
  }

void element2D::advecRHS(const Vector2d a)
  {
    VectorXd du((itsN+1)*4);
    EdgeMap B(itsU.data(),itsN+1,InnerStride<>(itsN+1));
    EdgeMap L(itsU.data(),itsN+1,InnerStride<>(1));
    EdgeMap T(itsU.data()+itsN,itsN+1,InnerStride<>(itsN+1));
    EdgeMap R(itsU.data()+itsN*(itsN+1),itsN+1,InnerStride<>(1));

    EdgeMap JJ(itsJ.data()+itsN*(itsN+1),itsN+1,InnerStride<>(1));

    du << ((L*a.transpose() - itsLeftFlux.transpose())*itsLeftNormal).array()/sqrt(JJ.array())/2,
          ((B*a.transpose() - itsBottomFlux.transpose())*itsBottomNormal).array()/sqrt(JJ.array())/2,
          ((R*a.transpose() - itsRightFlux.transpose())*itsRightNormal).array()/sqrt(JJ.array())/2,
          ((T*a.transpose() - itsTopFlux.transpose())*itsTopNormal).array()/sqrt(JJ.array())/2;

    itsRHSU = -a[0]*itsDYDS.array()/itsJ.array()*operators2D::DDr(itsN,itsU)-a[0]*itsDYDR.array()/(-itsJ.array())*operators2D::DDs(itsN,itsU)
              -a[1]*itsDXDS.array()/(-itsJ.array())*operators2D::DDr(itsN,itsU)-a[1]*itsDXDR.array()/itsJ.array()*operators2D::DDs(itsN,itsU)
              +operators2D::Lift2D(itsN,du);

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

VectorXd element2D::getY()
  {
    return itsY;
  }
