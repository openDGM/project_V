/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "edge.h"
#include <cmath>

using namespace Eigen;
using namespace std;


edge2D::edge2D() {} 

edge2D::edge2D(unique_ptr<EdgeMap> theLValue, unique_ptr<EdgeMap> theRValue, Vector2d theNormal) : 
    itsLValue(move(theLValue)), 
    itsRValue(move(theRValue)), 
    itsNormal(theNormal.normalized()) {}

void edge2D::evaluateFlux(const Vector2d a)
{
    double alpha = 1; //upwind flux

    auto mean = a*(*itsLValue + *itsRValue).transpose()/2;
    auto jump = itsNormal*a.cwiseAbs().transpose()*(itsNormal*itsLValue->transpose()-itsNormal*itsRValue->transpose())*(1-alpha)/2;
    itsFlux = (mean+jump).transpose()*itsNormal;
}

VectorXd edge2D::Flux(Vector2d theNormal)
{
    if (theNormal == itsNormal)
    {
      return itsFlux;
    }
    else
    {
      return -itsFlux;
    }
}

