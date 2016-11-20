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

edge2D::edge2D(unique_ptr<EdgeMap> theLValue, unique_ptr<EdgeMap> theRValue, Vector2d theLNormal, Vector2d theRNormal) : 
    itsLValue(move(theLValue)), 
    itsRValue(move(theRValue)), 
    itsLNormal(theLNormal.normalized()),
    itsRNormal(theRNormal.normalized()) {}

void edge2D::evaluateFlux2D(const Vector2d a)
{
    double alpha = 1; //upwind flux

    auto mean = a*(*itsLValue + *itsRValue).transpose()/2;
    auto jump = abs(a.transpose()*itsLNormal)*(itsLNormal*itsLValue->transpose()+itsRNormal*itsRValue->transpose())*(1-alpha)/2;

    itsFlux2D = mean+jump;
}

MatrixXd edge2D::Flux2D()
{
    return itsFlux2D;
}
