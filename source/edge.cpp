/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "edge.h"
#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;

edge::edge() {}

edge::edge(double* theLValue, double* theRValue, double theNormal) : itsLValue(theLValue), itsRValue(theRValue), itsNormal(theNormal) {}

void edge::evaluateFlux(const double a)
{
  double alpha = 0; //upwind flux
  itsFlux = a*(*itsLValue + *itsRValue)/2 + std::abs(a)*(1-alpha)/double(2)*(*itsLValue*itsNormal - *itsRValue*itsNormal);
}

double edge::Flux(double theNormal)
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

