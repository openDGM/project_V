/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#define EIGEN_MPL2_ONLY

#ifndef GLOBALS_H
#define GLOBALS_H

using namespace Eigen;

namespace globals 
{
  typedef Map<VectorXd,0,InnerStride<>> EdgeMap;

  extern MatrixXd V;           // Vandermonde matrix
  extern MatrixXd invV;        // Inverse Vandermonde matrix
  extern MatrixXd DVr;         // Vandermonde gradient matrix
  extern MatrixXd Dr;          // Differentiation operator
  extern MatrixXd Lift;        // Surface flux integral operator

  extern MatrixXd V2D;
  extern MatrixXd DVs;
  extern SparseMatrix<double> Dr2D;
  extern SparseMatrix<double> Ds2D;
  extern SparseMatrix<double> Lift2D;

  // 4th order Runge-Kutta coefficients
  const double rk4a[5] = {0.0,
                           -567301805773.0/1357537059087.0,
                          -2404267990393.0/2016746695238.0,
                          -3550918686646.0/2091501179385.0,
                          -1275806237668.0/842570457699.0};
  const double rk4b[5] = {1432997174477.0/9575080441755.0,
                          5161836677717.0/13612068292357.0,
                          1720146321549.0/2090206949498.0,
                          3134564353537.0/4481467310338.0,
                          2277821191437.0/14882151754819.0}; 
  const double rk4c[5] = {0.0,
                          1432997174477.0/9575080441755.0,
                          2526269341429.0/6820363962896.0,
                          2006345519317.0/3224310063776.0,
                          2802321613138.0/2924317926251.0};
}
#endif
