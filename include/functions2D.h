/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * A collection of functions needed to create the Nth order basis for two dimensional nodal elements in the DGM discretization.
 *
 * N : refers to the order of the polynomial in first dimension
 * M : refers to the order of the polynomial in second dimension
 * r : refers to location of node on 1st dimension in the reference element on interval -1 <= r <= 1,
 * s : refers to location of node on 2nd dimension in the reference element on interval -1 <= s <= 1
 
 */

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <utility>
#define EIGEN_MPL2_ONLY

#ifndef FUNCTIONS2D_H
#define FUNCTIONS2D_H

using namespace Eigen;
using namespace std;

namespace functions2D{

// return jacobi polynomial of Nth order evaluated at locations r
ArrayXd jacobiP2D(const ArrayXd& r, const ArrayXd& s, int alpha, int beta, int N, int M);

// return gradient in r-direction of 2D jacobi polynomial of Nth order evaluated at locations r,s
VectorXd GradRJacobiP2D(const VectorXd& r, const ArrayXd& s, int alpha, int beta, int N, int M);

// return gradient in s-direction of 2D jacobi polynomial of Nth order evaluated at locations r,s
VectorXd GradSJacobiP2D(const VectorXd& r, const ArrayXd& s, int alpha, int beta, int N, int M);

// return the Vandermonde matrix for the Nth order reference element 
MatrixXd Vandermonde2D(int N, const VectorXd& r, const VectorXd& s);

// return gradient of Vandermonde matrix for the Nth order reference element
MatrixXd GradRVandermonde2D(int N, const VectorXd& r, const VectorXd& s);

// return gradient of Vandermonde matrix for the Nth order reference element
MatrixXd GradSVandermonde2D(int N, const VectorXd& r, const VectorXd& s);

// return Differentiation matrix/operator for the Nth order reference element
SparseMatrix<double> DRMatrix2D(int N, const VectorXd& r, const VectorXd& s);

// return Differentiation matrix/operator for the Nth order reference element
SparseMatrix<double> DSMatrix2D(int N, const VectorXd& r, const VectorXd& s);

// Transformed generic element locations x through affine mapping coordinates of r
VectorXd LocalX2D(const ArrayXd& r, pair<double,double> BottomLeft, pair<double,double> TopRight);

// Transformed generic element locations y through affine mapping coordinates of s
VectorXd LocalY2D(const ArrayXd& s, pair<double,double> BottomLeft, pair<double,double> TopRight);

// Surface integral operator of the Nth order reference element
SparseMatrix<double> Lift2D(int N);

}
#endif
