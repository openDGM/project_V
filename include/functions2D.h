/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * A collection of functions needed to create the Nth order basis in for two dimensional nodal elements in the DGM discretization.
 *
 * N : refers to the order of the polynomial base
 * r : refers to location of a single node in the reference element on interval -1 <= r <= 1
 * x : refers to location of a single node in the generic element on interval XLeft <= x <= XRight
 */

#include <Eigen/Dense>
#define EIGEN_MPL2_ONLY

#ifndef FUNCTIONS2D_H
#define FUNCTIONS2D_H

using namespace Eigen;

namespace functions2D{

// return jacobi polynomial of Nth order evaluated at locations r
ArrayXd jacobiP2D(const ArrayXd& r, const ArrayXd& s, int alpha, int beta, int N);

// return gradient in r-direction of 2D jacobi polynomial of Nth order evaluated at locations r,s
VectorXd GradRJacobiP2D(const VectorXd& r, const ArrayXd & s, int alpha, int beta, int N);

// return gradient in s-direction of 2D jacobi polynomial of Nth order evaluated at locations r,s
VectorXd GradSJacobiP2D(const VectorXd& r, const ArrayXd & s, int alpha, int beta, int N);

// return the Vandermonde matrix for the Nth order reference element 
MatrixXd Vandermonde2D(int N, const VectorXd& r, const VectorXd& s);

// return gradient of Vandermonde matrix for the Nth order reference element
MatrixXd GradRVandermonde2D(int N, const VectorXd& r, const VectorXd& s);

// return gradient of Vandermonde matrix for the Nth order reference element
MatrixXd GradSVandermonde2D(int N, const VectorXd& r, const VectorXd& s);

// return Differentiation matrix/operator for the Nth order reference element
MatrixXd DRMatrix2D(int N, const VectorXd& r, const VectorXd& s, const MatrixXd& V2D);

// return Differentiation matrix/operator for the Nth order reference element
MatrixXd DSMatrix2D(int N, const VectorXd& r, const VectorXd& s, const MatrixXd& V2D);
}
#endif
