/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * A collection of functions needed to create the Nth order basis for the nodal elements in the DGM discretization.
 * Implementation based on the book "Nodal Discontinous Galerkin Methods: Algorithms, Analysis, and Applications",
 * ISBN 978-0-387-72065-4
 *
 * N : refers to the order of the polynomial base
 * r : refers to location of a single node in the reference element on interval -1 <= r <= 1
 * x : refers to location of a single node in the generic element on interval XLeft <= x <= XRight
 */

#include <Eigen/Dense>
#define EIGEN_MPL2_ONLY
#include <tuple>

using namespace Eigen;

namespace functions1D{

// return jacobi polynomial of Nth order evaluated at locations r
ArrayXd jacobiP(const ArrayXd& r, int alpha, int beta, int N);

// return gradient of jacobi polynomial of Nth order evaluated at locations r
VectorXd GradJacobiP(const VectorXd& r, int alpha, int beta, int N);

// return Gauss-Quadrature points of Nth order jacobi polynomial
std::pair<VectorXd, VectorXd> jacobiGQ(int alpha, int beta, int N);

// return Gauss-Lobatto points of Nth order jacobi polynomial
VectorXd jacobiGL(int alpha, int beta, int N);

// return the Vandermonde matrix for the Nth order reference element 
MatrixXd Vandermonde1D(int N, const VectorXd& r);

// return gradient of Vandermonde matrix for the Nth order reference element
MatrixXd GradVandermonde1D(int N, const VectorXd& r);

// return Differentiation matrix/operator for the Nth order reference element
MatrixXd DMatrix1D(int N, const VectorXd& r, const MatrixXd& V);

// return Jacobian between reference element and generic element
VectorXd Jacobian(const VectorXd& x, const MatrixXd& Dr); 

// Surface integral operator of the Nth order reference element
MatrixXd Lift1D(int N, const MatrixXd& V);

// Transformed generic element locations x through affine mapping coordinates of r
VectorXd LocalX1D(const ArrayXd& r, double XLeft, double XRight);

// inlined helper functions
inline double a(int n, int alpha, int beta);
inline double b(int n, int alpha, int beta);
}
