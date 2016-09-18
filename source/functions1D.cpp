/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "functions1D.h"
#include <Eigen/Eigenvalues>
#define EIGEN_MPL2_ONLY

#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace functions1D;

ArrayXd functions1D::jacobiP(const ArrayXd& r, int alpha, int beta, int N)
{
    // Zero order polynomial
    //-------------------------------------------------------------------------------------------------------
    double gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);
    if (N==0)
    {
        return ArrayXd::Ones(r.size())/sqrt(gamma0);
    }
    //-------------------------------------------------------------------------------------------------------

    // First order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    double gamma1 = (alpha+1)*(beta+1)/double(alpha+beta+3)*gamma0;
    if (N==1)
    {
       return ((alpha+beta+2)*r/2+(alpha-beta)/2)/sqrt(gamma1);
    }
    //-------------------------------------------------------------------------------------------------------
    
    // Recursion until first order is reached for higher order polynomials using recurrence equation
    // x*P_n(x) = a_n*P_(n-1)(x) + b_n*P_n(x) + a_(n+1)*P_(n+1)(x)
    // a_n = 2/(2n+alpha+beta)*sqrt[(n*(n+alpha+beta)*(n+alpha)*(n+beta))/((2n+alpha+beta-1)*(2n+alpha+beta+1))]
    // b_n = -(alpha^2-beta^2)/((2n+alpha+beta)*(2n+alpha+beta+2))
    //-------------------------------------------------------------------------------------------------------
    return (r*jacobiP(r,alpha,beta,N-1)-a(N-1,alpha,beta)*jacobiP(r,alpha,beta,N-2)-b(N-1,alpha,beta)*
            jacobiP(r,alpha,beta,N-1))/a(N,alpha,beta);
    //-------------------------------------------------------------------------------------------------------
}

VectorXd functions1D::GradJacobiP(const VectorXd& r, int alpha, int beta, int N)
{
    // Zero order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    if (N==0)
    {
        return VectorXd::Zero(r.size());
    }
    //-------------------------------------------------------------------------------------------------------

    // Nth order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    return sqrt(N*(N+alpha+beta+1))*jacobiP(r,alpha+1,beta+1,N-1);
    //-------------------------------------------------------------------------------------------------------
}

pair<VectorXd,VectorXd> functions1D::jacobiGQ(int alpha, int beta, int N)
{
    VectorXd r(N+1); // Quadrature locations
    VectorXd w(N+1); // Quadrature weights

    // Gauss quadrature points for Zero order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    if (N==0)
    {
        r[0] = (alpha-beta)/double(alpha+beta+2);
        w[0] = 2;
        return pair<VectorXd,VectorXd>(r,w);
    }
    //-------------------------------------------------------------------------------------------------------


    // Define diagonal matrix from recurrence equation (see jacobiP)
    //-------------------------------------------------------------------------------------------------------

    MatrixXd J = MatrixXd::Zero(N+1,N+1);

    Diagonal<MatrixXd,0> Jdiag0(J);
    Diagonal<MatrixXd,1> Jdiag1(J);

    for (int n=0; n<N+1; ++n)
    {
        Jdiag0[n] = b(n,alpha,beta)/2;
    }

    for (int n=1; n<N+1; ++n)
    {
        Jdiag1[n-1] = a(n,alpha,beta);
    }


    if (alpha + beta == 0)
    {
        J(0,0) = 0.0;
    }

    MatrixXd Jt = J;

    J += Jt.transpose();
    //-------------------------------------------------------------------------------------------------------

    // Solve Eigenvalue problem
    //-------------------------------------------------------------------------------------------------------
    EigenSolver<MatrixXd> es(J);
    //-------------------------------------------------------------------------------------------------------

    // extract Quadrature points and weights from eigenvalues and eigenvectors
    //-------------------------------------------------------------------------------------------------------
    r = es.eigenvalues().real();

    ArrayXd V;
    V = es.eigenvectors().real().row(0);
    w = V*V*pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)
        /tgamma(alpha+beta+1);
    //-------------------------------------------------------------------------------------------------------

    // Sort points and weights in ascending order
    //-------------------------------------------------------------------------------------------------------
    vector<pair<double,double>> sort_temp(N+1);
    for (int i=0; i<N+1; ++i)
    {
        sort_temp[i] = pair<double,double>(r[i],w[i]);
    }

    sort(sort_temp.begin(),sort_temp.end());

    for (int i=0; i<N+1; ++i)
    {
        r[i]=sort_temp[i].first;
        w[i]=sort_temp[i].second;
    }
    //-------------------------------------------------------------------------------------------------------
    return pair<VectorXd,VectorXd>(r,w);
}

VectorXd functions1D::jacobiGL(int alpha, int beta, int N)
{
    VectorXd r(N+1);
    // First order Gauss Lombetto points (boundaries only)
    //-------------------------------------------------------------------------------------------------------
    if (N == 1)
    {
        r << -1, 1;
        return r;
    }
    //-------------------------------------------------------------------------------------------------------

    // Nth order Gauss Lombetto points are Quadrature points surrounded by boundary points 
    //-------------------------------------------------------------------------------------------------------
    r << -1, jacobiGQ(alpha+1,beta+1,N-2).first,1;
    //-------------------------------------------------------------------------------------------------------
    return r;
}

MatrixXd functions1D::Vandermonde1D(int N, const VectorXd& r)
{
    MatrixXd V1D(r.size(), N+1);
    // Create Vandermonde matrix that contains the base functions of the reference element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<r.size(); ++i)
    {
        V1D.col(i) = jacobiP(r,0,0,i);
    }
    //-------------------------------------------------------------------------------------------------------
    return V1D;
}

MatrixXd functions1D::GradVandermonde1D(int N, const VectorXd& r)
{
    MatrixXd DVr(r.size(),N+1);
    // Create the GradVandermonde matrix that contains the gradients of base functions of the reference 
    // element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<r.size(); i++)
    {
        DVr.col(i) = GradJacobiP(r,0,0,i);
    }
    //-------------------------------------------------------------------------------------------------------
    return DVr;
}

MatrixXd functions1D::DMatrix1D(int N, const VectorXd& r, const MatrixXd& V)
{
    auto DVr = GradVandermonde1D(N,r);
    // Create the differentiation operator of the reference element
    //-------------------------------------------------------------------------------------------------------
    return DVr*V.inverse();
    //-------------------------------------------------------------------------------------------------------
}

VectorXd functions1D::Jacobian(const VectorXd& x, const MatrixXd& Dr)
{
    // Calculate Jacobian to transform between element and reference element
    //-------------------------------------------------------------------------------------------------------
    return Dr*x;
    //-------------------------------------------------------------------------------------------------------
}

MatrixXd functions1D::Lift1D(int N, const MatrixXd& V)
{
    // Create operator for surface flux integral on the reference element
    //-------------------------------------------------------------------------------------------------------
    MatrixXd Emat = MatrixXd::Zero(N+1,2);
    Emat(0,0)   = 1.0;
    Emat(N,1) = 1.0;

    return V*(V.transpose()*Emat);
    //-------------------------------------------------------------------------------------------------------
}

VectorXd functions1D::LocalX1D(const ArrayXd& r, double XLeft, double XRight)
{
    // Coordinate transformation from reference element on interval -1 <= x <= 1 to arbitrary interval
    // XLeft <= x <= XRight using affine mapping
    //-------------------------------------------------------------------------------------------------------
    return XLeft + 0.5*(r+1)*(XRight-XLeft);
    //-------------------------------------------------------------------------------------------------------
}

inline double functions1D::a(int n, int alpha, int beta)
{
    return 2/double(2*n+alpha+beta)*sqrt(n*(n+alpha+beta)*(n+alpha)*(n+beta)
           /double((2*n+alpha+beta-1)*(2*n+alpha+beta+1)));
}

inline double functions1D::b(int n, int alpha, int beta)
{
    return -(alpha*alpha-beta*beta)/double((2*n+alpha+beta)*(2*n+alpha+beta+2));
}
