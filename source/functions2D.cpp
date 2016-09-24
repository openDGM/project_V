/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "functions1D.h"
#include "functions2D.h"
#define EIGEN_MPL2_ONLY

#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace functions1D;
using namespace functions2D;

ArrayXd functions2D::jacobiP2D(const ArrayXd& r, const ArrayXd& s, int alpha, int beta, int N)
{
    ArrayXd jp2D((N+1)*(N+1));

    ArrayXd jpr = jacobiP(r,alpha,beta,N);
    ArrayXd jps = jacobiP(s,alpha,beta,N);

    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            jp2D[i*N+j] = jpr[i]*jps[j];
        }
    }

    return jp2D;
}

VectorXd functions2D::GradRJacobiP2D(const VectorXd& r, const ArrayXd& s, int alpha, int beta, int N)
{
    // Zero order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    if (N==0)
    {
        return VectorXd::Zero(r.size()*s.size());
    }
    //-------------------------------------------------------------------------------------------------------

    // Nth order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    ArrayXd jp2drgrad((N+1)*(N+1));

    ArrayXd jprgrad = sqrt(N*(N+alpha+beta+1))*jacobiP(r,alpha+1,beta+1,N-1);
    ArrayXd jps = jacobiP(s,alpha,beta,N);

    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            jp2drgrad[i*N+j] = jprgrad[i]*jps[j];
        }
    }

    return jp2drgrad;
    //-------------------------------------------------------------------------------------------------------
}

VectorXd functions2D::GradSJacobiP2D(const VectorXd& r, const ArrayXd& s, int alpha, int beta, int N)
{
    // Zero order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    if (N==0)
    {
        return VectorXd::Zero(r.size()*s.size());
    }
    //-------------------------------------------------------------------------------------------------------

    // Nth order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    ArrayXd jp2dsgrad((N+1)*(N+1));

    ArrayXd jpsgrad = sqrt(N*(N+alpha+beta+1))*jacobiP(s,alpha+1,beta+1,N-1);
    ArrayXd jpr = jacobiP(r,alpha,beta,N);

    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            jp2dsgrad[i*N+j] = jpr[i]*jpsgrad[j];
        }
    }

    return jp2dsgrad;
    //-------------------------------------------------------------------------------------------------------
}

MatrixXd functions2D::Vandermonde2D(int N, const VectorXd& r, const VectorXd& s)
{
    MatrixXd V2D(r.size()*s.size(), N+1);
    // Create Vandermonde matrix that contains the base functions of the reference element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<N+1; ++i)
    {
        V2D.col(i) = jacobiP2D(r,s,0,0,i);
    }
    //-------------------------------------------------------------------------------------------------------
    return V2D;
}

MatrixXd functions2D::GradRVandermonde2D(int N, const VectorXd& r, const VectorXd& s)
{
    MatrixXd DVr(r.size()*s.size(),N+1);
    // Create the GradVandermonde matrix that contains the gradients of base functions of the reference 
    // element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<N+1; i++)
    {
        DVr.col(i) = GradRJacobiP2D(r,s,0,0,i);
    }
    //-------------------------------------------------------------------------------------------------------
    return DVr;
}

MatrixXd functions2D::GradSVandermonde2D(int N, const VectorXd& r, const VectorXd& s)
{
    MatrixXd DVs(r.size()*s.size(),N+1);
    // Create the GradVandermonde matrix that contains the gradients of base functions of the reference 
    // element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<N+1; i++)
    {
        DVs.col(i) = GradSJacobiP2D(r,s,0,0,i);
    }
    //-------------------------------------------------------------------------------------------------------
    return DVs;
}

MatrixXd functions2D::DRMatrix2D(int N, const VectorXd& r, const VectorXd& s, const MatrixXd& V2D)
{
    auto DVr = GradRVandermonde2D(N,r,s);
    // Create the differentiation operator of the reference element
    //-------------------------------------------------------------------------------------------------------
    return DVr*V2D.inverse();
    //-------------------------------------------------------------------------------------------------------
}

MatrixXd functions2D::DSMatrix2D(int N, const VectorXd& r, const VectorXd& s, const MatrixXd& V2D)
{
    auto DVs = GradSVandermonde2D(N,r,s);
    // Create the differentiation operator of the reference element
    //-------------------------------------------------------------------------------------------------------
    return DVs*V2D.inverse();
    //-------------------------------------------------------------------------------------------------------
}

