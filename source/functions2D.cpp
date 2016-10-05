/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "functions1D.h"
#include "functions2D.h"
#define EIGEN_MPL2_ONLY
#define epsilon 1e-10 

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace functions1D;
using namespace functions2D;

ArrayXd functions2D::jacobiP2D(const ArrayXd& r, const ArrayXd& s, int alpha, int beta, int N, int M)
{
    // 2D Jacobi polynomial constructed from orthorgonal 1D Jacobi polynomials
    ArrayXd jacobi2D(r.size()*s.size());

    // 1D polynomials in r and s direction
    ArrayXd jacobiPr = jacobiP(r,alpha,beta,N);
    ArrayXd jacobiPs = jacobiP(s,alpha,beta,M);

    // Tensor product of the 1D polynomials construct 2D polynomial.
    for (int i=0; i<r.size(); ++i)
    {
        for(int j=0; j<s.size(); ++j)
        {
            jacobi2D[i*(r.size())+j] = jacobiPr[i]*jacobiPs[j];
        }
    }

    return jacobi2D;
}

VectorXd functions2D::GradRJacobiP2D(const VectorXd& r, const ArrayXd& s, int alpha, int beta, int N, int M)
{
    // N Zero order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    if (N==0)
    {
        return VectorXd::Zero(r.size()*s.size());
    }
    //-------------------------------------------------------------------------------------------------------

    // N,Mth order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    ArrayXd jacobi2DGradR(r.size()*s.size());

    ArrayXd jacobiPGradR = sqrt(N*(N+alpha+beta+1))*jacobiP(r,alpha+1,beta+1,N-1);
    ArrayXd jacobiPS = jacobiP(s,alpha,beta,M);

    // Tensor product of gradient polynomial in r and polynomial in s form 2D gradient in r of order N,M
    // polynomial 
    for (int i=0; i<r.size(); ++i)
    {
        for (int j=0; j<s.size(); ++j)
        {
            jacobi2DGradR[i*(r.size())+j] = jacobiPGradR[i]*jacobiPS[j];
        }
    }

    return jacobi2DGradR;
    //-------------------------------------------------------------------------------------------------------
}

VectorXd functions2D::GradSJacobiP2D(const VectorXd& r, const ArrayXd& s, int alpha, int beta, int N, int M)
{
    // M Zero order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    if (M==0)
    {
        return VectorXd::Zero(r.size()*s.size());
    }
    //-------------------------------------------------------------------------------------------------------

    // N,Mth order jacobi polynomial
    //-------------------------------------------------------------------------------------------------------
    ArrayXd jacobi2DGradS(r.size()*s.size());

    ArrayXd jacobiPGradS = sqrt(M*(M+alpha+beta+1))*jacobiP(s,alpha+1,beta+1,M-1);
    ArrayXd jacobiPR = jacobiP(r,alpha,beta,N);


    // Tensor product of gradient polynomial in r and polynomial in s form 2D gradient in r of order N,M
    // polynomial 
    for (int i=0; i<r.size(); ++i)
    {
        for (int j=0; j<s.size(); ++j)
        {
            jacobi2DGradS[i*(r.size())+j] = jacobiPR[i]*jacobiPGradS[j];
        }
    }

    return jacobi2DGradS;
    //-------------------------------------------------------------------------------------------------------
}

MatrixXd functions2D::Vandermonde2D(int N, const VectorXd& r, const VectorXd& s)
{

    MatrixXd V2D(r.size()*s.size(), (N+1)*(N+1));
    // Create Vandermonde matrix that contains the base functions of the reference element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            V2D.col(i*(N+1)+j) = jacobiP2D(r,s,0,0,i,j);
        }
    }
    //-------------------------------------------------------------------------------------------------------
    return V2D;
}

MatrixXd functions2D::GradRVandermonde2D(int N, const VectorXd& r, const VectorXd& s)
{
    MatrixXd DVr(r.size()*s.size(),(N+1)*(N+1));
    // Create the GradVandermonde matrix that contains the gradients of base functions of the reference 
    // element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            DVr.col(i*(N+1)+j) = GradRJacobiP2D(r,s,0,0,i,j);
        }
    }
    //-------------------------------------------------------------------------------------------------------
    return DVr;
}

MatrixXd functions2D::GradSVandermonde2D(int N, const VectorXd& r, const VectorXd& s)
{
    MatrixXd DVs(r.size()*s.size(),(N+1)*(N+1));
    // Create the GradVandermonde matrix that contains the gradients of base functions of the reference 
    // element as columns
    //-------------------------------------------------------------------------------------------------------
    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            DVs.col(i*(N+1)+j) = GradSJacobiP2D(r,s,0,0,i,j);
        }
    }
    //-------------------------------------------------------------------------------------------------------
    return DVs;
}

SparseMatrix<double> functions2D::DRMatrix2D(int N, const VectorXd& r, const VectorXd& s)
{
    int S = (N+1)*(N+1);
  
    // Create the differentiation operator of the reference element
    auto DVr = GradRVandermonde2D(N,r,s);
    auto V2D = Vandermonde2D(N,r,s);
    auto Dr = DVr*V2D.inverse();

    // Make operator sparse
    //-------------------------------------------------------------------------------------------------------
    vector<Triplet<double>> triplets;
    SparseMatrix<double> ret(S,S);

    // Find non-zero entries
    for (int i=0; i<S; ++i)
    {
        for (int j=0; j<S; ++j)
        {
            if(abs(Dr(i,j)) > epsilon) triplets.push_back(Triplet<double>(i,j,Dr(i,j)));
        }
    }
 
    // Set sparse matrix values
    ret.setFromTriplets(triplets.begin(),triplets.end());

    return ret;
    //-------------------------------------------------------------------------------------------------------
}

SparseMatrix<double> functions2D::DSMatrix2D(int N, const VectorXd& r, const VectorXd& s)
{
    int S = (N+1)*(N+1);

    // Create the differentiation operator of the reference element
    auto DVs = GradSVandermonde2D(N,r,s);
    auto V2D = Vandermonde2D(N,r,s);
    auto Ds = DVs*V2D.inverse();

    // Make operator sparse
    //-------------------------------------------------------------------------------------------------------
    vector<Triplet<double>> triplets;
    SparseMatrix<double> ret(S,S);

    // Find non-zero entries
    for (int i=0; i<S; ++i)
    {
        for (int j=0; j<S; ++j)
        {
            if(abs(Ds(i,j)) > epsilon) triplets.push_back(Triplet<double>(i,j,Ds(i,j)));
        }
    }

    // Set sparse matrix values
    ret.setFromTriplets(triplets.begin(),triplets.end());

    return ret;
    //-------------------------------------------------------------------------------------------------------
}
