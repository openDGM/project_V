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

using namespace std;
using namespace Eigen;
using namespace functions1D;
using namespace functions2D;

//===========================================================================================================
// Due to expensive sparse matrix creation operation in Eigen these functions are not actively used at the
// moment but will remain here as they might be useful at later stages.
//===========================================================================================================

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
    SparseMatrix<double> ret(S,S);
    ret.reserve((N+1)*S);

    // Find non-zero entries
    for (int i=0; i<N+1; ++i)
    {
        for (int j=0; j<N+1; ++j)
        {
            for (int k=0; k<N+1; ++k)
            {
                ret.insert(i*(N+1)+k, j*(N+1)+k) = Dr(i*(N+1)+k,j*(N+1)+k);
            }
        }
    }
 
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
    SparseMatrix<double> ret(S,S);
    ret.reserve((N+1)*S);

    // Find non-zero entries
    for (int k=0; k<N+1; ++k)
    {
        for (int i=k*(N+1); i<(k+1)*(N+1); ++i)
        {
            for (int j=k*(N+1); j<(k+1)*(N+1); ++j)
            {
                ret.insert(i,j) = Ds(i,j);
            }
        }
    }

    return ret;
    //-------------------------------------------------------------------------------------------------------
}

VectorXd functions2D::LocalX2D(const ArrayXd& r, pair<double,double> BottomLeft, pair<double,double> TopRight)
{
    int N = r.size();
    VectorXd x = LocalX1D(r,BottomLeft.first,TopRight.first);
    VectorXd ret(N*N);

    for (int i=0;i<N;++i)
    {
        for (int j=0;j<N;++j)
        {
            ret[i*N+j] = x[i];
        }
    }

    return ret;
}

VectorXd functions2D::LocalY2D(const ArrayXd& s, pair<double,double> BottomLeft, pair<double,double> TopRight)
{
    int N = s.size();
    VectorXd y = LocalX1D(s,BottomLeft.second,TopRight.second);
    VectorXd ret(N*N);

    for (int i=0;i<N;++i)
    {
        for (int j=0;j<N;++j)
        {
            ret[i*N+j] = y[j];
        }
    }

    return ret;
}

SparseMatrix<double> functions2D::Lift2D(int N)
{
    ArrayXd r    = jacobiGL(0,0,N);
    MatrixXd V1D = Vandermonde1D(N,r);
    MatrixXd V2D = Vandermonde2D(N,r,r);
    MatrixXd M   = (V1D*V1D.transpose()).inverse();

    // Create operator for 2D surface flux integral on the reference element
    //-------------------------------------------------------------------------------------------------------
    MatrixXd Emat = MatrixXd::Zero((N+1)*(N+1),4*(N+1));


    Map<MatrixXd,0,Stride<Dynamic,Dynamic>> Left(Emat.data(),N+1,N+1,Stride<Dynamic,Dynamic>((N+1)*(N+1),1));
    Map<MatrixXd,0,Stride<Dynamic,Dynamic>> Top(Emat.data()+(N+1)*(N+1)*(N+1)+N,N+1,N+1,Stride<Dynamic,Dynamic>((N+1)*(N+1),N+1));
    Map<MatrixXd,0,Stride<Dynamic,Dynamic>> Right(Emat.data()+2*(N+1)*(N+1)*(N+1)+N*(N+1),N+1,N+1,Stride<Dynamic,Dynamic>((N+1)*(N+1),1));
    Map<MatrixXd,0,Stride<Dynamic,Dynamic>> Bottom(Emat.data()+3*(N+1)*(N+1)*(N+1),N+1,N+1,Stride<Dynamic,Dynamic>((N+1)*(N+1),N+1));

    Left = M;
    Top = M;
    Right = M;
    Bottom = M;

    int R = (N+1)*(N+1);
    int S = 4*(N+1);
    auto tmp = V2D*(V2D.transpose()*Emat);

    // Make operator sparse
    //-------------------------------------------------------------------------------------------------------
    SparseMatrix<double> ret(R,S);
    ret.reserve(4*R);

    // Find non-zero entries
    for (int i=0; i<N+1; ++i)
    {
        for (int k=0; k<N+1; ++k)
        {
            ret.insert(i*(N+1)+k,k) = tmp(i*(N+1)+k,k);
            ret.insert(i*(N+1)+k,i+N+1) = tmp(i*(N+1)+k,i+N+1);
            ret.insert(i*(N+1)+k,k+2*(N+1)) = tmp(i*(N+1)+k,k+2*(N+1));
            ret.insert(i*(N+1)+k,i+3*(N+1)) = tmp(i*(N+1)+k,i+3*(N+1));
        }
    }

    return ret;
    //-------------------------------------------------------------------------------------------------------
}
