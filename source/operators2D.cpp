/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "functions1D.h"
#include "functions2D.h"
#include "operators2D.h"

#define EIGEN_MPL2_ONLY

using namespace std;
using namespace Eigen;
using namespace operators2D;

ArrayXd operators2D::DDr(int N, ArrayXd& U)
{
    ArrayXd ret = VectorXd::Zero((N+1)*(N+1));

    for(int i=0; i<N+1; ++i)
    {
        EdgeMap A(ret.data()+i,N+1,InnerStride<>(N+1));
        EdgeMap B(U.data()+i,N+1,InnerStride<>(N+1));

        A = globals::Dr*B;
    }

    return ret;
}

ArrayXd operators2D::DDs(int N, ArrayXd& U)
{
    ArrayXd ret = VectorXd::Zero((N+1)*(N+1));

    for(int i=0; i<N+1; ++i)
    {
        EdgeMap A(ret.data()+i*(N+1),N+1,InnerStride<>(1));
        EdgeMap B(U.data()+i*(N+1),N+1,InnerStride<>(1));

        A = globals::Dr*B;
    }

    return ret;
}
