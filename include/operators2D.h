/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * Operators than can be applied to an 2D element data
 *
 * N : refers to the order of the polynomial in first dimension
 
 */

#include <Eigen/Dense>
#include "globals.h"

#define EIGEN_MPL2_ONLY

#ifndef OPERATORS2D_H
#define OPERATORS2D_H

using namespace Eigen;
using namespace std;

typedef Map<VectorXd,0,InnerStride<>> EdgeMap;

namespace operators2D{

// Partial derivative in r direction for rectangular element
ArrayXd DDr(int N, ArrayXd& U);

// Partial derivative in s direction for rectangular element
ArrayXd DDs(int N, ArrayXd& U);
}
#endif
