/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

/*
 * Description:
 * The edge class is used to solve the
 * Riemann problem on the interface between
 * adjecent elements and return the flux.
 */

#include <Eigen/Dense>
#include <memory>
#define EIGEN_MPL2_ONLY

#ifndef EDGE_H
#define EDGE_H

using namespace std;
using namespace Eigen;

typedef Map<VectorXd,0,InnerStride<>> EdgeMap;

class edge
{
    private:
    double itsFlux;
    double* itsLValue;
    double* itsRValue;
  
    double itsNormal;

    public:
    edge();

    // Construct a new edge from pointers to left and right value and a edge normal
    edge(double* theLValue, double* theRValue, double theNormal);

    // Solve Riemann problem at the edge
    void evaluateFlux(const double a);

    // Return flux
    double Flux(double theNormal);
};

class edge2D
{
    private:
    MatrixXd itsFlux2D;
    unique_ptr<EdgeMap> itsLValue;
    unique_ptr<EdgeMap> itsRValue;

    Vector2d itsLNormal;
    Vector2d itsRNormal;

    public:
    edge2D();

    // Construct a new edge from map to left and right value and a edge normal
    edge2D(unique_ptr<EdgeMap> theLValue, unique_ptr<EdgeMap> theRValue, Vector2d theLNormal, Vector2d theRNormal);

    // Solve Riemann problem at the edge
    void evaluateFlux2D(const Vector2d a);

    // Return flux
    MatrixXd Flux2D();
};
#endif 
