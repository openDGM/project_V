/*
 *  Copyright (C) Andreas Tack, 2016
 *  This Source Code is subject to the 
 *  terms of the Mozilla Public License, 
 *  v. 2.0.
 */

#include "edge.h"
#include "element.h"
#include "globals.h"
#include "functions1D.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

#define PI 3.14159265

using namespace std;
using namespace Eigen;
using namespace functions1D;

  //declare globals
  MatrixXd globals::V;
  MatrixXd globals::invV;
  MatrixXd globals::DVr;
  MatrixXd globals::Dr;
  MatrixXd globals::Lift;

double AskInput(string s)
{
  double value;
  cout << s << endl;
  cin >> value;
  return value;
}

int main()
{

  // output file
  ofstream ofile;
  ofile.open ("output.dat");

  // User input
  //---------------------------------------------------------------------------------------------------------
  int K       = AskInput("Number of elements: ");
  int N       = AskInput("Polynomial order of elements: ");
  double a    = AskInput("Enter advection speed: ");
  double dt   = AskInput("Enter timestep: ");
  //---------------------------------------------------------------------------------------------------------
  
  // Initial time
  double time = 0;

  // Init globals
  //---------------------------------------------------------------------------------------------------------
  globals::V     = Vandermonde1D(N,jacobiGL(0,0,N));
  globals::invV  = globals::V.inverse();
  globals::DVr   = GradVandermonde1D(N,jacobiGL(0,0,N));
  globals::Lift  = Lift1D(N,globals::V);
  globals::Dr    = DMatrix1D(N, jacobiGL(0,0,N), globals::V);
  //---------------------------------------------------------------------------------------------------------

  // Set up a DGM grid consisting of K elements of order N
  vector<element> theElements(K);
  vector<edge>    theEdges(K);

  for (int i=0; i<K; ++i)
  {
    theElements[i] = element(i*2*PI/K,(i+1)*2*PI/K,-1,1,N);
    // write grid to output file
    ofile << theElements[i].getX().transpose() << " ";
  }
  ofile << endl;

  for (int i=0; i<K; ++i)
  {
    // periodic boundary
    if (i==0)
    {
      theEdges[0] = edge(theElements[K-1].connectToRightBound(), theElements[0].connectToLeftBound(),1);
    }
    else
    {
      theEdges[i] = edge(theElements[i-1].connectToRightBound(), theElements[i].connectToLeftBound(),1);
    }
  }

  for (int i=0; i<K; ++i)
  {
    // periodic boundary
    if (i==K-1)
    {
      theElements[K-1].setLEdge(&theEdges[K-1]);
      theElements[K-1].setREdge(&theEdges[0]);
    }
    else
    {
      theElements[i].setLEdge(&theEdges[i]);
      theElements[i].setREdge(&theEdges[i+1]);
    }
  }
  //---------------------------------------------------------------------------------------------------------

  // Calculate initial condition U=sin(x)
  //---------------------------------------------------------------------------------------------------------
  for (int i=0; i<K; ++i)
  {
    theElements[i].setU(sin(theElements[i].getX().array()));
  }
  //---------------------------------------------------------------------------------------------------------

  // Integrate solution in time
  //---------------------------------------------------------------------------------------------------------
  while (time < 2*PI)
  {
    for (int INTRK=0; INTRK<5; ++INTRK)
    {
      // Update fluxes
      for (int i=0; i<K; i++)
      {
        theEdges[i].evaluateFlux(a);
      }

      // Perform RK step
      for (int i=0; i<K; i++)
      {
        theElements[i].updateFluxes();
        theElements[i].advecRHS(a);
        theElements[i].advecRK1D(INTRK,dt);
      }
    }
  time += dt;
  }
  //---------------------------------------------------------------------------------------------------------

  // Output solution to output file
  //---------------------------------------------------------------------------------------------------------
  for (int i=0; i<K-1; i++)
  {
    ofile << theElements[i].getU().transpose() << " ";
  }
  ofile << theElements[K-1].getU().transpose() << endl;
  //---------------------------------------------------------------------------------------------------------

  // Close output file
  ofile.close();

  return 0;
}
