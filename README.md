# project_V
This has been a small learning project of mine
to study the Nodal Discontinuous Galerking scheme
to solve initial-boundary value problems. The
implementation follows chapter 3 of the book
Nodal Discontinuous Galerkin Methods Algorithms, 
Analysis, and Applications by Hesthaven and 
Warburton. It's a simple advection problem with
periodic boundaries.

Requirements:  
c++11 compiler,  
Eigen library (eigen.tuxfamily.org/)

Set path to Eigen in Makefile and type make.
An executable will be compiled to build/bin when
you run the program you will be asked to input
number of Nodal elements, Polynom order of the 
Element, advection speed and timestep.

The domain is set to the interval 0,2PI and time
as well from 0,2PI. You can estimate the timestep
dt < 1/(2N+1)dx/|a|
where N is the order of the element, dx is width
of the element and a is advection speed. E.g. for
N=4 and number of elements=4 and a=1 the time step
would need to be dt < PI/18 ~= 0.15s

The program will write an output to "output.dat"
with the first line containing the coordinates and 
the second line containing the function values at 
the final time step.
