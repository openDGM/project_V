#include <iostream>
#include <Eigen/Dense>
#include "jacobiP.h"

int main()
{
  VectorXd x(7);
  
  x << -1.0,-8.3022e-01,-4.6885e-01,0,4.6885e-01,8.3022e-01,1.0;

  auto y = jacobiGL(0,0,6);
  auto V1D = Vandermonde1D(6,y);

  std::cout << y << std::endl;
 /* std::cout << jacobiP(x,0,0,0) << std::endl<< std::endl;
  std::cout << jacobiP(x,0,0,1) << std::endl<< std::endl;
  std::cout << jacobiP(x,0,0,2) << std::endl<< std::endl;
  std::cout << jacobiP(x,0,0,3) << std::endl<< std::endl;
  std::cout << jacobiP(x,0,0,4) << std::endl<< std::endl;
*/
  return 0;
}
