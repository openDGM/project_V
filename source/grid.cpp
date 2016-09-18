#include "grid.h"

  grid::grid(int N)
  {
    itsGrid.resize(N);
  }

  grid::~grid()
  {

  }

  grid::buildNodes()
  {
    for (int i=0; i<itsGrid.size(); ++i)
      {
        itsGrid[i]=node(itsVertices[itsEdge[i][0]],itsVertices[itsEdge[i][1]]);
      }
  }
