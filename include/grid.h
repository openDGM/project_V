#include <vector>
#include "element"

class grid
{
    public:
    std::vector<double> itsVertices;
    std::vector<int[2]> itsEdges;
    std::vector<node> itsGrid;
    
    void readGrid(std::string filename);
    void buildNodes();
    void connectNodes(node& leftNode, node& rightNode)
}
