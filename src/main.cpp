#include    <fstream>
#include    <iostream>
#include    <vector>

#include    "node.h"
#include    "element.h"

int main(){
    std::vector<Node>   globalNodes;
    std::vector<FourNodeQuadrilateralElement>    elements;

    std::cout << "reading nodes ... ";
    Node node;
    std::ifstream    iFile("nodes.txt");
    if (!iFile) std::cerr << "unknown node file" << std::endl;
    while   (!iFile.eof()){
        iFile >> node.x >> node.y >> node.bc >> node.bcX >> node.bcY;
        globalNodes.push_back(node);
    }
    iFile.close();
    std::cout << globalNodes.size() << " read" << std::endl;

    std::cout << "reading elements ... ";
    FourNodeQuadrilateralElement element;
    iFile.open("elements.txt");
    if (!iFile) std::cerr << "unknown element file" << std::endl;
    while   (!iFile.eof()){
        iFile >> element;
        elements.push_back(element);
    }
    iFile.close();
    std::cout << elements.size() << " read" << std::endl;

    FourNodeQuadrilateralElement::typeK    outK;
    elements.at(0).K(globalNodes, outK);

    std::cout << outK << std::endl;

    return 0;
}
