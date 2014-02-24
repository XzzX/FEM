#include    <fstream>
#include    <iostream>
#include    <vector>

#include    "node.h"
#include    "element.h"

int main(){
    std::vector<Node>   nodes;
    std::vector<Element>    elements;

    std::cout << "reading nodes ... ";
    Node node;
    std::ifstream    iFile("nodes.txt");
    if (!iFile) std::cerr << "unknown node file" << std::endl;
    while   (!iFile.eof()){
        iFile >> node.x >> node.y >> node.bc >> node.bcX >> node.bcY;
        nodes.push_back(node);
    }
    iFile.close();
    std::cout << nodes.size() << " read" << std::endl;

    std::cout << "reading elements ... ";
    Element element;
    iFile.open("elements.txt");
    if (!iFile) std::cerr << "unknown element file" << std::endl;
    while   (!iFile.eof()){
        iFile >> element.ll >> element.lr >> element.ur >> element.ul;
        elements.push_back(element);
    }
    iFile.close();
    std::cout << elements.size() << " read" << std::endl;

    return 0;
}
