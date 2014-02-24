#include    <iostream>
#include    "node.h"

std::istream& operator >> (std::istream& stream, BC& bc){
    char    ch;
    stream >> ch;
    if (ch=='D')
        bc = BC::D;
    else if (ch=='N')
        bc = BC::N;
    else{
        std::cerr << "unknown boundary condition: " << ch << std::endl;
        std::cerr << "defaulting to N" << std::endl;
    }

    return stream;
}
