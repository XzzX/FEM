#include    <iostream>
#include    "boundary_condition.h"

const   BC  BC::DEFAULT;

std::istream& operator >> (std::istream& stream, BCType& type){
    char    ch;
    stream >> ch;
    if (ch=='D')
        type = BCType::D;
    else if (ch=='N')
        type = BCType::N;
    else{
        std::cerr << "unknown boundary condition: " << ch << std::endl;
        std::cerr << "defaulting to N" << std::endl;
        type = BCType::N;
    }

    return stream;
}

std::istream& operator >> (std::istream& stream, BC& bc){
    stream >> bc.mType >> bc.mX >> bc.mY;
    return stream;
}
