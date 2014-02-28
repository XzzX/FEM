#ifndef BOUNDARY_CONDITION_H_INCLUDED
#define BOUNDARY_CONDITION_H_INCLUDED

#include    <sstream>

enum class BCType {N, D};
std::istream& operator >> (std::istream& stream, BCType& type);

class   BC{
public:
    BCType  mType;
    double  mX;
    double  mY;
};
std::istream& operator >> (std::istream& stream, BC& bc);

#endif // BOUNDARY_CONDITION_H_INCLUDED
