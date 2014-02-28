#ifndef BOUNDARY_CONDITION_H_INCLUDED
#define BOUNDARY_CONDITION_H_INCLUDED

#include    <sstream>

enum class BCType {N, D};
std::istream& operator >> (std::istream& stream, BCType& type);

class   BC{
public:
    static const BC DEFAULT;
    BCType  mType = BCType::N;
    double  mX = 0;
    double  mY = 0;
};
std::istream& operator >> (std::istream& stream, BC& bc);

#endif // BOUNDARY_CONDITION_H_INCLUDED
