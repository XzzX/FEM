#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include    <sstream>

enum class BC {N, D};
std::istream& operator >> (std::istream& stream, BC& bc);

struct   Node{
    double  x;
    double  y;
    BC      bc;
    double  bcX;
    double  bcY;
    };

#endif // NODE_H_INCLUDED
