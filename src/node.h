#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED

#include    "boundary_condition.h"

class   Node{
public:
    double  mX;
    double  mY;
    BC      mBC;
    int     mPosition;  ///< storage for reordered position

    double  NodeDistance(const Node& nd2) const;
    };



std::istream& operator >> (std::istream& stream, Node& node);

#endif // NODE_H_INCLUDED
