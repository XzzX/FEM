#include    "node.h"

#include    <cmath>

double  Node::NodeDistance(const Node& nd2) const {
    return sqrt((mX - nd2.mX)*(mX - nd2.mX) + (mY - nd2.mY)*(mY - nd2.mY));
}

std::istream& operator >> (std::istream& stream, Node& node){
    stream >> node.mX >> node.mY >> node.mBC.mType>> node.mBC.mX >> node.mBC.mY;
    return  stream;
}
