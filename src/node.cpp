#include    "node.h"

#include    <cmath>

double  Node::NodeDistance(const Node& nd2) const {
    return sqrt((mX - nd2.mX)*(mX - nd2.mX) + (mY - nd2.mY)*(mY - nd2.mY));
}

bool	Node::IsEqual(const Node& nd2) const {
	const double	eps = 1e-7;
	if (fabs(mX - nd2.mX) > eps) return false;
	if (fabs(mY - nd2.mY) > eps) return false;
	if (mBC.mType!=nd2.mBC.mType) return false;
	if (fabs(mBC.mX - nd2.mBC.mX) > eps) return false;
	if (fabs(mBC.mX - nd2.mBC.mX) > eps) return false;
	return true;
}

std::istream& operator >> (std::istream& stream, Node& node){
    stream >> node.mX >> node.mY >> node.mBC.mType>> node.mBC.mX >> node.mBC.mY;
    return  stream;
}
