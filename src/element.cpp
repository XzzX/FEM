#include    "element.h"

#include    <cmath>
#include    <iostream>
#include    <vector>

#include    "node.h"
#include    "gauss_quadrature.h"
#include    "boundary_condition.h"

FourNodeQuadrilateralElement::FourNodeQuadrilateralElement() : mNodes(4,-1), mBCs(4), mD(2,2), mS(6) {
    mD.setZero();
    mD(0,0) = 5;
    mD(1,1) = 5;
}

FourNodeQuadrilateralElement& FourNodeQuadrilateralElement::operator=(const FourNodeQuadrilateralElement& rhs) {
	// Now, swap the data members with the temporary:
	mNodes = rhs.mNodes;
	mBCs = rhs.mBCs;
	mD(0, 0) = rhs.mD(0, 0);
	mD(0, 1) = rhs.mD(0, 1);
	mD(1, 0) = rhs.mD(1, 0);
	mD(1, 1) = rhs.mD(1, 1);
	mS = rhs.mS;

	return *this;
}

////!!!!!!!!!!!!!!!!!!!!!!IMPLEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bool    FourNodeQuadrilateralElement::CheckBoundaries(const std::vector<Node>& globalNodes) const{
    //check if face boundaries match node boundaries
    std::cout << "not yet implemented" << std::endl;
	return true;
}

Eigen::Matrix<double, 4, 1>&    FourNodeQuadrilateralElement::N(const double xi, const double eta, Eigen::Matrix<double, 4, 1>& outN) const {
    outN(0,0) = N1(xi, eta);
    outN(1,0) = N2(xi, eta);
    outN(2,0) = N3(xi, eta);
    outN(3,0) = N4(xi, eta);
    return outN;
}

Eigen::Matrix<double, 2, 2>&    FourNodeQuadrilateralElement::J(const std::vector<Node>& globalNodes, const double xi, const double eta, Eigen::Matrix<double, 2, 2>& outJ) const {
    Eigen::Matrix<double, 2, 4> G;
    G <<    dN1xi(xi, eta), dN2xi(xi, eta), dN3xi(xi, eta), dN4xi(xi, eta),
            dN1eta(xi, eta), dN2eta(xi, eta), dN3eta(xi, eta), dN4eta(xi, eta);

    Eigen::Matrix<double, 4, 2> X;
    X <<    globalNodes.at(mNodes.at(0)).mX, globalNodes.at(mNodes.at(0)).mY,
            globalNodes.at(mNodes.at(1)).mX, globalNodes.at(mNodes.at(1)).mY,
            globalNodes.at(mNodes.at(2)).mX, globalNodes.at(mNodes.at(2)).mY,
            globalNodes.at(mNodes.at(3)).mX, globalNodes.at(mNodes.at(3)).mY;

    outJ = G*X;
    return outJ;
}

Eigen::Matrix<double, 2, 4>&    FourNodeQuadrilateralElement::B(const std::vector<Node>& globalNodes, const double xi, const double eta, Eigen::Matrix<double, 2, 4>& outB) const {
    Eigen::Matrix<double, 2, 4> G;
    G <<    dN1xi(xi, eta), dN2xi(xi, eta), dN3xi(xi, eta), dN4xi(xi, eta),
            dN1eta(xi, eta), dN2eta(xi, eta), dN3eta(xi, eta), dN4eta(xi, eta);

    Eigen::Matrix<double, 4, 2> X;
    X <<    globalNodes.at(mNodes.at(0)).mX, globalNodes.at(mNodes.at(0)).mY,
            globalNodes.at(mNodes.at(1)).mX, globalNodes.at(mNodes.at(1)).mY,
            globalNodes.at(mNodes.at(2)).mX, globalNodes.at(mNodes.at(2)).mY,
            globalNodes.at(mNodes.at(3)).mX, globalNodes.at(mNodes.at(3)).mY;

    outB = (G*X).inverse()*G;
    return outB;
}

Eigen::Matrix<double, 4, 4>&    FourNodeQuadrilateralElement::K(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 4>& outK) const {
    Eigen::Matrix<double, 2, 4> outB;
    Eigen::Matrix<double, 2, 2> outJ;

    outK.setZero();

    for (unsigned int i=0; i<GAUSS_POINTS.size(); i++)
        for (unsigned int j=0; j<GAUSS_POINTS.size(); j++){
            B(globalNodes, GAUSS_POINTS.at(i), GAUSS_POINTS.at(j), outB);
            J(globalNodes, GAUSS_POINTS.at(i), GAUSS_POINTS.at(j), outJ);
            outK += GAUSS_WEIGHTS.at(i) * GAUSS_WEIGHTS.at(j) * outB.transpose() * mD * outB * fabs(outJ.determinant());
        }

    return outK;
}

Eigen::Matrix<double, 4, 1>&    FourNodeQuadrilateralElement::fs(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 1>& outfs) const {
    Eigen::Matrix<double, 4, 1> outN;
    Eigen::Matrix<double, 2, 2> outJ;

    outfs.setZero();

    for (unsigned int i=0; i<GAUSS_POINTS.size(); i++)
        for (unsigned int j=0; j<GAUSS_POINTS.size(); j++){
            J(globalNodes, GAUSS_POINTS.at(i), GAUSS_POINTS.at(j), outJ);
            N(GAUSS_POINTS.at(i), GAUSS_POINTS.at(j), outN);
            outfs += GAUSS_WEIGHTS.at(i) * GAUSS_WEIGHTS.at(j) * outN * mS * fabs(outJ.determinant());
        }

    return outfs;
}

Eigen::Matrix<double, 4, 1>&    FourNodeQuadrilateralElement::fq(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 1>& outfq) const {
    double  xi = 0;
    double  eta = 0;
    double  d = 0;

    Eigen::Matrix<double, 4, 1> outN;

    outfq.setZero();
    for (unsigned int i=0; i<mBCs.size(); i++){
        if (mBCs.at(i).mType == BCType::N){
            for (unsigned int j=0; j<GAUSS_POINTS.size(); j++){
                switch (i){
                    case    0: eta = -1; xi = GAUSS_POINTS.at(j);
                            d = globalNodes.at(mNodes.at(0)).NodeDistance(globalNodes.at(mNodes.at(1))); break;
                    case    1: eta = GAUSS_POINTS.at(j); xi = 1;
                            d = globalNodes.at(mNodes.at(1)).NodeDistance(globalNodes.at(mNodes.at(2))); break;
                    case    2: eta =  1; xi = GAUSS_POINTS.at(j);
                            d = globalNodes.at(mNodes.at(2)).NodeDistance(globalNodes.at(mNodes.at(3))); break;
                    case    3: eta = GAUSS_POINTS.at(j); xi = -1;
                            d = globalNodes.at(mNodes.at(3)).NodeDistance(globalNodes.at(mNodes.at(0))); break;
                }
                outfq -= GAUSS_WEIGHTS.at(j) * mBCs.at(i).mX * N(xi, eta, outN) * d * 0.5;
            }
        }
    }

    for (unsigned int i=0; i<mNodes.size(); i++){
        if (globalNodes.at(mNodes.at(i)).mBC.mType == BCType::N){
            outfq(i,0) -= globalNodes.at(mNodes.at(i)).mBC.mX;
        }
    }

    return outfq;
}

Eigen::Matrix<double, 4, 1>&    FourNodeQuadrilateralElement::f(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 1>& outf) const {
    Eigen::Matrix<double, 4, 1> outfq;
    fq(globalNodes, outfq);
    fs(globalNodes, outf);
    outf += outfq;

    return outf;
}

void    FourNodeQuadrilateralElement::Refine(std::vector<Node>& globalNodes, std::vector<FourNodeQuadrilateralElement>& elements){
	std::vector<int>	newNodes(5, -1);
    Node    node;
    FourNodeQuadrilateralElement    element;

    for (int i=0; i<4; i++){
        node.mX = (globalNodes.at(mNodes.at(i)).mX + globalNodes.at(mNodes.at((i+1)%4)).mX) * 0.5;
        node.mY = (globalNodes.at(mNodes.at(i)).mY + globalNodes.at(mNodes.at((i+1)%4)).mY) * 0.5;
        node.mBC = mBCs.at(i);
        if (node.mBC.mType == BCType::N){
            node.mBC.mX = 0; node.mBC.mY = 0;
        }
		for (unsigned int j = 0; j < globalNodes.size(); j++) {
			if (node.IsEqual(globalNodes.at(j))) {
				newNodes.at(i) = j;
				break;
			}
		}
		if (newNodes.at(i) == -1) {
			globalNodes.push_back(node);
			newNodes.at(i) = globalNodes.size() - 1;
		}
    }
    node.mX = 0; node.mY = 0;
    for (int i=0; i<4; i++){
        node.mX += globalNodes.at(mNodes.at(i)).mX;
        node.mY += globalNodes.at(mNodes.at(i)).mY;
    }
    node.mX *= 0.25; node.mY *= 0.25;
    node.mBC = BC::DEFAULT;

	for (unsigned int j = 0; j < globalNodes.size(); j++) {
		if (node.IsEqual(globalNodes.at(j))) {
			newNodes.at(4) = j;
			break;
		}
	}
	if (newNodes.at(4) == -1) {
		globalNodes.push_back(node);
		newNodes.at(4) = globalNodes.size() - 1;
	}

    element = *this;
	element.mNodes.at(1) = newNodes.at(0);
	element.mNodes.at(2) = newNodes.at(4);
	element.mNodes.at(3) = newNodes.at(3);
    element.mBCs.at(1) = BC::DEFAULT;
    element.mBCs.at(2) = BC::DEFAULT;
	elements.push_back(element);

	element = *this;
	element.mNodes.at(0) = newNodes.at(0);
	element.mNodes.at(2) = newNodes.at(1);
	element.mNodes.at(3) = newNodes.at(4);
    element.mBCs.at(2) = BC::DEFAULT;
    element.mBCs.at(3) = BC::DEFAULT;
	elements.push_back(element);

    element = *this;
	element.mNodes.at(0) = newNodes.at(4);
	element.mNodes.at(1) = newNodes.at(1);
	element.mNodes.at(3) = newNodes.at(2);
    element.mBCs.at(0) = BC::DEFAULT;
    element.mBCs.at(3) = BC::DEFAULT;
	elements.push_back(element);

	mNodes.at(0) = newNodes.at(3);
	mNodes.at(1) = newNodes.at(4);
	mNodes.at(2) = newNodes.at(2);
    mBCs.at(0) = BC::DEFAULT;
    mBCs.at(1) = BC::DEFAULT;
}

std::istream& operator >> (std::istream& stream, FourNodeQuadrilateralElement& element){
    for (auto& node : element.mNodes)
        stream >> node;
    for (unsigned int i = 0; i < element.mBCs.size(); i++)
        stream >> element.mBCs.at(i);

    return stream;
}
