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
////!!!!!!!!!!!!!!!!!!!!!!IMPLEMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bool    FourNodeQuadrilateralElement::CheckBoundaries(const std::vector<Node>& globalNodes) const{
    //check if face boundaries match node boundaries
    std::cout << "not yet implemented" << std::endl;
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
                outfq -= GAUSS_WEIGHTS.at(j) * mBCs.at(i).mX * N(xi, eta, outN) * 2.0 / d;
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

std::istream& operator >> (std::istream& stream, FourNodeQuadrilateralElement& element){
    for (auto& node : element.mNodes)
        stream >> node;
    for (unsigned int i = 0; i < element.mBCs.size(); i++)
        stream >> element.mBCs.at(i);

    return stream;
}
