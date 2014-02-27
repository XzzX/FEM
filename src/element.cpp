#include    "element.h"

#include    <cmath>
#include    <vector>

#include    "node.h"
#include    "gauss_quadrature.h"

FourNodeQuadrilateralElement::FourNodeQuadrilateralElement() : mNodes(4,-1), mD(2,2) {
    mD.setZero();
    mD(0,0) = 5;
    mD(1,1) = 5;
}

Eigen::Matrix<double, 2, 2>&    FourNodeQuadrilateralElement::J(const std::vector<Node>& globalNodes, const double xi, const double eta, Eigen::Matrix<double, 2, 2>& outJ) const {
    Eigen::Matrix<double, 2, 4> G;
    G <<    dN1xi(xi, eta), dN2xi(xi, eta), dN3xi(xi, eta), dN4xi(xi, eta),
            dN1eta(xi, eta), dN2eta(xi, eta), dN3eta(xi, eta), dN4eta(xi, eta);

    Eigen::Matrix<double, 4, 2> X;
    X <<    globalNodes.at(mNodes.at(0)).x, globalNodes.at(mNodes.at(0)).y,
            globalNodes.at(mNodes.at(1)).x, globalNodes.at(mNodes.at(1)).y,
            globalNodes.at(mNodes.at(2)).x, globalNodes.at(mNodes.at(2)).y,
            globalNodes.at(mNodes.at(3)).x, globalNodes.at(mNodes.at(3)).y;

    outJ = G*X;
    return outJ;
}

Eigen::Matrix<double, 2, 4>&    FourNodeQuadrilateralElement::B(const std::vector<Node>& globalNodes, const double xi, const double eta, Eigen::Matrix<double, 2, 4>& outB) const {
    Eigen::Matrix<double, 2, 4> G;
    G <<    dN1xi(xi, eta), dN2xi(xi, eta), dN3xi(xi, eta), dN4xi(xi, eta),
            dN1eta(xi, eta), dN2eta(xi, eta), dN3eta(xi, eta), dN4eta(xi, eta);

    Eigen::Matrix<double, 4, 2> X;
    X <<    globalNodes.at(mNodes.at(0)).x, globalNodes.at(mNodes.at(0)).y,
            globalNodes.at(mNodes.at(1)).x, globalNodes.at(mNodes.at(1)).y,
            globalNodes.at(mNodes.at(2)).x, globalNodes.at(mNodes.at(2)).y,
            globalNodes.at(mNodes.at(3)).x, globalNodes.at(mNodes.at(3)).y;

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

std::istream& operator >> (std::istream& stream, FourNodeQuadrilateralElement& element){
    for (unsigned int i=0; i<element.mNodes.size(); i++)
        stream >> element.mNodes.at(i);
    return stream;
}
