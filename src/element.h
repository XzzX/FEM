#ifndef ELEMENT_H_INCLUDED
#define ELEMENT_H_INCLUDED

#include    <vector>
#include    <Eigen/Dense>

class   Node;

class   FourNodeQuadrilateralElement{
public:
    typedef Eigen::Matrix<double, 2, 2> typeJ;
    typedef Eigen::Matrix<double, 2, 4> typeB;
    typedef Eigen::Matrix<double, 4, 4> typeK;

    std::vector<int32_t>    mNodes; ///< local->global nodes conversion table
    Eigen::MatrixXd mD; ///< conductivity matrix

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    FourNodeQuadrilateralElement();

    inline  static   double  N1(const double xi, const double eta) {
        return  0.25 * (1-xi) * (1-eta);
    }
    inline  static   double  N2(const double xi, const double eta) {
        return  0.25 * (1+xi) * (1-eta);
    }
    inline  static   double  N3(const double xi, const double eta) {
        return  0.25 * (1+xi) * (1+eta);
    }
    inline  static   double  N4(const double xi, const double eta) {
        return  0.25 * (1-xi) * (1+eta);
    }

    inline  static   double  dN1xi(const double xi, const double eta){
        return  0.25 * (-1+eta);
    }
    inline  static   double  dN2xi(const double xi, const double eta){
        return  0.25 * (1-eta);
    }
    inline  static   double  dN3xi(const double xi, const double eta){
        return  0.25 * (1+eta);
    }
    inline  static   double  dN4xi(const double xi, const double eta){
        return  0.25 * (-1-eta);
    }

    inline  static   double  dN1eta(const double xi, const double eta){
        return  0.25 * (-1+xi);
    }
    inline  static   double  dN2eta(const double xi, const double eta){
        return  0.25 * (-1-xi);
    }
    inline  static   double  dN3eta(const double xi, const double eta){
        return  0.25 * (1+xi);
    }
    inline  static   double  dN4eta(const double xi, const double eta){
        return  0.25 * (1-xi);
    }

    Eigen::Matrix<double, 2, 2>&    J(const std::vector<Node>& globalNodes, const double xi, const double eta, Eigen::Matrix<double, 2, 2>& outJ) const ;
    Eigen::Matrix<double, 2, 4>&    B(const std::vector<Node>& globalNodes, const double xi, const double eta, Eigen::Matrix<double, 2, 4>& outB) const ;
    Eigen::Matrix<double, 4, 4>&    K(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 4>& outK) const ;
};

std::istream& operator >> (std::istream& stream, FourNodeQuadrilateralElement& element);

#endif // ELEMENT_H_INCLUDED
