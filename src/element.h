#ifndef ELEMENT_H_INCLUDED
#define ELEMENT_H_INCLUDED

#include    <vector>
#include    <Eigen/Dense>

class   Node;
class   BC;

class   FourNodeQuadrilateralElement{
private:
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
public:
    typedef Eigen::Matrix<double, 4, 1> typeN;
    typedef Eigen::Matrix<double, 2, 2> typeJ;
    typedef Eigen::Matrix<double, 2, 4> typeB;
    typedef Eigen::Matrix<double, 4, 4> typeK;

    std::vector<int32_t>    mNodes; ///< local->global nodes conversion table
    std::vector<BC>         mBCs;
    Eigen::MatrixXd         mD; ///< conductivity matrix
    double                  mS; ///< global source

    FourNodeQuadrilateralElement();

    typeN&    N(const double xi, const double eta, typeN& outN) const ;
    typeJ&    J(const std::vector<Node>& globalNodes, const double xi, const double eta, typeJ& outJ) const ;
    typeB&    B(const std::vector<Node>& globalNodes, const double xi, const double eta, typeB& outB) const ;
    typeK&    K(const std::vector<Node>& globalNodes, typeK& outK) const ;

    Eigen::Matrix<double, 4, 1>&    fs(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 1>& outfs) const ;
    Eigen::Matrix<double, 4, 1>&    fq(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 1>& outfg) const ;
    Eigen::Matrix<double, 4, 1>&    f(const std::vector<Node>& globalNodes, Eigen::Matrix<double, 4, 1>& outf) const ;

    bool    CheckBoundaries(const std::vector<Node>& globalNodes) const;
};

std::istream& operator >> (std::istream& stream, FourNodeQuadrilateralElement& element);

#endif // ELEMENT_H_INCLUDED
