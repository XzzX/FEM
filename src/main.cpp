#include    <cassert>
#include    <fstream>
#include    <iostream>
#include    <vector>

#include    "node.h"
#include    "element.h"

int main(){
    std::vector<Node>                           globalNodes; ///< all nodes in global numbering
    std::vector<FourNodeQuadrilateralElement>   elements;    ///< list of all elements

    std::cout << "reading nodes ... ";
    Node node;
    std::ifstream    iFile("nodes.txt");
    if (!iFile) std::cerr << "unknown node file" << std::endl;
    while   (!iFile.eof()){
        iFile >> node;
        globalNodes.push_back(node);
    }
    iFile.close();
    std::cout << globalNodes.size() << " read" << std::endl;

    std::cout << "reading elements ... ";
    FourNodeQuadrilateralElement element;
    iFile.open("elements.txt");
    if (!iFile) std::cerr << "unknown element file" << std::endl;
    while   (!iFile.eof()){
        iFile >> element;
        elements.push_back(element);
    }
    iFile.close();
    std::cout << elements.size() << " read" << std::endl;

    std::cout << "sorting nodes ... ";
    //E stands for essential nodes (dirichlet boundary)
    //reordering node list so that essential nodes are first
    //to not corrupt global node indices in elements the new ordering is just saved inside the nodes
    int nextE = 0;
    int nextF = globalNodes.size() - 1;
    for (Node& nd : globalNodes){
        if (nd.mBC.mType == BCType::D){
            nd.mPosition = nextE;
            nextE++;
        } else {
            nd.mPosition = nextF;
            nextF--;
        }
    }
    assert(nextE-nextF == 1);
    std::cout << "done" << std::endl;

    std::cout << "partition at " << nextE << std::endl;
    std::cout << "partitioning and assembling ... ";
    //dE consists of prescribed boundaries
    Eigen::MatrixXd d(globalNodes.size(), 1);
    for (const Node& nd : globalNodes){
        if (nd.mBC.mType == BCType::D){
            //if first sorting did not fail this should never be false
            assert(nd.mPosition < nextE);
            //for scalar problems only mX is considered boundary value
            d(nd.mPosition, 0) = nd.mBC.mX;
        }
    }

    FourNodeQuadrilateralElement::typeK outK;
    Eigen::Matrix<double, 4, 1>         outf;
    Eigen::MatrixXd KFE(globalNodes.size() - nextE, nextE);
    Eigen::MatrixXd KF(globalNodes.size() - nextE, globalNodes.size() - nextE);
    Eigen::MatrixXd fF(globalNodes.size() - nextE, 1);
    KF.setZero();
    fF.setZero();
    for (const FourNodeQuadrilateralElement& ele : elements){
        ele.K(globalNodes, outK);
        ele.f(globalNodes, outf);
        for (unsigned int i = 0; i < ele.mNodes.size(); i++){
            if (globalNodes.at(ele.mNodes.at(i)).mPosition < nextE) continue;
            for (unsigned int j = 0; j < ele.mNodes.size(); j++){
                if (globalNodes.at(ele.mNodes.at(j)).mPosition < nextE){
                    KFE(globalNodes.at(ele.mNodes.at(i)).mPosition - nextE, globalNodes.at(ele.mNodes.at(j)).mPosition) += outK(i,j);
                } else {
                    KF(globalNodes.at(ele.mNodes.at(i)).mPosition - nextE, globalNodes.at(ele.mNodes.at(j)).mPosition - nextE) += outK(i, j);
                }
            }

            fF(globalNodes.at(ele.mNodes.at(i)).mPosition - nextE, 0) = outf(i,0);
        }
    }

    std::cout << "done" << std::endl;

    std::cout << "solving!!!" << std::endl;
    //perhaps looking for a more suitable solver
    d.bottomRows(globalNodes.size() - nextE) = KF.colPivHouseholderQr().solve(fF - KFE * d.topRows(nextE));
    std::cout << d << std::endl;

    return 0;
}
