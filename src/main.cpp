#include    <cassert>
#include    <fstream>
#include    <iostream>
#include    <vector>

#include	"configuration.h"
#include    "node.h"
#include    "element.h"


///!!!!!!!!!!!!ATTENTION!!!!!!!!!!!!!!!!
///Pay attention to std::vector relocation after push_back especially for refinement loop -> kills this pointer!!!!!!!
///!!!!!!!!!!!!ATTENTION!!!!!!!!!!!!!!!!
int main(unsigned int argc, char** argv){
	gConfig.ReadCommandLineParameters(argc, argv);

    std::vector<Node>                           globalNodes; ///< all nodes in global numbering
	globalNodes.reserve(1000);
    std::vector<FourNodeQuadrilateralElement>   elements;    ///< list of all elements
	elements.reserve(1000);

    std::cout << "reading nodes ... ";
    Node node;
    std::ifstream    iFile(gConfig.mFilename + "_nodes.txt");
    if (!iFile) std::cerr << "unknown node file" << std::endl;
    while   (!iFile.eof()){
        iFile >> node;
        globalNodes.push_back(node);
    }
    iFile.close();
    std::cout << globalNodes.size() << " read" << std::endl;

    std::cout << "reading elements ... ";
    FourNodeQuadrilateralElement element;
    iFile.open(gConfig.mFilename + "_elements.txt");
    if (!iFile) std::cerr << "unknown element file" << std::endl;
    while   (!iFile.eof()){
        iFile >> element;
        elements.push_back(element);
    }
    iFile.close();
    std::cout << elements.size() << " read" << std::endl;

    std::cout << "refining mesh ... ";
	for (unsigned int n = 0; n < gConfig.mNumberOfRefinements; n++){
		unsigned int	maxElements = elements.size();
		for (unsigned int i = 0; i < maxElements; i++){
			elements.at(i).Refine(globalNodes, elements);
		}
	}
    std::cout << "done" << std::endl;

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
    Eigen::MatrixXd globalD(globalNodes.size(), 1);
    for (const Node& nd : globalNodes){
        if (nd.mBC.mType == BCType::D){
            //if first sorting did not fail this should never be false
            assert(nd.mPosition < nextE);
            //for scalar problems only mX is considered boundary value
            globalD(nd.mPosition, 0) = nd.mBC.mX;
        }
    }

    FourNodeQuadrilateralElement::typeK outK;
    Eigen::Matrix<double, 4, 1>         outf;
    Eigen::MatrixXd KFE(globalNodes.size() - nextE, nextE);
	KFE.setZero();
    Eigen::MatrixXd KF(globalNodes.size() - nextE, globalNodes.size() - nextE);
	KF.setZero();
    Eigen::MatrixXd fF(globalNodes.size() - nextE, 1);
    fF.setZero();
	for (unsigned int k = 0; k < elements.size(); k++) {
        elements.at(k).K(globalNodes, outK);
		elements.at(k).f(globalNodes, outf);

		for (unsigned int i = 0; i < elements.at(k).mNodes.size(); i++) {
			if (globalNodes.at(elements.at(k).mNodes.at(i)).mPosition < nextE) continue;
			for (unsigned int j = 0; j < elements.at(k).mNodes.size(); j++) {
				if (globalNodes.at(elements.at(k).mNodes.at(j)).mPosition < nextE) {
					KFE(globalNodes.at(elements.at(k).mNodes.at(i)).mPosition - nextE, globalNodes.at(elements.at(k).mNodes.at(j)).mPosition) += outK(i, j);
                } else {
					KF(globalNodes.at(elements.at(k).mNodes.at(i)).mPosition - nextE, globalNodes.at(elements.at(k).mNodes.at(j)).mPosition - nextE) += outK(i, j);
                }
            }

			fF(globalNodes.at(elements.at(k).mNodes.at(i)).mPosition - nextE, 0) = outf(i, 0);
        }
    }

    std::cout << "done" << std::endl;

    std::cout << "solving!!!" << std::endl;
    //perhaps looking for a more suitable solver
	globalD.bottomRows(globalNodes.size() - nextE) = KF.colPivHouseholderQr().solve(fF - KFE * globalD.topRows(nextE));

	std::cout << "outputting temperature map ... ";
	std::ofstream   ofile(gConfig.mFilename + "_temperaturemap.txt");
	for (const Node& nd : globalNodes) {
		ofile << nd.mX << "\t" << nd.mY << "\t" << globalD(nd.mPosition) << std::endl;
	}
	ofile.close();
	std::cout << "done" << std::endl;

	std::cout << "post processing ... ";
	ofile.open(gConfig.mFilename + "_flux.txt");
	for (const auto& ele : elements) {
		Eigen::Matrix<double, 4, 4> temp;
		ele.PostProcess(globalNodes, globalD, temp);
		for (unsigned int i = 0; i < 4; i++)
			ofile << temp(i, 0) << "\t" << temp(i, 1) << "\t" << temp(i, 2) << "\t" << temp(i, 3) << std::endl;
	}
	ofile.close();
	std::cout << "done" << std::endl;

    return 0;
}
