/*
 * Macro for creating a randomized dxd covariance matrix
 * in d dimensions, as well as a corresponding random
 * vector of means in d dimensions
 */

#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "Covariance.h"
#include "CommandHandler.h"

#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

    CommandHandler handler(argc, argv);
    covariance* covmat = new covariance();
    std::cout << "Generating test vector..." << std::endl;
    TVectorD* test_vec = covmat->GetRandomMeanVec(10);
    for (int i = 0; i < 10; ++i) {
        std::cout << "Entry " << i << ": " << (*test_vec)(i) << std::endl;
    }
    
    std::cout << "\n\nHello World!" << std::endl;

    return 0;
}

