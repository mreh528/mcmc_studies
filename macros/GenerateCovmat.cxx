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
    
    std::cout << "Hello World!" << std::endl;

    return 0;
}

