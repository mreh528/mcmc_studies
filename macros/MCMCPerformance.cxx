/*
 * Macro that runs studies on the MCMC performance
 */

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "mcmc.h"
#include "Covariance.h"
#include "CommandHandler.h"

#include <iostream>
#include <string>
#include <vector>
using namespace std;


int main(int argc, char* argv[]) {

    CommandHandler handler(argc, argv);

    std::cout << "Hello World!" << std::endl;

    return 0;
}
