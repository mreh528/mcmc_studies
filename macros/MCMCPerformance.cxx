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
#include "ConfigManager.h"

#include <iostream>
#include <string>
#include <vector>
using namespace std;


int main(int argc, char* argv[]) {

    CommandHandler handler(argc, argv);
    ConfigManager configs(handler.fname_config.Data());

    std::cout << "Creating MCMC using configs specified in " << handler.fname_config.Data() << std::endl;
    mcmc* markov_chain = new mcmc(&configs, handler.GetRunNumber(), handler.GetBranchNumber());

    std::cout << "Starting MCMC" << std::endl;
    markov_chain->RunMCMC();

    std::cout << "Finished!\n" << std::endl;
    return 0;
}
