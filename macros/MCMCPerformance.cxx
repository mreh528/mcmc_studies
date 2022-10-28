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
    ConfigManager configs;
    if (handler.fname_config) { configs.readConfig(handler.fname_config.Data()); }

    std::cout << "Hello World!" << std::endl;
    mcmc* markov_chain = new mcmc(&configs);
    delete markov_chain;

    return 0;
}
