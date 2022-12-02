/*
 * Macro to find a good MCMC starting point using some tuning algorithm
 */

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "MCMCTune.h"
#include "CommandHandler.h"
#include "ConfigManager.h"

#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

    // Handle configurations
    CommandHandler handler(argc, argv);
    ConfigManager configs(handler.fname_config.Data());

    std::cout << "Tuning MCMC using configs specified in " << handler.fname_config.Data() << std::endl;
    MCMCTune* tuner = new MCMCTune(&configs);
    tuner->TuneMCMC();

    // Load tuned values into objects for writing
    TVectorD* starting_pars = tuner->GetStartingPars();
    TVectorD* starting_steps = tuner->GetStartingSteps();
    TMatrixDSym* starting_cov = new TMatrixDSym(starting_steps->GetNoElements());
    starting_cov->Zero();
    // Starting steps become diagonal covariance matrix
    for (int i = 0; i < starting_steps->GetNoElements(); ++i) {
        (*starting_cov)(i,i) = (*starting_steps)(i);
    }

    TFile* fout = new TFile(Form("%s%s_npars%d_branch%d_tuned.root",\
                                 configs.GetProposalCovDir().Data(),\
                                 configs.GetProposalCovFileBase().Data(),\
                                 configs.GetNPars(),\
                                 configs.GetBranchNumber()),\
                                 "RECREATE");
    std::cout << "Writing outputs to " << fout->GetName() << std::endl;
    starting_pars->Write("mean_vec");
    starting_cov->Write("cov_mat");
    fout->Close();

    std::cout << "Finished!\n" << std::endl;
    return 0;
}

