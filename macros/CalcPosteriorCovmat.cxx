/*
 * Macro for updating the proposal covariance matrix for a MCMC
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
#include "ConfigManager.h"

#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

    // Process command line and read config file
    CommandHandler handler(argc, argv);
    ConfigManager configs(handler.fname_config.Data());

    if (handler.GetRunNumber() < 0 || handler.GetBranchNumber() < 0) {
        std::cout << "ERROR: Run number or branch number not set at cmd line" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load our inputs based on the config file
    TString input_mcmc_fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                    configs.GetChainDir().Data(),\
                                    configs.GetMCMCFileBase().Data(),\
                                    configs.GetNSteps(),\
                                    configs.GetNPars(),\
                                    handler.GetBranchNumber(),\
                                    handler.GetRunNumber());
    TString prev_cov_fname = "";
    if (handler.GetRunNumber() > 0) {
        prev_cov_fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                              configs.GetProposalCovDir().Data(),\
                              configs.GetProposalCovFileBase().Data(),\
                              configs.GetNSteps(),\
                              configs.GetNPars(),\
                              handler.GetBranchNumber(),\
                              handler.GetRunNumber()-1);
    }

    // Create the covariance object and load the inputs into it
    covariance* covmat = new covariance(configs.GetNPars());
    covmat->LoadChain(input_mcmc_fname);
    if (handler.GetRunNumber() > 0) { covmat->LoadCovPrev(prev_cov_fname); }

    // Prepare our writing objects
    std::cout << "Calculating posterior covariance matrix..." << std::endl;
    TMatrixDSym* post_cov = NULL;
    TVectorD* post_means = NULL;

    // If it's a new chain, calc covariance from scratch
    if (handler.GetRunNumber() == 0) {
        covmat->CalcCovariance(); // Means calc'd implicitly
        post_means = covmat->GetMeanVec();
        post_cov = covmat->GetCovMat();
    }
    // If not a new chain, update covariance from previous run
    else {
        post_means = covmat->UpdateMeans(handler.GetRunNumber());
        post_cov = covmat->UpdateCovariance(handler.GetRunNumber());
    }

    // Write covariance matrix to file
    TFile* output_cov_file = new TFile(Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                            configs.GetProposalCovDir().Data(),\
                                            configs.GetProposalCovFileBase().Data(),\
                                            configs.GetNSteps(),\
                                            configs.GetNPars(),\
                                            handler.GetBranchNumber(),\
                                            handler.GetRunNumber()),\
                                            "RECREATE");
    std::cout << "Covariance matrix and mean vector generated. Saving output to "
              << output_cov_file->GetName() << "..." << std::endl;
    post_cov->Write("cov_mat");
    post_means->Write("mean_vec");
    output_cov_file->Close();

    std::cout << "Done!\n" << std::endl;

    return 0;
}

