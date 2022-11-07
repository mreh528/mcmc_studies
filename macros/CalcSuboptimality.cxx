/*
 * Macro for calculating the "suboptimality factor" of an MCMC
 * following the direction of Roberts and Rosenthal (2006)
 */

#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
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

    // Want to be able to see how the suboptimality factor changes with each iteration
    if (handler.GetRunNumber() < 0 || handler.GetBranchNumber() < 0) {
        std::cout << "ERROR: Run number or branch number not set at cmd line" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load the posterior and target covariance files based on the config file
    TString post_cov_fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                  configs.GetProposalCovDir().c_str(),\
                                  configs.GetProposalCovFileBase().c_str(),\
                                  configs.GetNSteps(),\
                                  configs.GetNPars(),\
                                  handler.GetBranchNumber(),\
                                  handler.GetRunNumber()-1);
    TString target_cov_fname = Form("%s%s_npars%d_branch%d_target.root",\
                                    configs.GetTargetCovDir().c_str(),\
                                    configs.GetTargetCovFileBase().c_str(),\
                                    configs.GetNPars(),\
                                    handler.GetBranchNumber());

    // Create covariance objects for the posteriors and target
    std::cout << "Reading covariance matrices from "
              << post_cov_fname << " and " << target_cov_fname << std::endl;
    covariance* post_covmat = new covariance(configs.GetNPars());
    covariance* target_covmat = new covariance(configs.GetNPars());
    post_covmat->LoadCovPrev(post_cov_fname);
    target_covmat->LoadCovPrev(target_cov_fname);

    // Covmat from file gets stored in "Prev" elements in Covariance class
    TMatrixDSym* post_cov = post_covmat->GetPrevCovMat();
    TMatrixDSym* target_cov = target_covmat->GetPrevCovMat();
    TVectorD* post_means = post_covmat->GetPrevMeanVec();
    TVectorD* target_means = target_covmat->GetPrevMeanVec();

    // Get eigenvalues and eigenvectors
    TMatrixDSymEigen* post_eigen = new TMatrixDSymEigen(*post_cov);
    TMatrixD* eigenvectors = new TMatrixD(post_eigen->GetEigenVectors());
    TVectorD* eigenvalues  = new TVectorD(post_eigen->GetEigenValues());

    // Write covariance matrix to file
    TFile* output_cov_file = new TFile(Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                            configs.GetProposalCovDir().c_str(),\
                                            configs.GetProposalCovFileBase().c_str(),\
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

    std::cout << "Done!" << std::endl;

    return 0;
}

