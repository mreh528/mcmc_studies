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
#include "ConfigManager.h"

#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

    // Process command line and read config file
    CommandHandler handler(argc, argv);
    ConfigManager configs(handler.fname_config.Data());
    covariance* covmat = new covariance();

    if (handler.GetBranchNumber() < 0) {
        std::cout << "ERROR: Branch number not set at cmd line" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Generate a random covariance matrix and mean vector with dimension specified in configs
    std::cout << "Generating target distribution for "
              << configs.GetNPars() << " pars..." << std::endl;
    TMatrixDSym* random_cov = covmat->GetRandomCovMat(configs.GetNPars());
    TVectorD* random_means = covmat->GetRandomMeanVec(configs.GetNPars());

    // Write covariance matrix to file
    std::cout << "Covariance matrix and mean vector generated. Saving output..." << std::endl;
    TFile* output = new TFile(Form("%s%s_npars%d_branch%d_target.root",\
                                   configs.GetTargetCovDir().c_str(),\
                                   configs.GetTargetCovFileBase().c_str(),\
                                   configs.GetNPars(),\
                                   handler.GetBranchNumber()),\
                                   "RECREATE");
    random_cov->Write("cov_mat");
    random_means->Write("mean_vec");
    output->Close();

    // Clean up memory
    delete covmat;
    delete random_cov;
    delete random_means;

    std::cout << "Done! Output written to " << output->GetName() << std::endl;

    return 0;
}

