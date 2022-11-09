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
#include "Metrics.h"

#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

    std::cout << "Analyzing posteriors..." << std::endl;

    // Process command line and read config file
    CommandHandler handler(argc, argv);
    ConfigManager configs(handler.fname_config.Data());

    // Want to be able to see how the suboptimality factor changes with each iteration
    if (handler.GetRunNumber() < 0 || handler.GetBranchNumber() < 0) {
        std::cout << "ERROR: Run number or branch number not set at cmd line" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Form the posterior and target covariance file names based on the config file
    TString post_cov_fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                  configs.GetProposalCovDir().c_str(),\
                                  configs.GetProposalCovFileBase().c_str(),\
                                  configs.GetNSteps(),\
                                  configs.GetNPars(),\
                                  handler.GetBranchNumber(),\
                                  handler.GetRunNumber());
    TString target_cov_fname = Form("%s%s_npars%d_branch%d_target.root",\
                                    configs.GetTargetCovDir().c_str(),\
                                    configs.GetTargetCovFileBase().c_str(),\
                                    configs.GetNPars(),\
                                    handler.GetBranchNumber());

    // Open the files
    std::cout << "Reading posterior covariance from " << post_cov_fname
              << ", and target covariance from" << target_cov_fname << std::endl;
    TFile* post_cov_file = new TFile(post_cov_fname.Data(), "READ");
    TFile* target_cov_file = new TFile(target_cov_fname.Data(), "READ");
    if (!post_cov_file->IsOpen() || !target_cov_file->IsOpen()) {
        std::cout << "ERROR: posterior file or target file could not be opened" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load the covariance matrices and mean vectors
    TMatrixD* post_cov = (TMatrixD*)post_cov_file->Get("cov_mat");
    TVectorD* post_means = (TVectorD*)post_cov_file->Get("mean_vec");
    TMatrixD* target_cov = (TMatrixD*)target_cov_file->Get("cov_mat");
    TVectorD* target_means = (TVectorD*)target_cov_file->Get("mean_vec");

    // Pre-calculate all metric values to help clear final output of error messages
    double RRSF = RRSuboptimalityFactor(post_cov, target_cov);
    double AIRM = AffineInvariantRiemannianMetric(post_cov, target_cov);
    double LERM = LogEuclideanRiemannMetric(post_cov, target_cov);
    double WM   = WassersteinMetric(post_cov, post_means, target_cov, target_means);
    double WM2  = WassersteinNoMeans(post_cov, target_cov);

    // Calculate and print out the various convergence metrics
    std::cout << std::endl;
    std::cout << "Roberts-Rosenthal Suboptimality factor: " << RRSF << std::endl;
    std::cout << "Affine Invariant Riemannian Metric:     " << AIRM << std::endl;
    std::cout << "Log Euclidean Riemann Metric:           " << LERM << std::endl;
    std::cout << "Wasserstein Metric:                     " << WM   << std::endl;
    std::cout << "Wasserstein Metric without Means:       " << WM2  << std::endl;

    std::cout << "\nDone!" << std::endl;
    return 0;
}

