/*
 * Macro for calculating MCMC convergence factors for the
 * posteriors (without knowledge of the true target)
 */

#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "CommandHandler.h"
#include "ConfigManager.h"
#include "MCMCConvergence.h"

#include <iostream>


int main(int argc, char* argv[]) {

    std::cout << "Analyzing posteriors..." << std::endl;

    // Process command line and read config file
    CommandHandler handler(argc, argv);
    ConfigManager configs(handler.fname_config.Data());

    // Form the posterior and target covariance file names based on the config file
    TString mcmc_fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                              configs.GetChainDir().Data(),\
                              configs.GetMCMCFileBase().Data(),\
                              configs.GetNSteps(),\
                              configs.GetNPars(),\
                              configs.GetBranchNumber(),\
                              configs.GetRunNumber());

    // Open the files
    std::cout << "Reading MCMC chain from " << mcmc_fname << std::endl;
    mcmcConvergence* analyzer = new mcmcConvergence();
    analyzer->SetDim(configs.GetNPars());
    analyzer->LoadChain(mcmc_fname);

    // Prepare output
    TString out_fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                             configs.GetOutDir().Data(),\
                             configs.GetDiagnosticFileBase().Data(),\
                             configs.GetNSteps(),\
                             configs.GetNPars(),\
                             configs.GetBranchNumber(),\
                             configs.GetRunNumber());
    analyzer->PrepOutput(out_fname);

    ////////////////////// TODO: NEED TO INTERFACE NLAGS FOR AUTOCORRELATIONS
    int nlags = 500;

    // Get autocorrelations and parameter traces
    std::cout << "Calculating autocorrelations..." << std::endl;
    analyzer->CalcAutocorrelations(nlags);
    std::cout << "Plotting parameter traces..." << std::endl;
    analyzer->TraceParams();

    std::cout << "\nDone! Outputs written to " << out_fname.Data() << std::endl;
    return 0;
}

