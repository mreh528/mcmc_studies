#ifndef _CONFIGMANAGER_H_
#define _CONFIGMANAGER_H_

#include "TString.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>

class ConfigManager {
public:
    ConfigManager();
    ConfigManager(const char* fname);
    ~ConfigManager();

    int readConfig(const char* fname);

    TString GetOutDir() { return main_output_directory; }
    TString GetChainDir() {
        TString return_string = main_output_directory;
        return_string += chain_directory;
        return return_string;
    }
    TString GetProposalCovDir() {
        TString return_string = main_output_directory;
        return_string += proposal_cov_directory;
        return return_string;
    }
    TString GetTargetCovDir() {
        TString return_string = main_output_directory;
        return_string += target_cov_directory;
        return return_string;
    }
    TString GetPlotDir() {
        TString return_string = main_output_directory;
        return_string += plot_directory;
        return return_string;
    }

    TString GetMCMCFileBase() { return mcmc_file_base; }
    TString GetTargetCovFileBase() { return covmat_file_base; }
    TString GetProposalCovFileBase() { return covmat_file_base; }

    int GetNSteps() { return nsteps; }
    int GetNPars() { return npars; }
    int GetBranchNumber() { return branch; }
    double GetEpsilon() { return epsilon; }
    bool GreedyAcceptance() { return greedy; }
    bool AdaptiveMetropolis() { return adaptive; }
private:
    void SetDefault();

    TString main_output_directory;
    TString chain_directory;
    TString proposal_cov_directory;
    TString target_cov_directory;
    TString plot_directory;

    TString mcmc_file_base;
    TString covmat_file_base;

    int nsteps;
    int npars;
    int branch;
    double epsilon;
    bool greedy;
    bool adaptive;
};

#endif
