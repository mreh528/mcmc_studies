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

    std::string GetOutDir() { return main_output_directory; }
    std::string GetChainDir() { return main_output_directory + chain_directory; }
    std::string GetProposalCovDir() { return main_output_directory + proposal_cov_directory; }
    std::string GetTargetCovDir() { return main_output_directory + target_cov_directory; }
    std::string GetPlotDir() { return main_output_directory + plot_directory; }

    std::string GetMCMCFileBase() { return mcmc_file_base; }
    std::string GetTargetCovFileBase() { return covmat_file_base; }
    std::string GetProposalCovFileBase() { return covmat_file_base; }

    int GetNSteps() { return nsteps; }
    int GetNPars() { return npars; }
    double GetEpsilon() { return epsilon; }
private:
    std::string main_output_directory;
    std::string chain_directory;
    std::string proposal_cov_directory;
    std::string target_cov_directory;
    std::string plot_directory;

    std::string mcmc_file_base;
    std::string covmat_file_base;

    int nsteps;
    int npars;
    double epsilon;
};

#endif
