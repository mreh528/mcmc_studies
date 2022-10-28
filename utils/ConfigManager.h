#ifndef _CONFIGMANAGER_H_
#define _CONFIGMANAGER_H_

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

    std::string GetOutputDirectory() { return output_directory; }
    std::string GetMCMCFileBase() { return mcmc_file_base; }
    std::string GetCovmatFileBase() { return covmat_file_base; }
    std::string GetInputDirectory() { return input_directory; }

    bool GetNewChainFlg() { return new_chain; }

    int GetNSteps() { return nsteps; }
    int GetNPars() { return npars; }
    double GetEpsilon() { return epsilon; }
    int GetRunNumber() { return run_number; }
    int GetNBranch() { return nbranch; }
private:
    std::string output_directory;
    std::string mcmc_file_base;
    std::string covmat_file_base;
    std::string input_directory;

    bool new_chain;

    int nsteps;
    int npars;
    double epsilon;
    int run_number;
    int nbranch;
};

#endif
