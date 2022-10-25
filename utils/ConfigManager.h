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

    std::string getOutputDirectory() { return output_directory; }
    std::string getMCMCFileBase() { return mcmc_file_base; }
    std::string getCovmatFileBase() { return covmat_file_base; }
    std::string getInputDirectory() { return input_directory; }

    bool getNewChainFlg() { return new_chain; }

    int getNSteps() { return nsteps; }
    int getNPars() { return npars; }
    double getEpsilon() { return epsilon; }
    int getRunNumber() { return run_number; }
    int getNBranches() { return nbranches; }
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
    int nbranches;
};

#endif
