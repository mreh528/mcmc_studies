#ifndef _COMMANDHANDLER_h_
#define _COMMANDHANDLER_h_

#include "TString.h"

#include <iostream>
#include <string>

// Class to handle command line arguments
class CommandHandler {
public:
    CommandHandler(int argc, char* argv[]);
    ~CommandHandler();
    TString fname_in;
    TString fname_out;
    TString fname_config;
    bool verbose;

    int GetRunNumber() { return nrun; }
private:
    // MCMC-specific variables
    int nrun;
};

#endif
