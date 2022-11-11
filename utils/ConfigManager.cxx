/*
 * Class for handling configuration file settings
 */

#include "ConfigManager.h"

// Default constructor
ConfigManager::ConfigManager() {
    SetDefault();
}


// Constructor using config file name
ConfigManager::ConfigManager(const char* fname) {
    // Set default values for variables in case they are not set in config
    SetDefault();

    // Now read set values from config
    std::cout << "\nReading configuration from " << fname << std::endl;
    if (!readConfig(fname)) {
        std::cout << "Something went wrong reading config file " << fname;
        std::cout << "... Exiting\n" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Successfully read configuration file!" << std::endl;
}


// Sets default values to empty or out of bounds
void ConfigManager::SetDefault() {
    main_output_directory = "";
    mcmc_file_base = "";
    covmat_file_base = "";
    chain_directory = "";
    proposal_cov_directory = "";
    target_cov_directory = "";
    nsteps = -1;
    npars = -1;
    epsilon = -1.;
}


// Destructor
ConfigManager::~ConfigManager() {}


// Main function that reads config file and holds data for later use
int ConfigManager::readConfig(const char* fname) {

    std::ifstream config_file(fname);
    if (!config_file.is_open()) {
        std::cout << "Could not open config file " << fname << std::endl;
        return 0;
    }

    std::string line_text;
    while (std::getline(config_file, line_text)) {
        std::istringstream line_stream(line_text);
        std::string key;
        if (std::getline(line_stream, key, ' ')) {
            if (key == "//") { continue; } // skip comments
            else if (key == "MAIN_OUTPUT_DIRECTORY") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                main_output_directory = key;
            }
            else if (key == "MCMC_FILE_BASE") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                mcmc_file_base = key;
            }
            else if (key == "COVMAT_FILE_BASE") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                covmat_file_base = key;
            }
            else if (key == "CHAIN_DIRECTORY") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                chain_directory = key;
            }
            else if (key == "PROPOSAL_COV_DIRECTORY") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                proposal_cov_directory = key;
            }
            else if (key == "PLOT_DIRECTORY") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                plot_directory = key;
            }
            else if (key == "TARGET_COV_DIRECTORY") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                target_cov_directory = key;
            }
            else if (key == "NSTEPS") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                nsteps = std::stoi(key);
            }
            else if (key == "NPARS") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                npars = std::stoi(key);
            }
            else if (key == "EPSILON") {
                std::getline(line_stream, key, ' '); // burn the '='
                std::getline(line_stream, key);
                epsilon = std::stod(key);
            }
        }
    }

    return 1;
}
