/*
 * File to handle command line arguments and pass configurations to other calling functions.
 * Work in progress until this analysis becomes more fleshed out.
 */

#include "CommandHandler.h"

CommandHandler::CommandHandler(int argc, char* argv[]) {

    fname_in = "test.root";
    fname_out = "test.root";
    fname_config = "";
    verbose = false;

    for (int iarg = 0; iarg < argc; ++iarg) {
        // Print usage message
        if (std::string(argv[iarg]) == "-h" || std::string(argv[iarg]) == "--help" || argc < 2) {
            std::cout << "****************************************\n";
            std::cout << "Run Inputs:\n";
            std::cout << "    -h || --help              Print help info\n";
            std::cout << "    -i || --input [fname]     Input file name\n";
            std::cout << "    -o || --output [fname]    Output file name\n";
            std::cout << "    -c || --config [fname]    Configuration input file name\n";
            std::cout << "    -v || --verbose           If set, prints more help messages while running\n";
            std::cout << "****************************************" << std::endl;
            exit(0);
        }
        // Input file
        else if (std::string(argv[iarg]) == "-i" || std::string(argv[iarg]) == "--input") {
            fname_in = argv[++iarg];
        }
        // Output file
        else if (std::string(argv[iarg]) == "-o" || std::string(argv[iarg]) == "--output") {
            fname_out = argv[++iarg];
        }
        // Config file
        else if (std::string(argv[iarg]) == "-c" || std::string(argv[iarg]) == "--config") {
            fname_config = argv[++iarg];
        }
        // Verbosity
        else if (std::string(argv[iarg]) == "-v" || std::string(argv[iarg]) == "--verbose") {
            verbose = true;
        }
        // Default
        else {
            if (iarg > 0) { std::cout << "Unrecognized argument: " << argv[iarg] << "\n"; }
        }
    }

}

CommandHandler::~CommandHandler() {}

