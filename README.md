# mcmc_studies
Project studying Markov Chain Monte Carlo (MCMC) convergence performance at various dimensionalities.
To build, must have ROOT 5 installed. Built natively using gcc v6.1.0

Recommend making a ./build/ folder and running "cmake ../" from there, with executables being sent to a ./bin/ folder.

utils/ is the main library I am building. Here, CommandHandler handles command line arguments, ConfigManager parses config files, Covariance contains functions for calculating and updating covariance matrices, and mcmc maintains the MCMC.

macros/ holds the source code for the executables. CalcPosteriorCovmat calculates a covariance matrix for MCMC posteriors, GenerateTargetCovmat is used to generate random target covariance matrices to fit against, MCMCPerformance runs a MCMC, and PlotTraces plots the MCMC posteriors as a function of the step number.

Configuration files are stored in inputs/configs/, and batch job scripts are in job_scripts/

To run executables in bin/, most executables require you to specify a config file, run number, and branch number. For example,
./bin/mcmcPerformance -c inputs/configs/config.cfg -r 0 -b 0
runs a fresh MCMC on branch 0, run 0, using the parameters specified in the config file.
