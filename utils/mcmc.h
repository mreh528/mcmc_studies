#ifndef _MCMC_H_
#define _MCMC_H_

#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TString.h"

// Class to run and manage MCMC
class mcmc {
public:
    mcmc();
    mcmc(const char* fname);
    ~mcmc();
    
    void RunMCMC();

private:    
    // MCMC running functions
    void LoadPrevChain(const char* fname);
    void InitCovMat();
    void ProposeStep();
    bool CheckProposedStep();
    void AcceptStep();
    
    // File I/O
    TFile* input_file;
    TFile* output_file;
    TChain* mcmc_chain;

    // Covariance matrices
    TMatrixDSym* input_covariance; // raw covmat from previous N steps
    TMatrixDSym* proposal_covariance; // modified version of input used for step proposal
    TMatrixDSym* proposal_cholesky_decomp;

    // Parameter vectors
    TVectorD* current_pars;
    TVectorD* proposed_pars;

    // Proposal function variables
    Float_t lnl_current;
    Float_t lnl_proposed;

    // General MCMC variables
    Float_t global_step_size; // global step size; should be 2.38/sqrt(npars)
    Float_t epsilon; // size of the regularization term applied to proposal covariance
    Int_t npars; // number of parameters to run with (number of dimensions)
    Int_t nsteps; // number of steps to run in this chain
    Int_t nstart; // starting step number (# of accumulated steps + 1)
};

#endif

