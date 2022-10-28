#ifndef _MCMC_H_
#define _MCMC_H_

#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TString.h"

#include "Covariance.h"
#include "ConfigManager.h"

// Class to run and manage MCMC
class mcmc {
public:
    mcmc();
    mcmc(ConfigManager* configs);
    ~mcmc();

    void RunMCMC(double _lnl_current = -1.);

private:
    void SetDefault(); // constructor function to set everything out-of-bounds

    // MCMC loading functions
    void LoadPrevChain();
    void InitNewChain();
    void InitCovMats(TString covmat_fname_base, bool new_chain);
    void PrepareOutput();
    void SaveChain();

    // MCMC Running Functions
    void ProposeStep();
    Float_t CalcPDF();
    bool CheckIfAccepted();
    void AcceptStep();
    void RejectStep();

    // File I/O
    TString input_dir;
    TString output_dir;
    TString input_fname;
    TString output_fname;
    TFile* input_file;
    TFile* output_file;
    TTree* mcmc_chain;

    // Covariance matrices
    TString target_cov_fname;
    TFile* target_cov_file;
    TMatrixDSym* target_covmat;
    TMatrixDSym* target_covmat_inverted;
    TVectorD* target_means;
    TString proposal_cov_fname;
    TFile* proposal_cov_file;
    TMatrixDSym* proposal_covmat; // modified version of empirical covmat used for step proposal
    TMatrixD* proposal_cholesky; // U^T from the cholesky decomposition
    TMatrixDSym* identity_matrix; // npars x npars identity matrix
    TMatrixD GetCholDecomp();

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
    Int_t nstep_current; // Current step number
    Int_t nrun_current; // Current run number
    Int_t nrun_previous; // Previous run number (if reading from file)
    Int_t branch; // branch index (for parallel chains)
    Int_t naccepted; // # of accepted steps in the MCMC

    // Others
    TRandom3* rng;
    TStopwatch clock;
};

#endif

