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
    mcmc(ConfigManager* configs, int _nrun=-1, int _nbranch=-1);
    ~mcmc();

    void RunMCMC();

private:
    // Constructor helper functions
    void SetDefault();
    void ReadConfigs(ConfigManager* configs, int _nrun, int _nbranch);

    // MCMC loading functions
    void LoadPrevChain();
    void InitCovMats();
    void PrepareOutput();
    void SaveChain();

    // MCMC Running Functions
    void ProposeStep();
    void CalcPDF();
    bool CheckIfAccepted();
    void AcceptStep();
    void RejectStep();

    // General File I/O
    TString chain_dir;
    TString input_chain_fname;
    TString output_chain_fname;
    TFile* input_chain_file;
    TFile* output_chain_file;
    TTree* mcmc_chain;

    // Covariance matrices for the target distribution
    TString target_cov_dir;
    TString target_cov_fname;
    TFile* target_cov_file;
    TMatrixDSym* target_covmat;
    TMatrixDSym* target_covmat_inverted;
    TVectorD* target_means;

    // Covariance matrices for the proposal function
    TString proposal_cov_dir;
    TString proposal_cov_fname;
    TFile* proposal_cov_file;
    TMatrixDSym* proposal_covmat; // modified version of empirical covmat used for step proposal
    TMatrixDSym* identity_matrix; // npars x npars identity matrix
    TMatrixD* proposal_cholesky; // U^T from the cholesky decomposition
    void GetCholDecomp();

    // MCMC Parameter vectors
    TVectorD* current_pars;
    TVectorD* proposed_pars;

    // Proposal function variables
    Float_t lnl_current;
    Float_t lnl_proposed;

    // General MCMC variables
    Float_t global_step_size; // global step size; should be 2.38^2/npars
    Float_t epsilon; // size of the regularization term applied to proposal covariance
    Int_t npars; // number of parameters to run with (number of dimensions)
    Int_t nsteps; // number of steps to run in this chain
    Int_t nstep_current; // Current step number
    Int_t nstep_start; // Starting step number
    Int_t nrun_current; // Current run number
    Int_t nrun_previous; // Previous run number (if reading from file)
    Int_t branch; // branch index (for parallel chains)
    Int_t naccepted; // # of accepted steps in the MCMC
    bool new_chain; // Are we starting a fresh chain?
    bool greedy; // Greedy only fills output on accepted step
    bool adaptive; // Running AM or non-adaptive?

    // Others
    TRandom3* rng;
    TStopwatch clock;
};

#endif

