#ifndef _MCMCCONVERGENCE_H_
#define _MCMCCONVERGENCE_H_

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"

#include <iostream>

// Class to calculated chain convergence metrics, like autocorrelations
class mcmcConvergence {
public:
    mcmcConvergence();
    mcmcConvergence(int npars);
    ~mcmcConvergence();

    void SetDim(int npars);
    void LoadChain(TString mcmc_fname);
    void PrepOutput(TString fname_out);

    void CalcAutocorrelations(int nlag);
    void TraceParams();
private:
    void SetDefault();

    void SetupBranches();
    void FillArrays();

    int ndims;
    int nsteps;

    // Holds MCMC root file info
    TString* branch_names;
    double* branch_values;
    TFile* mcmc_file;
    TTree* mcmc_chain;

    // Holds MCMC steps in memory instead of disk
    double** par_vals;
    Float_t* lnls;
    Float_t lnl_current;

    // Output plots
    TH1D** param_trace_plots;
    TH1D** autocorrelation_plots;

    // Output file info
    TString output_name;
    TFile* output_file;
    TTree* output_tree;
};

#endif

