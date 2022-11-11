#ifndef _COVARIANCE_H_
#define _COVARIANCE_H_

#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include <iostream>

// Class to take care of covariance matrix handling
class covariance {
public:
    covariance();
    covariance(int npars, bool adapt=true);
    ~covariance();

    void SetDim(int npars);
    void LoadChain(TString mcmc_fname);
    void LoadCovPrev(TString cov_fname);
    void CalcCovariance();
    TMatrixDSym* UpdateCovariance();
    void CalcMeans();
    TVectorD* UpdateMeans();

    TMatrixDSym* GetRandomCovMat(int npars);
    TVectorD* GetRandomMeanVec(int npars);

    TVectorD* GetMeanVec() { return mean_vec; }
    TMatrixDSym* GetCovMat() { return cov_mat; }
    TVectorD* GetPrevMeanVec() { return mean_vec_prev; }
    TMatrixDSym* GetPrevCovMat() { return cov_mat_prev; }
    TVectorD* GetStepVec() { return step_vec; }
    int GetNDims() { return ndims; }
    Double_t GetNSteps() { return step_vec->Sum(); }

    void PrintMatrix(TMatrixDSym* mat);

private:
    void SetDefault();

    TString* branch_names;
    double* branch_values;
    TFile* mcmc_file;
    TFile* cov_file;
    TTree* mcmc_chain;
    bool chain_loaded;
    bool adaptive;
    int nsteps_current;
    int nsteps_previous;

    int ndims;
    TVectorD* mean_vec;
    TMatrixDSym* cov_mat;
    TVectorD* mean_vec_prev;
    TMatrixDSym* cov_mat_prev;
    TVectorD* step_vec; // step_vec(0) = previous nsteps, step_vec(1) = current nsteps
};

#endif

