#ifndef _COVARIANCE_H_
#define _COVARIANCE_H_

#include "TMath.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TObjArray.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include <iostream>

// Class to take care of covariance matrix handling
class covariance {
public:
    covariance();
    covariance(int npars, char* fname);
    ~covariance();

    void SetDim(int npars);
    void LoadChain(char* fname);
    void CalcCovariance();
    TMatrixDSym* UpdateCovariance(TMatrixDSym* cov_mat_prev, TVectorD* mean_vec_prev, int nsteps_prev);
    void CalcMeans();
    TVectorD* UpdateMeans(TVectorD* mean_vec_prev, int nsteps_prev);

    TMatrixDSym* GetRandomCovMat(int npars);
    TVectorD* GetRandomMeanVec(int npars);

    TVectorD* GetMeanVec() { return mean_vec; }
    TMatrixDSym* GetCovMat() { return cov_mat; }
    int GetNDims() { return ndims; }

private:
    TString* branch_names;
    double* branch_values;
    TChain* mcmc_chain;
    bool chain_loaded;
    int nsteps;

    int ndims;
    TVectorD* mean_vec;
    TMatrixDSym* cov_mat;
};

#endif

