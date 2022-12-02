#ifndef _MCMCTUNE_H_
#define _MCMCTUNE_H_

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "ConfigManager.h"

#include <iostream>

// Tuning threshold
static const Double_t TUNING_LNL_THRESHOLD = 1.;
static const Double_t TUNING_PAR_THRESHOLD = 0.000001;
static const Double_t TUNING_STEP_THRESHOLD = 16.;
static const Double_t STARTING_STEP_SIZE = 10.;
static const Double_t GOLDEN_RATIO = 0.5*(1.+TMath::Sqrt(5.));

// Class to help find a good MCMC starting point
class MCMCTune {
public:
    MCMCTune();
    MCMCTune(ConfigManager* configs);
    ~MCMCTune();

    void TuneMCMC();

    TVectorD* GetStartingPars() { return parameters; }
    TVectorD* GetStartingSteps() { return step_sizes; }
    Float_t GetStartingLnl() { return lnl_current; }
private:
    void SetDefault();
    void ReadConfigs(ConfigManager* configs);
    void InitCovMats();

    void MinimizeParam(int ipar);
    void TuneStep(int ipar);
    Double_t BisectRoot(int ipar, Double_t par_left, Double_t par_right, Double_t lnl_target);

    void TunePars();
    void TuneStepSizes();

    Double_t CalcPDF(TVectorD* par_vec, int ipar=-1, Double_t parval=0.);

    TString target_cov_dir;
    TString target_cov_fname;
    TFile* target_cov_file;
    TMatrixD* target_covmat;
    TMatrixD* target_covmat_inverted;
    TVectorD* target_means;

    TVectorD* parameters;
    TVectorD* step_sizes;
    Int_t npars;
    Int_t branch;
    Double_t lnl_current;
    Double_t lnl_previous;
};

#endif

