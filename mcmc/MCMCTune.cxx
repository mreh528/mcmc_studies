/*
 * Module for finding a good MCMC starting point
 */

#include "MCMCTune.h"


// Default Constructor
MCMCTune::MCMCTune() {
    SetDefault();
}


// Constructor from a config file
MCMCTune::MCMCTune(ConfigManager* configs) {
    SetDefault();
    ReadConfigs(configs);
    InitCovMats();
}


// Destructor
MCMCTune::~MCMCTune() {
    if (target_cov_file) { target_cov_file->Close(); }
    if (target_covmat) { delete target_covmat; }
    if (target_covmat_inverted) { delete target_covmat_inverted; }
    if (target_means) { delete target_means; }
    if (parameters) { delete parameters; }
    if (step_sizes) { delete step_sizes; }
}


// Sets default values of variables
void MCMCTune::SetDefault() {
    target_cov_dir = "";
    target_cov_fname = "";
    target_cov_file = NULL;

    target_covmat = NULL;
    target_covmat_inverted = NULL;
    target_means = NULL;

    parameters = NULL;
    step_sizes = NULL;

    npars = -1;
    branch = -1;

    lnl_current = -1.;
    lnl_previous = -1.;

    return;
}


// Reads config file and sets up variables accordingly
void MCMCTune::ReadConfigs(ConfigManager* configs) {
    std::cout << "Reading configs..." << std::endl;

    // Set up parameter vectors depending on the number of parameters
    if (configs->GetNPars() > 0) {
        npars = configs->GetNPars();
        parameters = new TVectorD(npars);
        step_sizes = new TVectorD(npars);
        // Set defaults to 1
        for (int ipar = 0; ipar < npars; ++ipar) {
            (*parameters)(ipar) = 1.;
            (*step_sizes)(ipar) = STARTING_STEP_SIZE;
        }
    } else {
        std::cout << "ERROR: number of parameters not specified." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Prepare target distribution file name for loading later
    if (configs->GetTargetCovDir().Length() > 0 &&
        configs->GetTargetCovFileBase().Length() > 0) {
        target_cov_dir = configs->GetTargetCovDir();
        target_cov_fname.Form("%s%s_npars%d_branch%d_target.root",
                              target_cov_dir.Data(),
                              configs->GetTargetCovFileBase().Data(),
                              npars, branch);
    } else {
        std::cout << "ERROR: No target distribution specified." << std::endl;
        exit(EXIT_FAILURE);
    }

    return;
}


// Initialize target distribution
void MCMCTune::InitCovMats() {
    std::cout << "Initializing target covariance matrix..." << std::endl;

    // Load up the target distribution covariance matrix
    if (target_cov_fname.Length() > 0) {
        target_cov_file = new TFile(target_cov_fname.Data(), "READ");
        if (!target_cov_file->IsOpen()) {
            std::cout << "ERROR: Invalid covmat file " << target_cov_fname << std::endl;
            exit(EXIT_FAILURE);
        }
        // Pull the matrices out of the file
        target_covmat = (TMatrixD*)target_cov_file->Get("cov_mat");
        target_means = (TVectorD*)target_cov_file->Get("mean_vec");
        if (target_covmat == NULL || target_means == NULL) {
            std::cout "ERROR: Target distribution not set from file." << std::endl;
            exit(EXIT_FAILURE);
        }
        // Get inverted target matrix for PDF calculation
        target_covmat_inverted = (TMatrixD*)target_covmat->Clone();
        target_covmat_inverted->Invert();
        std::cout << "Got Target distribution from "
                  << target_cov_file->GetName() << std::endl;
    } else {
        std::cout << "ERROR: No target file specified." << std::endl;
        exit(EXIT_FAILURE);
    }

    return;
}


// Golden Ratio minimization search in one dimension
void MCMCTune::MinimizeParam(int ipar) {
    Double_t step_size = STARTING_STEP_SIZE;
    // Init parameter guesses
    Double_t par_mid = (*parameters)(ipar);
    Double_t par_left = (*parameters)(ipar) - step_size;
    Double_t par_right = (*parameters)(ipar) + step_size;

    // Init lnls for all starting guesses
    Double_t lnl_mid = CalcPDF(parameters, ipar, par_mid);
    Double_t lnl_left = CalcPDF(parameters, ipar, par_left);
    Double_t lnl_right = CalcPDF(parameters, ipar, par_right);

    // Widen search until we bracket the minimum
    while ((lnl_left < lnl_mid) || (lnl_right < lnl_mid)) {
        step_size *= 2.;
        par_left = par_mid - step_size;
        lnl_left = CalcPDF(parameters, ipar, par_left);
        par_right = par_mid + step_size;
        lnl_right = CalcPDF(parameters, ipar, par_right);
    }

    // Initialize pars needed for golden ratio search
    Double_t par_midleft = par_left + (par_right - par_left)/GOLDEN_RATIO;
    Double_t par_midright = par_right - (par_right - par_left)/GOLDEN_RATIO;
    Double_t lnl_midleft = CalcPDF(parameters, ipar, par_midleft);
    Double_t lnl_midright = CalcPDF(parameters, ipar, par_midright);

    // Golden Ratio search until there's no significant change to be found
    while ((TMath::Abs(lnl_left-lnl_right) > TUNING_LNL_THRESHOLD
            || TMath::Abs(lnl_left-lnl_midleft) > TUNING_LNL_THRESHOLD
            || TMath::Abs(lnl_midleft-lnl_midright) > TUNING_LNL_THRESHOLD
            || TMath::Abs(lnl_midright-lnl_right) > TUNING_LNL_THRESHOLD)
           && TMath::Abs(par_midright-par_midleft) > TUNING_PAR_THESHOLD) {

        // If ML > MR, min lies between ML and R
        if (lnl_midleft > lnl_midright) {
            par_left = par_midleft;
            lnl_left = lnl_midleft;
            par_midleft = par_midright;
            lnl_midleft = lnl_midright;
            par_midright = par_left + (par_right - par_left)/GOLDEN_RATIO;
            lnl_midright = CalcPDF(parameters, ipar, par_midright);
        }
        // Otherwise, min lies between L and MR
        else {
            par_right = par_midright;
            lnl_right = lnl_midright;
            par_midright = par_midleft;
            lnl_midright = lnl_midleft;
            par_midleft = par_right - (par_right - par_left)/GOLDEN_RATIO;
            lnl_midleft = CalcPDF(parameters, ipar, par_midleft);
        }
    }

    (*parameters)(ipar) = (par_midleft + par_midright)/2.;
    return;
}


// Finds the 1D width of the distribution at LNL_MIN + TUNING_STEP_THRESHOLD
void MCMCTune::TuneStep(int ipar) {
    Double_t step_size = STARTING_STEP_SIZE;
    // Init parameter guesses
    Double_t par_mid = (*parameters)(ipar); // mid = min if par tuning has run already
    Double_t par_left = par_mid - step_size;
    Double_t par_right = par_mid + step_size;

    // Init lnls for all starting guesses
    Double_t lnl_mid = CalcPDF(parameters, ipar, par_mid);
    Double_t lnl_left = CalcPDF(parameters, ipar, par_left);
    Double_t lnl_right = CalcPDF(parameters, ipar, par_right);
    Double_t lnl_target = lnl_mid + TUNING_STEP_THRESHOLD;

    // Widen distribution until we bracket the lnl threshold
    while (lnl_left < lnl_target || lnl_right < lnl_target) {
        step_size *= 2.;
        par_left = par_mid - step_size;
        lnl_left = CalcPDF(parameters, ipar, par_left);
        par_right = par_mid + step_size;
        lnl_right = CalcPDF(parameters, ipar, par_right);
    }

    // Find left and right par values that match the lnl threshold
    par_left = BisectRoot(ipar, par_left, par_mid, lnl_target);
    par_right = BisectRoot(ipar, par_mid, par_right, lnl_target);

    // Set step size to the width of the distribution
    (*step_sizes)(ipar) = par_right - par_left;

    return;
}


// Use the bisection method to find a parameter value such that
// PDF(parval) ~= lnl_min + TUNING_STEP_THRESHOLD
Double_t MCMCTune::BisectRoot(int ipar, Double_t par_left, Double_t par_right, Double_t lnl_target) {
    Double_t par_mid = (par_left + par_right) / 2.;
    Double_t lnl_mid = CalcPDF(parameters, ipar, par_mid);
    Double_t lnl_left = CalcPDF(parameters, ipar, par_left);
    Double_t lnl_right = CalcPDF(parameters, ipar, par_right);

    // Loop until we hit the target precision
    while (TMath::Abs(lnl_mid - lnl_target) > TUNING_LNL_THRESHOLD) {
        // If mid is on the same side of target as left, left = mid
        if ((lnl_mid-lnl_target)*(lnl_left-lnl_target) > 0.) {
            par_left = par_mid;
            lnl_left = CalcPDF(parameters, ipar, par_left);
        }
        // Otherwise, if mid is on the same side as right, right = mid
        else {
            par_right = par_mid;
            lnl_right = CalcPDF(parameters, ipar, par_right);
        }
        par_mid = (par_left + par_right) / 2.;
        lnl_mid = CalcPDF(parameters, ipar, par_mid);
    }

    return par_mid;
}


// Finds an optimal starting position for the MCMC
void MCMCTune::TunePars() {
    lnl_current = CalcPDF(parameters);
    lnl_previous = lnl_current + 2.*TUNING_LNL_THRESHOLD; // make sure we at least take one pass

    // Iterate until no significant changes are seen after moving all parameters
    while (TMath::Abs(lnl_current-lnl_previous) > TUNING_LNL_THRESHOLD) {
        lnl_previous = lnl_current;
        for (int ipar = 0; ipar < npars; ++ipar) {
            MinimizeParam(ipar);
        }
        lnl_current = CalcPDF(parameters);
    }

    return;
}


// Finds an optimal set of starting step sizes for the MCMC
// Must be run after parameter tuning or else it won't work
void MCMCTune::TuneStepSizes() {
    for (int ipar = 0; ipar < npars; ++ipar) {
        TuneStep(ipar);
    }
    return;
}


// Finds optimal starting parameters and step sizes for the MCMC
void MCMCTune::TuneMCMC() {
    TunePars();
    TuneStepSizes();
    return;
}


// Calculate multivariate normal PDF in log form for the test parameters
Double_t MCMCTune::CalcPDF(TVectorD* par_vec, int ipar, Double_t parval) {

    TVectorD* diff_vec = (TVectorD*)par_vec->Clone();
    if (ipar > -1) { (*diff_vec)(ipar) = parval; }
    (*diff_vec) -= (*target_means);
    Double_t PDF = 0.5 * ( (*diff_vec) * ((*target_covmat_inverted) * (*diff_vec)) );

    return PDF;
}

