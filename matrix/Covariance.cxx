/*
 * Module to handle mcmc covariance matrix
 */

#include "Covariance.h"


// Default constructor
covariance::covariance() {
    SetDefault();
    std::cout << "Warning: Covariance default constructor used. Most variables will be uninitialized." << std::endl;
}


// Constructor with only the parameter count specified
covariance::covariance(int npars) {
    // Make sure we don't get bad input
    if (npars < 1) {
        std::cout << "ERROR: invalid number of parameters specified: " << npars << std::endl;
        exit(EXIT_FAILURE);
    }

    // Initialize
    SetDefault();
    SetDim(npars);

    // Ensure all values are initialized at 0 instead of memory junk
    mean_vec->Zero();
    cov_mat->Zero();
}


// Destructor
covariance::~covariance() {
    if (branch_names) { delete branch_names; }
    if (branch_values) { delete branch_values; }
    if (mcmc_chain) { delete mcmc_chain; }
    if (mean_vec) { delete mean_vec; }
    if (cov_mat) { delete cov_mat; }
    if (mean_vec_prev) { delete mean_vec_prev; }
    if (cov_mat_prev) { delete cov_mat_prev; }
    if (cov_file) { cov_file->Close(); }
    if (mcmc_file) { mcmc_file->Close(); }
}


// Set default values of data
void covariance::SetDefault() {
    ndims = -1;
    nsteps_current = 0;
    nsteps_previous = 0;
    branch_names = NULL;
    branch_values = NULL;
    mcmc_file = NULL;
    cov_file = NULL;
    mcmc_chain = NULL;
    mean_vec = NULL;
    cov_mat = NULL;
    mean_vec_prev = NULL;
    cov_mat_prev = NULL;
    step_vec = new TVectorD(2);
    step_vec->Zero();
    chain_loaded = false;
    return;
}


// Setup covariance parameters
void covariance::SetDim(int npars) {
    // Watch out for bad input
    if (npars < 1) {
        std::cout << "ERROR: invalid number of dimensions specified: " << npars << std::endl;
        exit(EXIT_FAILURE);
    }

    // Setup everything that depends on the number of dimensions
    ndims = npars;
    branch_names = new TString[ndims];
    branch_values = new double[ndims];
    mean_vec  = new TVectorD(ndims);
    cov_mat  = new TMatrixDSym(ndims);

    // Wait to see if we have an input covariance file for these
    mean_vec_prev = NULL;
    cov_mat_prev = NULL;

    return;
}


// Load MCMC chain from file
void covariance::LoadChain(TString mcmc_fname) {
    // Watch out for bad input
    if (ndims < 1) {
        std::cout << "ERROR: invalid number of dimensions specified: " << ndims << std::endl;
        std::cout << "Did you use the correct constructor?" << std::endl;
        exit(EXIT_FAILURE);
    }

    chain_loaded = false;

    // Load mcmc posteriors chain from file
    mcmc_file = new TFile(mcmc_fname.Data(), "READ");
    if (!mcmc_file->IsOpen()) {
        std::cout << "ERROR: Invalid input file name " << mcmc_file->GetName() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Prepare branches for loading
    std::cout << "Loading MCMC chain from " << mcmc_file->GetName() << std::endl;
    mcmc_chain = (TTree*)mcmc_file->Get("posteriors");
    TObjArray* branch_list = (TObjArray*)mcmc_chain->GetListOfBranches();
    int nbranches = branch_list->GetEntries();

    // Set branch addresses so we can load them easily
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        TBranch* br = (TBranch*)branch_list->At(ibr);
        TString bname = br->GetName();
        for (int idim = 0; idim < ndims; ++idim) {
            if (bname.CompareTo(Form("mcmc_par_%d",idim))==0) {
                branch_names[idim] = br->GetName();
                mcmc_chain->SetBranchAddress(branch_names[idim].Data(), &branch_values[idim]);
                break;
            }
        }
    }

    nsteps_current = mcmc_chain->GetEntries();
    (*step_vec)(1) = (Double_t)nsteps_current;
    chain_loaded = true;

    return;
}


// Load a previous covariance matrix file
void covariance::LoadCovPrev(TString cov_fname) {

    cov_file = new TFile(cov_fname.Data(), "READ");
    if (!cov_file->IsOpen()) {
        std::cout << "ERROR: could not open input covariance file" << std::endl;
        exit(EXIT_FAILURE);
    }

    mean_vec_prev = (TVectorD*)cov_file->Get("mean_vec");
    cov_mat_prev = (TMatrixDSym*)cov_file->Get("cov_mat");
    TVectorD* step_vec_prev = (TVectorD*)cov_file->Get("step_vec");
    (*step_vec)(0) = (*step_vec_prev)(0);
}


// Calculates the covariance matrix of some input mcmc chain
void covariance::CalcCovariance() {
    // Make sure we're not running uninitialized
    if (ndims < 1) {
        std::cout << "ERROR: invalid number of dimensions specified: " << ndims << std::endl;
        std::cout << "Did you use the correct constructor?" << std::endl;
        exit(EXIT_FAILURE);
    } else if (!chain_loaded) {
        std::cout << "ERROR: no chain loaded to read from!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // First need means in order to calculate covariance
    CalcMeans();

    // Make sure everything is reset
    cov_mat->Zero();

    // Loop through all entries in the chain to calculate the covariance matrix
    std::cout << "Building new covariance matrix with " << nsteps_current
              << " mcmc steps." << std::endl;
    for (int ientry = 0; ientry < nsteps_current; ++ientry) {
        mcmc_chain->GetEntry(ientry);
        for (int ibr = 0; ibr < ndims; ++ibr) {
            for (int jbr = 0; jbr < ndims; ++jbr) {
                (*cov_mat)(ibr,jbr) += (branch_values[ibr] - (*mean_vec)(ibr)) \
                                     * (branch_values[jbr] - (*mean_vec)(jbr)) \
                                     / (*step_vec)(1);
            }
        }
    }

    return;
}


// Updates the value of an input covariance matrix using some loaded set of data.
// Must have previously loaded the data using LoadChain()
TMatrixDSym* covariance::UpdateCovariance() {
    if (!chain_loaded) {
        std::cout << "ERROR: Chain of new data not loaded. "
                  << "Cannot update input covmat until LoadChain() is called." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Step count has to be stored in a float-point vector unfortunately
    nsteps_previous = (int)((*step_vec)(0));
    std::cout << "Updating previous covariance matrix (built with " << nsteps_previous
              << " steps) with " << nsteps_current << " new steps." << std::endl;

    CalcCovariance(); // calc covariance for new data loaded in chain; stored in cov_mat
    TMatrixDSym* new_cov = (TMatrixDSym*)cov_mat->Clone();
    TMatrixDSym* tmp_cov = (TMatrixDSym*)cov_mat_prev->Clone();

    // Updated covariance matrix is a weighted average of the previous two
    (*new_cov) *= (((*step_vec)(1)-1.) / (step_vec->Sum()-1.));
    (*tmp_cov) *= (((*step_vec)(0)-1.) / (step_vec->Sum()-1.));
    // Fill temp matrix holding outer product of the two mean vector differences
    TMatrixDSym* mean_mat = new TMatrixDSym(ndims);
    for (int i = 0; i < ndims; ++i) {
        for (int j = 0; j < ndims; ++j) {
            (*mean_mat)(i,j) = ((*mean_vec)(i)-(*mean_vec_prev)(i)) \
                             * ((*mean_vec)(j)-(*mean_vec_prev)(j)) \
                             * (((*step_vec)(1)*(*step_vec)(0)) \
                                / ((step_vec->Sum()-1.) * step_vec->Sum()));
        }
    }

    // Add together the three matrices to finish the weighted average
    (*new_cov) += (*tmp_cov);
    (*new_cov) += (*mean_mat);

    return new_cov;
}


// Calculates the mean vector of some input data set
void covariance::CalcMeans() {
    // Make sure we're not running uninitialized
    if (ndims < 1) {
        std::cout << "ERROR: invalid number of dimensions specified: " << ndims << std::endl;
        std::cout << "Did you use the correct constructor?" << std::endl;
        exit(EXIT_FAILURE);
    } else if (!chain_loaded) {
        std::cout << "ERROR: No chain loaded to read from!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Make sure everything is reset
    mean_vec->Zero();

    // Loop through all entries in the chain to calculate the mean of each par
    for (int ientry = 0; ientry < nsteps_current; ++ientry) {
        mcmc_chain->GetEntry(ientry);
        for (int ibr = 0; ibr < ndims; ++ibr) {
            (*mean_vec)(ibr) += branch_values[ibr]/(*step_vec)(1);
        }
    }

    return;
}


// Updates the values of an input mean vector using some loaded set of data.
// Must have loaded the data previously using LoadChain()
TVectorD* covariance::UpdateMeans() {
    if (!chain_loaded) {
        std::cout << "ERROR: Chain of new data not loaded. "
                  << "Cannot update input means until LoadChain() is called." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Step count has to be stored in a float-point vector unfortunately
    nsteps_previous = (int)((*step_vec)(0));

    CalcMeans(); // Calculate means for new data loaded in chain
    TVectorD* new_means = (TVectorD*)mean_vec->Clone();
    TVectorD* tmp_vec = (TVectorD*)mean_vec_prev->Clone();

    // New mean vector is the weighted average of the old + new
    *new_means *= (*step_vec)(1) / step_vec->Sum();
    *tmp_vec *= (*step_vec)(0) / step_vec->Sum();
    *new_means += *tmp_vec;

    return new_means;
}


// Generate a random covariance matrix of a specified dimension
TMatrixDSym* covariance::GetRandomCovMat(int npars) {
    TMatrixDSym* rand_cov = new TMatrixDSym(npars);
    TRandom3* rng = new TRandom3();
    rng->SetSeed(); // Defaults to a random seed
    for (int i = 0; i < npars; ++i) {
        for (int j = i; j < npars; ++j) {
            (*rand_cov)(i,j) = rng->Gaus(0., 5.);
        }
    }

    // Have to multiply by transpose to guarantee positive definite
    TMatrixDSym* tmp = (TMatrixDSym*)rand_cov->Clone();
    tmp->T();
    rand_cov->TMult(*tmp);
    return rand_cov;
}


// Generate a random mean vector of a specified dimension
TVectorD* covariance::GetRandomMeanVec(int npars) {
    TVectorD* rand_means = new TVectorD(npars);
    TRandom3* rng = new TRandom3();
    rng->SetSeed(); // Defaults to random seed
    for (int i = 0; i < npars; ++i) {
        (*rand_means)(i) = rng->Gaus(0., 50.);
    }

    return rand_means;
}


// Prints the values of the input matrix to the screen
void covariance::PrintMatrix(TMatrixDSym* mat) {

    for (int i = 0; i < mat->GetNrows(); ++i) {
        std::cout << "  ";
        for (int j = 0; j < mat->GetNcols(); ++j) {
            std::cout << (*mat)(i,j) << "  ";
        }
        std::cout << std::endl;
    }

    return;
}

