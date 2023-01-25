/*
 * Module to calculate convergence metrics on the mcmc
 */

#include "MCMCConvergence.h"


// Default constructor
mcmcConvergence::mcmcConvergence() {
    SetDefault();
}


// Default constructor
mcmcConvergence::mcmcConvergence(int npars) {
    SetDefault();
    SetDim(npars);
}


// Destructor
mcmcConvergence::~mcmcConvergence() {
    if (branch_names) { delete branch_names; }
    if (branch_values) { delete branch_values; }
    if (mcmc_chain) { delete mcmc_chain; }
    if (mcmc_file) { mcmc_file->Close(); }

    if (lnls) { delete lnls; }
    if (par_vals) {
        for (int i = 0; i < nsteps; ++i) {
            delete[] par_vals[i];
        }
        delete[] par_vals;
    }

    if (output_file) { output_file->Close(); }
    if (output_tree) { delete output_tree; }
}


// Set default values of data
void mcmcConvergence::SetDefault() {
    ndims = -1;
    nsteps = -1;

    branch_names = NULL;
    branch_values = NULL;
    mcmc_file = NULL;
    mcmc_chain = NULL;

    par_vals = NULL;
    lnls = NULL;
    lnl_current = -1.;

    param_trace_plots = NULL;
    autocorrelation_plots = NULL;

    output_file = NULL;
    output_tree = NULL;

    return;
}


// Setup chain parameters
void mcmcConvergence::SetDim(int npars) {
    // Watch for bad input
    if (npars < 1) {
        std::cout << "ERROR: invalid number of dimensions specified: " << npars << std::endl;
        exit(EXIT_FAILURE);
    }

    ndims = npars;

    branch_names = new TString[ndims];
    branch_values = new double[ndims];

    return;
}


// Load MCMC chain from file
void mcmcConvergence::LoadChain(TString mcmc_fname) {
    // Watch for bad input
    if (ndims < 1) {
        std::cout << "ERROR: branch arrays not initialized. Must specify a number" << std::endl;
        std::cout << "       of dimensions using SetDim() before loading MCMC chain." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load mcmc posteriors chain from file
    mcmc_file = new TFile(mcmc_fname.Data(), "READ");
    if (!mcmc_file->IsOpen()) {
        std::cout << "ERROR: File " << mcmc_fname.Data() << " could not be opened" << std::endl;
        exit(EXIT_FAILURE);
    }

    SetupBranches();
    FillArrays();

    return;
}


// Sets up the output file
void mcmcConvergence::PrepOutput(TString fname_out) {

    output_name = fname_out;
    output_file = new TFile(output_name.Data(), "RECREATE");
    output_file->cd();

    return;
}


// Sets memory addresses for input mcmc root file
void mcmcConvergence::SetupBranches() {
    // Prepare branches for loading
    std::cout << "Loading MCMC chain from " << mcmc_file->GetName() << std::endl;
    mcmc_chain = (TTree*)mcmc_file->Get("posteriors");
    TObjArray* branch_list = (TObjArray*)mcmc_chain->GetListOfBranches();
    int nbranches = branch_list->GetEntries();

    // Set branch addresses so we can dynamically load mcmc posteriors
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
        if (bname.CompareTo("lnl_current")==0) {
            mcmc_chain->SetBranchAddress("lnl_current", &lnl_current);
        }
    }

    return;
}


// Loads disk root data into memory for fast access
void mcmcConvergence::FillArrays() {

    std::cout << "Filling arrays for chain in " << mcmc_file->GetName() << std::endl;

    // Fill arrays with root data
    nsteps = mcmc_chain->GetEntries();
    lnls = new Float_t[nsteps];
    par_vals = new double*[nsteps];
    for (int istep = 0; istep < nsteps; ++istep) {
        par_vals[istep] = new double[ndims];
        // Fill with chain info
        mcmc_chain->GetEntry(istep);
        for (int idim = 0; idim < ndims; ++idim) {
            par_vals[istep][idim] = branch_values[idim];
        }
        lnls[istep] = lnl_current;
    }

    return;
}


// Fill Parameter trace plots
void mcmcConvergence::TraceParams() {
    // Check for bad input
    if (!mcmc_chain || ndims < 1) {
        std::cout << "ERROR: no mcmc file loaded for param traces" << std::endl;
        exit(EXIT_FAILURE);
    } else if (!output_file) {
        std::cout << "ERROR: Must set up output file before calling TraceParams" << std::endl;
        exit(EXIT_FAILURE);
    }
    TDirectory* trace_dir = output_file->mkdir("ParTraces");
    trace_dir->cd();

    // Initialize histograms
    param_trace_plots = new TH1D*[ndims+1]; // trace for [ndims] parameters + [1] lnl
    for (int idim = 0; idim < ndims; ++idim) {
        TString plt_name = branch_names[idim] + "_trace";
        param_trace_plots[idim] = new TH1D(plt_name.Data(), plt_name.Data(), nsteps, 0., nsteps);
        param_trace_plots[idim]->GetXaxis()->SetTitle("Step");
        param_trace_plots[idim]->GetYaxis()->SetTitle(branch_names[idim].Data());
    }
    param_trace_plots[ndims] = new TH1D("Lnl", "Lnl", nsteps, 0, nsteps);
    param_trace_plots[ndims]->GetXaxis()->SetTitle("Step");
    param_trace_plots[ndims]->GetYaxis()->SetTitle("Lnl");

    // Fill histograms
    for (int istep = 0; istep < nsteps; ++istep) {
        for (int idim = 0; idim < ndims; ++idim) {
            param_trace_plots[idim]->SetBinContent(istep+1, par_vals[istep][idim]);
        }
        param_trace_plots[ndims]->SetBinContent(istep+1, lnls[istep]);
    }

    for (int idim = 0; idim < ndims+1; ++idim) {
        param_trace_plots[idim]->Write();
    }

    return;
}


// Calculates Autocorrelation vector for a given lag value
void mcmcConvergence::CalcAutocorrelations(int nlag) {
    // Check for bad input
    if (!mcmc_chain || ndims < 1) {
        std::cout << "ERROR: no mcmc file loaded for autocorrelation calculation" << std::endl;
        exit(EXIT_FAILURE);
    } else if (nlag >= nsteps) {
        std::cout << "ERROR: nlag >= nsteps in the MCMC chain. Please use smaller nlag" << std::endl;
        exit(EXIT_FAILURE);
    } else if (!output_file) {
        std::cout << "ERROR: Must set up output file before calling CalcAutocorrelations" << std::endl;
        exit(EXIT_FAILURE);
    }

    TDirectory* autocorrelation_dir = output_file->mkdir("Autocorrelations");
    autocorrelation_dir->cd();

    autocorrelation_plots = new TH1D*[ndims];
    for (int idim = 0; idim < ndims; ++idim) {
        TString plt_name = branch_names[idim] + "_AC";
        autocorrelation_plots[idim] = new TH1D(plt_name.Data(), plt_name.Data(), nlag, 0, nlag);
        autocorrelation_plots[idim]->GetXaxis()->SetTitle("Lag");
        autocorrelation_plots[idim]->GetYaxis()->SetTitle("Autocorrelation function");
    }

    // Calculate at each lag value
    for (int ilag = 0; ilag < nlag; ++ilag) {

        // Prep variables to hold sums
        double* mult_sum = new double[ndims];
        double* par_sum = new double[ndims];
        double* lag_sum = new double[ndims];
        double* par_sum_sq = new double[ndims];
        double* lag_sum_sq = new double[ndims];
        // Initialize to 0
        for (int idim = 0; idim < ndims; ++idim) {
            mult_sum[idim] = 0;
            par_sum[idim] = 0;
            lag_sum[idim] = 0;
            par_sum_sq[idim] = 0;
            lag_sum_sq[idim] = 0;
        }

        // Loop through chain steps
        for (int istep = ilag; istep < nsteps; ++istep) {
            // Calculate for each parameter separately
            for (int idim = 0; idim < ndims; ++idim) {
                mult_sum[idim] += par_vals[istep][idim]*par_vals[istep-ilag][idim];
                par_sum[idim] += par_vals[istep][idim];
                lag_sum[idim] += par_vals[istep-ilag][idim];
                par_sum_sq[idim] += par_vals[istep][idim]*par_vals[istep][idim];
                lag_sum_sq[idim] += par_vals[istep-ilag][idim]*par_vals[istep-ilag][idim];
            }
        }

        // Finalize calculation
        double norm = (double)(nsteps-ilag);
        for (int idim = 0; idim < ndims; ++idim) {
            double numerator = norm*mult_sum[idim] - par_sum[idim]*lag_sum[idim];
            double denominator = (norm*par_sum_sq[idim] - par_sum[idim]*par_sum[idim]) \
                                 *(norm*lag_sum_sq[idim] - lag_sum[idim]*lag_sum[idim]);
            autocorrelation_plots[idim]->SetBinContent(ilag+1, numerator/TMath::Sqrt(denominator));
        }

        delete[] mult_sum;
        delete[] par_sum;
        delete[] lag_sum;
        delete[] par_sum_sq;
        delete[] lag_sum_sq;
    }

    for (int idim = 0; idim < ndims; ++idim) {
        autocorrelation_plots[idim]->Write();
    }

    return;
}


