/*
 * Module to run and manage mcmc
 */

#include "mcmc.h"


// Default constructor; sets everything to null or out-of-bounds
mcmc::mcmc() {
    SetDefault();
}


// Constructor using input configuration file
mcmc::mcmc(ConfigManager* configs) {
    // Start by setting everything out-of-bounds so we can tell if something is missing
    SetDefault();

    // Set necessary starting parameters
    if (configs->GetNPars() > 0) {
        npars = configs->GetNPars();
        // Setup par vectors if there is no input to go off
        if (configs->GetNewChainFlg()) {
            current_pars = new TVectorD(npars);
            current_pars->Zero();
            proposed_pars = new TVectorD(npars);
            proposed_pars->Zero();
        }
        // Identity matrix needed for proposal function
        identity_matrix = new TMatrixDSym(npars);
        for (int i = 0; i < npars; ++i) {
            (*identity_matrix)(i,i) = 1.;
            for (int j = i+1; j < npars; ++j) {
                (*identity_matrix)(i,j) = 0.;
            }
        }
        global_step_size = (2.38*2.38)/(Float_t)npars;
    }
    if (configs->GetEpsilon() > -1.) {
        epsilon = configs->GetEpsilon();
    }
    if (configs->GetNSteps() > -1) {
        nsteps = configs->GetNSteps();
    }
    if (configs->GetNBranch() > -1) {
        branch = configs->GetNBranch();
    }
    if (configs->GetRunNumber() > -1) {
        nrun_current = configs->GetRunNumber();
        if (configs->GetNewChainFlg() && nrun_current > 0) {
            std::cout << "ERROR: New chains should start at run number 0." << std::endl;
            std::cout << "       You input " << nrun_current << " for the current run number." << std::endl;
            std::cout << "       Check your config file." << std::endl;
            exit(EXIT_FAILURE);
        } else {
            nrun_previous = nrun_current - 1;
        }
    }
    // If any of the above is uninitialized, tell user to check config file
    if (epsilon < 0 || npars < 1 || nsteps < 0 || branch < 0 || nrun_current < 0) {
        std::cout << "ERROR: One or more critical MCMC variables uninitialized." << std::endl;
        std::cout << "       Please check your config file." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Check for I/O file naming
    if (configs->GetOutputDirectory().length() > 0) {
        output_dir = configs->GetOutputDirectory();
    }
    if (!configs->GetNewChainFlg() && configs->GetInputDirectory().length() > 0) {
        input_dir = configs->GetInputDirectory();
    }
    if (configs->GetMCMCFileBase().length() > 0) {
        if (output_dir.Length() > 0) {
            output_fname.Form("%s%s_branch%d_run%d.root",output_dir.Data(),configs->GetMCMCFileBase().c_str(),branch,nrun_current);
        }
        if (!configs->GetNewChainFlg() && input_dir.Length() > 0) {
            input_fname.Form("%s%s_branch%d_run%d.root",input_dir.Data(),configs->GetMCMCFileBase().c_str(),branch,nrun_previous);
            LoadPrevChain(); // Load the chain from this file
        }
    }

    // Prepare covariance matrices needed for step proposal and lnl evaluation
    InitCovMats(configs->GetCovmatFileBase(), configs->GetNewChainFlg());

    // If we're starting fresh, initialize the chain from scratch
    if (configs->GetNewChainFlg()) { InitNewChain(); }

}


// Sets default values of all variables to be default or out-of-bounds
void mcmc::SetDefault() {
    // File I/O
    input_dir = "";
    output_dir = "";
    input_fname = "";
    output_fname = "";
    input_file = NULL;
    output_file = NULL;
    mcmc_chain = NULL;

    // Covariance
    target_cov_fname = "";
    target_cov_file = NULL;
    target_covmat = NULL;
    proposal_cov_fname = "";
    proposal_cov_file = NULL;
    proposal_covmat = NULL;
    proposal_cholesky = NULL;
    identity_matrix = NULL;

    // Par vecs
    current_pars = NULL;
    proposed_pars = NULL;

    // Proposal function vars
    lnl_current = -1.;
    lnl_proposed = -1.;

    // MCMC run configs
    global_step_size = -1.;
    epsilon = -1.;
    npars = -1;
    nsteps = -1;
    nstep_current = 0;
    nrun_current = -1;
    nrun_previous = -1;
    branch = -1;
    naccepted = 0;

    // Setup RNG
    rng = new TRandom3();
    rng->SetSeed();

    return;
}


// Destructor
mcmc::~mcmc() {
    delete input_file;
    delete output_file;
    delete mcmc_chain;
    delete target_cov_file;
    delete target_covmat;
    delete proposal_cov_file;
    delete proposal_covmat;
    delete proposal_cholesky;
    delete current_pars;
    delete proposed_pars;
}


// Load an MCMC chain, reading from the end of a previous chain
void mcmc::LoadPrevChain() {

    std::cout << "Starting MCMC from previous chain in " << input_fname << std::endl;

    input_file = new TFile(input_fname.Data(),"READ");
    if (!input_file->IsOpen()) {
        std::cout << "ERROR: Invalid input file name " << input_fname << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load data from input file
    TTree* mcmc_chain_prev = (TTree*)input_file->Get("posteriors");
    TObjArray* branch_list = (TObjArray*)mcmc_chain_prev->GetListOfBranches();
    int nbranches = branch_list->GetEntries();
    TString branch_names[nbranches];
    double branch_vals[nbranches];
    int nstep_previous = -1;

    // Setup branches
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        TBranch* br = (TBranch*)branch_list->At(ibr);
        branch_names[ibr] = br->GetName();
        if (branch_names[ibr].CompareTo("nstep_current")) {
            mcmc_chain_prev->SetBranchAddress(branch_names[ibr], &nstep_previous);
        } else {
            mcmc_chain_prev->SetBranchAddress(branch_names[ibr], &branch_vals[ibr]);
        }
    }

    // Load last entry in previous chain for new starting point
    mcmc_chain_prev->GetEntry(mcmc_chain_prev->GetEntries()-1);

    // Set up the new MCMC chain
    mcmc_chain = new TTree("posteriors", "Posterior_Distributions");
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        // Load miscellaneous pars first
        if (branch_names[ibr].CompareTo("lnl_current")) {
            lnl_current = branch_vals[ibr];
            mcmc_chain->Branch("lnl_current", &lnl_current, "lnl_current/D");
            break;
        } else if (branch_names[ibr].CompareTo("nstep_current")) {
            nstep_current = nstep_previous + 1;
            mcmc_chain->Branch("nstep_current", &nstep_current, "nstep_current/I");
            break;
        }

        // If none of those, look through mcmc pars
        for (int ipar = 0; ipar < npars; ++ipar) {
            if (branch_names[ibr].CompareTo(Form("mcmc_par_%d",ipar))) {
                current_pars[ipar] = branch_vals[ibr];
                mcmc_chain->Branch(Form("mcmc_par_%d",ipar), &(*current_pars)(ipar), Form("mcmc_par_%d/D",ipar));
                break;
            }
        }
    }

    input_file->Close();

    return;
}


// Initializes a new MCMC chain
void mcmc::InitNewChain() {
    if (npars < 0) {
        std::cout << "ERROR: Forgot to set a number of parameters." << std::endl;
        std::cout << "       Check your config file." << std::endl;
        exit(EXIT_FAILURE);
    }

    mcmc_chain = new TTree("posteriors", "Posterior_Distributions");
    // Start with mcmc chain parameters
    for (int ipar = 0; ipar < npars; ++ipar) {
        mcmc_chain->Branch(Form("mcmc_par_%d",ipar), &(*current_pars)(ipar), Form("mcmc_par_%d/D",ipar));
    }
    // Setup other misc branches
    mcmc_chain->Branch("lnl_current", &lnl_current, "lnl_current/D");
    mcmc_chain->Branch("nstep_current", &nstep_current, "nstep_current/I");

    return;
}


void mcmc::InitCovMats(TString covmat_fname_base, bool new_chain) {

    // Check for input covariance matrix
    if (covmat_fname_base.Length() > 0) {
        if (input_dir.Length() > 0) {
            // Load proposal covmat if not a new chain
            if (!new_chain) {
                proposal_cov_fname.Form("%s%s_branch%d_run%d.root",input_dir.Data(),covmat_fname_base.Data(),branch,nrun_previous);
                proposal_cov_file = new TFile(proposal_cov_fname.Data());
                proposal_covmat = (TMatrixDSym*)proposal_cov_file->Get("cov_mat");
            }
            // Load target distribution covmat
            target_cov_fname.Form("%s%s_branch%d_target.root",input_dir.Data(),covmat_fname_base.Data(),branch);
            target_cov_file = new TFile(target_cov_fname.Data());
            if (!target_cov_file->IsOpen()) {
                std::cout << "ERROR: Invalid covmat file " << target_cov_fname << std::endl;
                exit(EXIT_FAILURE);
            }
            target_covmat = (TMatrixDSym*)target_cov_file->Get("cov_mat");
            (*target_covmat_inverted) = target_covmat->Invert();
            target_means = (TVectorD*)target_cov_file->Get("mean_vec");
        }
    }

    // Make sure we actually have a target covmat to fit to
    if (target_covmat == NULL) {
        std::cout << "ERROR: Target distribution not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Default proposal covmat is identity matrix
    if (new_chain) {
        proposal_covmat = (TMatrixDSym*)identity_matrix->Clone();
    }

    // Initialize proposal covmat using Haario et al
    (*identity_matrix) *= epsilon;
    (*proposal_covmat) += (*identity_matrix);
    (*proposal_covmat) *= global_step_size;
    (*proposal_cholesky) = GetCholDecomp();

    return;
}


// Returns the Cholesky Decomposition of proposal input covariance matrix
TMatrixD mcmc::GetCholDecomp() {
    TDecompChol* chol_decomp = new TDecompChol(*proposal_covmat);
    if (!chol_decomp->Decompose()) {
        std::cout << "ERROR: Cholesky Decomposition Failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
    TMatrixD chol_decomp_u = chol_decomp->GetU();
    return chol_decomp_u.T();
}


// Prepare the output file
void mcmc::PrepareOutput() {

    output_file = new TFile(output_fname.Data(),"RECREATE");
    if (!output_file->IsOpen()) {
        std::cout << "ERROR: Invalid output file name " << output_fname << std::endl;
        exit(EXIT_FAILURE);
    }

    return;
}


// Main function for running the MCMC
void mcmc::RunMCMC(double _lnl_current) {

    PrepareOutput();

    // If it's a new chain, get a random starting position
    if (_lnl_current < 0.) {
        ProposeStep();
        lnl_current = CalcPDF();
    } else {
        lnl_current = _lnl_current;
    }

    clock.Start();

    while (nstep_current < nstep_current+nsteps) {
        // Progress update
        if (nstep_current % (nsteps/10) == 0) {
            std::cout << "  MCMC on step " << nstep_current << " / " << nstep_current+nsteps << std::endl;
            std::cout << "    Accepted " << naccepted << " steps so far." << std::endl;
        }

        ProposeStep();
        lnl_proposed = CalcPDF();
        if (CheckIfAccepted()) { AcceptStep(); }
        else { RejectStep(); }

        ++nstep_current;
    }

    clock.Stop();
    SaveChain();

    return;
}


// Propose next step in the MCMC using the Cholesky Decomposition of some covariance matrix
void mcmc::ProposeStep() {

    // Start with a vector of random gaussian throws
    TVectorD* randomizer = new TVectorD(npars);
    for (int i = 0; i < npars; ++i) {
        (*randomizer)(i) = rng->Gaus(0.,1.);
    }

    // Multiply by cholesky decomposed covmat to get throw deltas
    (*proposed_pars) = (*proposal_cholesky) * (*randomizer);

    // Add to current parameter vector to get new step proposal
    (*proposed_pars) += (*current_pars);

    return;
}


// Calculate the PDF using the multivariate normal distribution
Float_t mcmc::CalcPDF() {
    Float_t PDF = 0.;

    // Calc difference vector
    TVectorD* diff_vec = proposed_pars;
    (*diff_vec) -= (*target_means);

    // Compute matrix sum to get total PDF
    //for (int i = 0; i < npars; ++i) {
    //    for (int j = 0; j < npars; ++j ) {
    //        PDF += (*diff_vec)(i) * (*target_covmat_inverted)(i,j) * (*diff_vec)(j);
    //    }
    //}

    // This can be computed more simply using the following:
    PDF = (*diff_vec) * ((*target_covmat_inverted) * (*diff_vec));

    return PDF;
}


// Check if the proposed MCMC step should be accepted using the Metropolis-Hastings algorithm
bool mcmc::CheckIfAccepted() {
    Double_t accept_prob = TMath::Min(1.,TMath::Exp(lnl_current-lnl_proposed));
    return rng->Rndm() < accept_prob;
}


// Accept the proposed step
void mcmc::AcceptStep() {
    // Update parameters and write output
    ++naccepted;
    lnl_current = lnl_proposed;
    for (int ipar = 0; ipar < npars; ++ipar) {
        (*current_pars)(ipar) = (*proposed_pars)(ipar);
    }

    mcmc_chain->Write();
    return;
}


// Reject the proposed step
void mcmc::RejectStep() {
    // Keep parameters where they are and just fill the output
    mcmc_chain->Write();
    return;
}


// Save final output and close output file
void mcmc::SaveChain() {

    // Warn if the MCMC has accepted no steps
    if (naccepted < 1) {
        std::cout << "\nWARNING: No steps accepted in the MCMC!" << std::endl;
    } else {
        std::cout << "\nChain of length " << nsteps << " took " << clock.RealTime()
                  << " seconds to complete. (" << clock.RealTime() / (double)nsteps
                  << " seconds/step)." << std::endl;
        std::cout << naccepted << " steps were accepted. ("
                  << (double)naccepted*100./(double)nsteps << " \%)" << std::endl;
    }

    // Save output and close up shop
    mcmc_chain->Write();
    output_file->Close();

    return;
}

