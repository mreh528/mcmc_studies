/*
 * Module to run and manage mcmc
 */

#include "mcmc.h"


// Default constructor; sets everything to null or out-of-bounds
mcmc::mcmc() {
    SetDefault();
}


// Constructor using input configuration file
mcmc::mcmc(ConfigManager* configs, int _nrun, int _nbranch) {
    // Start by setting everything out-of-bounds so we can tell if something is missing
    SetDefault();

    // Read config file and set necessary variables
    ReadConfigs(configs, _nrun, _nbranch);

    // Prepare covariance matrices needed for step proposal and lnl evaluation
    InitCovMats();

    // If we're starting fresh, initialize the chain from scratch
    if (new_chain) { InitNewChain(); }
    else { LoadPrevChain(); }

}


// Destructor
mcmc::~mcmc() {
    delete input_chain_file;
    delete output_chain_file;
    delete mcmc_chain;
    delete target_cov_file;
    delete target_covmat;
    delete target_covmat_inverted;
    delete proposal_cov_file;
    delete proposal_covmat;
    delete proposal_cholesky;
    delete current_pars;
    delete proposed_pars;
}


// Sets default values of all variables to be default or out-of-bounds
void mcmc::SetDefault() {
    std::cout << "Setting default values..." << std::endl;

    // File I/O
    chain_dir = "";
    input_chain_fname = "";
    output_chain_fname = "";
    input_chain_file = NULL;
    output_chain_file = NULL;
    mcmc_chain = NULL;

    // Covariance
    target_cov_dir = "";
    target_cov_fname = "";
    target_cov_file = NULL;
    target_covmat = NULL;
    proposal_cov_dir = "";
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
    epsilon = 0.;
    npars = -1;
    nsteps = -1;
    nstep_current = 0;
    nstep_start = 0;
    nrun_current = -1;
    nrun_previous = -1;
    branch = -1;
    naccepted = 0;
    new_chain = true;

    // Setup RNG
    rng = new TRandom3();
    rng->SetSeed();

    return;
}


// Constructor helper to read config file and set needed variables
void mcmc::ReadConfigs(ConfigManager* configs, int _nrun, int _nbranch) {
    std::cout << "Reading configs..." << std::endl;

    // Check for bad input
    if (_nrun < 0 || _nbranch < 0) {
        std::cout << "ERROR: Branch number or run number not specified.\n"
                  << "       These should be specified at the command line\n"
                  << "       using \'-r\' and \'-b\'" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Set base configs
    branch = _nbranch;
    nrun_current = _nrun;
    nrun_previous = nrun_current - 1;
    new_chain = _nrun == 0;

    // Set necessary starting parameters
    if (configs->GetNPars() > 0) {
        npars = configs->GetNPars();
        // Setup things that depend on npars
        current_pars = new TVectorD(npars);
        proposed_pars = new TVectorD(npars);
        current_pars->Zero();
        proposed_pars->Zero();
        global_step_size = (2.38*2.38)/(Float_t)npars;
    }
    // Epsilon for regularization term
    if (configs->GetEpsilon() > 0.) {
        epsilon = configs->GetEpsilon();
    }
    // Number of steps in the chain
    if (configs->GetNSteps() > -1) {
        nsteps = configs->GetNSteps();
    }
    // If any of the above is uninitialized, tell user to check config file
    if (epsilon < 0. || npars < 1 || nsteps < 0) {
        std::cout << "ERROR: One or more critical MCMC variables uninitialized." << std::endl;
        std::cout << "       Please check your config file." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Check for I/O file naming
    // First check for input proposal covmat
    if (!new_chain) {
        if (configs->GetProposalCovDir().length() > 0 &&
            configs->GetProposalCovFileBase().length() > 0) {
            proposal_cov_dir = configs->GetProposalCovDir();
            proposal_cov_fname.Form("%s%s_npars%d_branch%d_run%d.root",\
                                    proposal_cov_dir.Data(),\
                                    configs->GetProposalCovFileBase().c_str(),\
                                    npars, branch, nrun_previous);
        } else {
            std::cout << "ERROR: Proposal cov not specified for continuing chain" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    // Check for new MCMC file output & previous chain input
    if (configs->GetMCMCFileBase().length() > 0 &&
        configs->GetChainDir().length() > 0) {
        chain_dir = configs->GetChainDir();
        output_chain_fname.Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                chain_dir.Data(),\
                                configs->GetMCMCFileBase().c_str(),\
                                nsteps, npars, branch, nrun_current);
        if (!new_chain) {
            input_chain_fname.Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                                   chain_dir.Data(),\
                                   configs->GetMCMCFileBase().c_str(),\
                                   nsteps, npars, branch, nrun_previous);
        }
    } else {
        std::cout << "ERROR: No output file specified" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Make sure there is a target to fit to
    if (configs->GetTargetCovDir().length() > 0 &&
        configs->GetTargetCovFileBase().length() > 0) {
        target_cov_dir = configs->GetTargetCovDir();
        target_cov_fname.Form("%s%s_npars%d_branch%d_target.root",
                              target_cov_dir.Data(),\
                              configs->GetTargetCovFileBase().c_str(),
                              npars, branch);
    } else {
        std::cout << "ERROR: No target distribution specified" << std::endl;
        exit(EXIT_FAILURE);
    }

    return;
}


// Load an MCMC chain, reading from the end of a previous chain
void mcmc::LoadPrevChain() {

    std::cout << "Starting MCMC from previous chain in " << input_chain_fname << std::endl;

    input_chain_file = new TFile(input_chain_fname.Data(), "READ");
    if (!input_chain_file->IsOpen()) {
        std::cout << "ERROR: Invalid input file name " << input_chain_fname << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load data from input file
    TTree* mcmc_chain_prev = (TTree*)input_chain_file->Get("posteriors");
    TObjArray* branch_list = (TObjArray*)mcmc_chain_prev->GetListOfBranches();
    int nbranches = branch_list->GetEntries();
    TString branch_names[nbranches];
    double branch_vals[nbranches];
    Float_t lnl_in = -1.;
    int nstep_previous = -1;

    // Setup branches
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        TBranch* br = (TBranch*)branch_list->At(ibr);
        branch_names[ibr] = br->GetName();
        if (branch_names[ibr].CompareTo("nstep_current")==0) {
            mcmc_chain_prev->SetBranchAddress(branch_names[ibr], &nstep_previous);
        } else if (branch_names[ibr].CompareTo("lnl_current")==0) {
            mcmc_chain_prev->SetBranchAddress(branch_names[ibr], &lnl_in);
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
        std::cout << "Loading branch " << branch_names[ibr] << std::endl;
        if (branch_names[ibr].CompareTo("lnl_current")==0) {
            lnl_current = lnl_in;
            mcmc_chain->Branch("lnl_current", &lnl_current, "lnl_current/F");
            continue;
        }
        if (branch_names[ibr].CompareTo("nstep_current")==0) {
            nstep_current = nstep_previous + 1;
            nstep_start = nstep_current;
            mcmc_chain->Branch("nstep_current", &nstep_current, "nstep_current/I");
            continue;
        }

        // If none of those, look through mcmc pars
        for (int ipar = 0; ipar < npars; ++ipar) {
            if (branch_names[ibr].CompareTo(Form("mcmc_par_%d",ipar))==0) {
                (*current_pars)(ipar) = branch_vals[ibr];
                mcmc_chain->Branch(branch_names[ibr].Data(), &(*current_pars)(ipar), Form("mcmc_par_%d/D",ipar));
                break;
            }
        }
    }

    std::cout << "\nStarting parameter values:" << std::endl;
    TObjArray* brs = (TObjArray*)mcmc_chain->GetListOfBranches();
    for (int ipar = 0; ipar < npars; ++ipar) {
        TBranch* br = (TBranch*)brs->At(ipar);
        std::cout << br->GetName() << " = " << (*current_pars)(ipar) << std::endl;
    }

    return;
}


// Initializes a new MCMC chain
void mcmc::InitNewChain() {
    std::cout << "Initializing new MCMC branches.." << std::endl;

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
    mcmc_chain->Branch("lnl_current", &lnl_current, "lnl_current/F");
    mcmc_chain->Branch("nstep_current", &nstep_current, "nstep_current/I");

    // Find a random starting position
    ProposeStep();
    CalcPDF();
    lnl_current = lnl_proposed;

    std::cout << "\nStarting parameter values:" << std::endl;
    TObjArray* brs = (TObjArray*)mcmc_chain->GetListOfBranches();
    for (int ipar = 0; ipar < npars; ++ipar) {
        TBranch* br = (TBranch*)brs->At(ipar);
        std::cout << br->GetName() << " = " << (*current_pars)(ipar) << std::endl;
    }

    return;
}


// Initialize proposal and target distributions
void mcmc::InitCovMats() {
    std::cout << "Initializing target and proposal covariance matrices..." << std::endl;

    // Identity matrix needed for proposal function
    identity_matrix = new TMatrixDSym(npars);
    for (int i = 0; i < npars; ++i) {
        (*identity_matrix)(i,i) = epsilon; // pre-apply epsilon scaling factor
        for (int j = i+1; j < npars; ++j) {
            (*identity_matrix)(i,j) = 0.;
        }
    }

    // Check for input covariance matrix
    if (proposal_cov_fname.Length() > 0 && !new_chain) {
        proposal_cov_file = new TFile(proposal_cov_fname.Data(), "READ");
        proposal_covmat = (TMatrixDSym*)proposal_cov_file->Get("cov_mat");
        std::cout << "Got Proposal covariance matrix from "
                  << proposal_cov_file->GetName() << std::endl;
    } else {
        proposal_covmat = (TMatrixDSym*)identity_matrix->Clone();
        proposal_covmat->Zero();
    }

    // Load target distribution covmat
    if (target_cov_fname.Length() > 0) {
        target_cov_file = new TFile(target_cov_fname.Data(), "READ");
        if (!target_cov_file->IsOpen()) {
            std::cout << "ERROR: Invalid covmat file " << target_cov_fname << std::endl;
            exit(EXIT_FAILURE);
        }
        target_covmat = (TMatrixDSym*)target_cov_file->Get("cov_mat");
        target_covmat_inverted = (TMatrixDSym*)target_covmat->Clone();
        target_covmat_inverted->Invert();
        target_means = (TVectorD*)target_cov_file->Get("mean_vec");
        std::cout << "Got Target distribution from "
                  << target_cov_file->GetName() << std::endl;
    }

    // Make sure we actually have a target covmat to fit to
    if (target_covmat == NULL) {
        std::cout << "ERROR: Target distribution not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Initialize proposal covmat using Haario et al
    (*proposal_covmat) += (*identity_matrix);
    (*proposal_covmat) *= global_step_size;
    GetCholDecomp();

    return;
}


// Computes the Cholesky Decomposition of proposal covariance matrix
void mcmc::GetCholDecomp() {
    TDecompChol* chol_decomp = new TDecompChol(*proposal_covmat);
    if (!chol_decomp->Decompose()) {
        std::cout << "ERROR: Cholesky Decomposition Failed!" << std::endl;
        exit(EXIT_FAILURE);
    }

    proposal_cholesky = new TMatrixD(chol_decomp->GetU());
    proposal_cholesky->T();
    return;
}


// Prepare the output file
void mcmc::PrepareOutput() {

    output_chain_file = new TFile(output_chain_fname.Data(), "RECREATE");
    if (!output_chain_file->IsOpen()) {
        std::cout << "ERROR: Invalid output file name " << output_chain_fname << std::endl;
        exit(EXIT_FAILURE);
    }

    return;
}


// Main function for running the MCMC
void mcmc::RunMCMC() {

    // Set up file I/O
    PrepareOutput();

    // Time the mcmc
    clock.Start();

    while (nstep_current < nstep_start+nsteps) {
        // Progress update
        if (nstep_current % (nsteps/10) == 0) {
            std::cout << "  MCMC on step " << nstep_current
                      << " / " << nstep_start+nsteps << std::endl;
            std::cout << "    Accepted " << naccepted << " steps so far." << std::endl;
        }

        ProposeStep();
        CalcPDF();
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
void mcmc::CalcPDF() {
    Float_t PDF = 0.;

    // Calc difference vector
    TVectorD* diff_vec = (TVectorD*)proposed_pars->Clone();
    (*diff_vec) -= (*target_means);

    // Compute matrix sum to get total PDF
    //for (int i = 0; i < npars; ++i) {
    //    for (int j = 0; j < npars; ++j ) {
    //        PDF += (*diff_vec)(i) * (*target_covmat_inverted)(i,j) * (*diff_vec)(j);
    //    }
    //}

    // This can be computed more simply using the following:
    PDF = (*diff_vec) * ((*target_covmat_inverted) * (*diff_vec));

    lnl_proposed = PDF;
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

    mcmc_chain->Fill();
    return;
}


// Reject the proposed step
void mcmc::RejectStep() {
    // Keep parameters where they are and just fill the output
    mcmc_chain->Fill();
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
        std::cout << "Saving output to " << output_chain_file->GetName() << std::endl;
    }

    // Save output and close up shop
    mcmc_chain->Write();
    output_chain_file->Close();

    return;
}

