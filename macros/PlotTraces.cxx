/*
 * Macro for plotting parameter and lnl traces from the MCMC
 */

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "CommandHandler.h"
#include "ConfigManager.h"

#include <iostream>
using namespace std;


int main(int argc, char* argv[]) {

    CommandHandler handler(argc, argv);
    ConfigManager* configs = new ConfigManager(handler.fname_config.Data());

    if (handler.GetRunNumber() < 0) {
        std::cout << "ERROR: Run number not set at cmd line" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Load MCMC file
    TString fname = Form("%s%s_nsteps%d_npars%d_branch%d_run%d.root",\
                         configs->GetChainDir().Data(),\
                         configs->GetMCMCFileBase().Data(),\
                         configs->GetNSteps(),\
                         configs->GetNPars(),\
                         configs->GetBranchNumber(),\
                         handler.GetRunNumber());
    std::cout << "Reading MCMC chain from " << fname.Data() << std::endl;
    TFile* fin = new TFile(fname.Data(), "READ");
    TTree* mcmc_chain = (TTree*)fin->Get("posteriors");

    // Set up branch objects
    TObjArray* branch_list = (TObjArray*)mcmc_chain->GetListOfBranches();
    int nbranches = branch_list->GetEntries();
    TString branch_names[nbranches];
    double branch_vals[nbranches];
    Float_t lnl_current = -1.;
    int nstep = -1;

    // Load branches to memory
    std::cout << "Loading branches to memory..." << std::endl;
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        TBranch* br = (TBranch*)branch_list->At(ibr);
        branch_names[ibr] = br->GetName();
        if (branch_names[ibr].CompareTo("nstep_current")==0) {
            mcmc_chain->SetBranchAddress(branch_names[ibr], &nstep);
        } else if (branch_names[ibr].CompareTo("lnl_current")==0) {
            mcmc_chain->SetBranchAddress(branch_names[ibr], &lnl_current);
        } else {
            mcmc_chain->SetBranchAddress(branch_names[ibr], &branch_vals[ibr]);
        }
    }

    // Set up graphs
    int nsteps = mcmc_chain->GetEntries();
    TGraph* graphs[nbranches];
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        if (branch_names[ibr].CompareTo("nstep_current")) { // Skips nstep
            graphs[ibr] = new TGraph(nsteps);
        }
    }

    // Fill up all the graphs
    std::cout << "Filling Graphs..." << std::endl;
    for (int istep = 0; istep < nsteps; ++istep) {
        mcmc_chain->GetEntry(istep);
        for (int ibr = 0; ibr < nbranches; ++ibr) {
            if (branch_names[ibr].CompareTo("nstep_current")==0) { // Skips nstep
                continue;
            } else if (branch_names[ibr].CompareTo("lnl_current")==0) {
                graphs[ibr]->SetPoint(istep, (double)nstep, (double)lnl_current);
            } else {
                graphs[ibr]->SetPoint(istep, (double)nstep, branch_vals[ibr]);
            }
        }
    }

    // Draw graphs
    std::cout << "Drawing Graphs..." << std::endl;
    for (int ibr = 0; ibr < nbranches; ++ibr) {
        if (branch_names[ibr].CompareTo("nstep_current")) { // Skips nstep
            TCanvas c1("c1","c1",1200,800);
            c1.cd();
            graphs[ibr]->SetLineColor(4);
            graphs[ibr]->SetTitle(Form("%s_trace",branch_names[ibr].Data()));
            graphs[ibr]->Draw();
            //if (branch_names[ibr].CompareTo("lnl_current")==0) {
            //    c1.SetLogy();
            //}
            c1.Print(Form("%s%s_trace_npars%d_branch%d_run%d.jpg",\
                          configs->GetPlotDir().Data(),\
                          branch_names[ibr].Data(),\
                          configs->GetNPars(),\
                          configs->GetBranchNumber(),\
                          handler.GetRunNumber()));
        }
    }

    fin->Close();
    return 0;
}

