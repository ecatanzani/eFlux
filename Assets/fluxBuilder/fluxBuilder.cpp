#include "TFile.h"
#include "TH1D.h"

#include <iostream>

void fluxBuilder(const char* electronMC, const char* protonMC, const char* data, const double livetime)
{
    // Extract electron acceptance
    TFile* electron_mc_file = TFile::Open(electronMC, "READ");
    if (electron_mc_file->IsZombie())
    {
        std::cerr << "\n\nError reading input electron MC reco file: " << electronMC << std::endl;
        exit(100);
    }

    TH1D* electron_acceptance = static_cast<TH1D*>(electron_mc_file->Get("Acceptance_histos/h_acceptance_all_cut"));
    electron_acceptance->SetDirectory(0);

    electron_mc_file->Close();

    // Extract proton contamination
    TFile* proton_mc_file = TFile::Open(protonMC, "READ");
    if (proton_mc_file->IsZombie())
    {
        std::cerr << "\n\nError reading input proton MC reco file: " << protonMC << std::endl;
        exit(100);
    }

    TH1D* proton_background_fraction = static_cast<TH1D*>(proton_mc_file->Get("mc_ancillary/proton_background_ratio"));
    proton_background_fraction->SetDirectory(0);

    proton_mc_file->Close();

    // Read data file
    TFile* data_file = TFile::Open(data, "READ");
    if (data_file->IsZombie())
    {
        std::cerr << "\n\nError reading input data reco file: " << data << std::endl;
        exit(100);
    }

    TH1D* h_electron_data_all_cut = static_cast<TH1D*>(data_file->Get("h_all_cut"));
    TH1D* h_electron_data_over_xtrl_cut = static_cast<TH1D*>(data_file->Get("mc_ancillary/h_background_over_xtrl_cut"));

    h_electron_data_all_cut->SetDirectory(0);
    h_electron_data_over_xtrl_cut->SetDirectory(0);

    data_file->Close();

    // Compute background protons
    TH1D* backgroung_protons = static_cast<TH1D*>(proton_background_fraction->Clone("backgroung_protons"));
    backgroung_protons->Multiply(h_electron_data_over_xtrl_cut);

    // Compute dN/dE
    TH1D* dNdE = static_cast<TH1D*>(h_electron_data_all_cut->Clone("dNdE"));
    dNdE->Add(backgroung_protons, -1);

    // Rescale by energy bin width
    for (int bIdx=1; bIdx <= dNdE->GetNbinsX(); ++bIdx)
        dNdE->SetBinContent(bIdx, dNdE->GetBinContent(bIdx)/dNdE->GetBinWidth(bIdx));

    // Build flux
    TH1D* all_electron_flux = static_cast<TH1D*>(dNdE->Clone("all_electron_flux"));
    all_electron_flux->Divide(electron_acceptance);
    all_electron_flux->Scale(1/livetime);

    TFile* all_electron_flux_file = TFile::Open("all_electron_flux.root", "RECREATE");
    if (all_electron_flux_file->IsZombie())
    {
        std::cerr << "\n\nError writing final flux ROOT file" << std::endl;
        exit(100);
    }

    electron_acceptance->Write();
    proton_background_fraction->Write();
    h_electron_data_all_cut->Write();
    h_electron_data_over_xtrl_cut->Write();
    backgroung_protons->Write();
    dNdE->Write();
    all_electron_flux->Write();

    all_electron_flux_file->Close();
}