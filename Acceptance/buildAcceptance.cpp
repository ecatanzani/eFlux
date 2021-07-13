#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"

const double generation_vertex_radius = 1.381976597885342;

void buildAcceptance(const char* full_acceptance_histos, const char* output_file = "acceptance.root")
{
    TFile* input_file = TFile::Open(full_acceptance_histos, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_acceptance_histos << "]\n\n";
        exit(100);
    }

    TH1D* h_gen = static_cast<TH1D*>(input_file->Get("h_gen"));
    TH1D* h_geometric = static_cast<TH1D*>(input_file->Get("h_geometric"));
    TH1D* h_bgo_fiducial = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial"));
    TH1D* h_all_cut = static_cast<TH1D*>(input_file->Get("h_all_cut"));

    h_gen->SetDirectory(0);
    h_geometric->SetDirectory(0);
    h_bgo_fiducial->SetDirectory(0);
    h_all_cut->SetDirectory(0);

    input_file->Close();

    TH1D* h_acc_geometric = static_cast<TH1D*>(h_geometric->Clone("h_acc_geometric"));
    TH1D* h_acc_bgo_fiducial = static_cast<TH1D*>(h_bgo_fiducial->Clone("h_acc_bgo_fiducial"));
    TH1D* h_acc_all_cut = static_cast<TH1D*>(h_all_cut->Clone("h_acc_all_cut"));

    h_acc_geometric->Divide(h_gen);
    h_bgo_fiducial->Divide(h_gen);
    h_all_cut->Divide(h_gen);

    double genSurface = 4 * TMath::Pi() * pow(generation_vertex_radius, 2) / 2;
    double scaleFactor = TMath::Pi() * genSurface;

    h_acc_geometric->Scale(scaleFactor);
    h_bgo_fiducial->Scale(scaleFactor);
    h_all_cut->Scale(scaleFactor);

    TFile* outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output acceptance file: [" << output_file << "]\n\n";
        exit(100);
    }

    h_acc_geometric->Write();
    h_bgo_fiducial->Write();
    h_all_cut->Write();

    outfile->Close();

    std::cout << "\n\nOutput file has been written: [" << output_file << "]\n\n";
}