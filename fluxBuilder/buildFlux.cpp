#include <iostream>

#include "TKey.h"
#include "TH1D.h"
#include "TFile.h"

inline std::shared_ptr<TH1D> extractHistoFromFile(const char* file, const bool verbose) {
    TFile* infile {TFile::Open(file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key {nullptr};
    std::shared_ptr<TH1D> myhisto;
    while ((key=static_cast<TKey*>(nextkey())))  {
        TObject *obj {key->ReadObj()};
        if (obj->IsA()->InheritsFrom(TH1D::Class())) {
            myhisto = std::shared_ptr<TH1D>(static_cast<TH1D*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound TTree in input file [" << myhisto->GetName() << "] --> [" << file << "]";
    
    return myhisto;
}

void buildFlux(
    const char* e_counts_input_file,
    const char* e_acc_input_file,
    const double exposure_time,
    const char* output_file,
    const bool verbose) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        // Extract info from input files
        auto h_e_counts = extractHistoFromFile(e_counts_input_file, verbose);
        auto h_e_acc = extractHistoFromFile(e_acc_input_file, verbose);

        // Divide by the acceptance
        h_e_counts->Divide(h_e_acc.get());

        // Scale by the exposure-time
        h_e_counts->Scale(1/exposure_time);

        // divide by the energy bin width
        for (int bIdx {0}; bIdx<=h_e_counts->GetNbinsX(); ++bIdx) {
            if(h_e_counts->GetBinContent(bIdx)) {
                h_e_counts->SetBinContent(bIdx, static_cast<double>(h_e_counts->GetBinContent(bIdx))/h_e_counts->GetBinWidth(bIdx));
                h_e_counts->SetBinError(bIdx, static_cast<double>(h_e_counts->GetBinError(bIdx))/h_e_counts->GetBinWidth(bIdx));
            }
        }

        // Build flux multiplied by E^3
        auto h_e_counts_E3 = static_cast<TH1D*>(h_e_counts->Clone("h_e_counts_E3"));
        for (int bIdx {0}; bIdx<=h_e_counts_E3->GetNbinsX(); ++bIdx) {
            if(h_e_counts_E3->GetBinContent(bIdx)) {
                h_e_counts_E3->SetBinContent(bIdx, h_e_counts_E3->GetBinContent(bIdx)*pow(h_e_counts_E3->GetXaxis()->GetBinCenter(bIdx), 3));
                h_e_counts_E3->SetBinError(bIdx, h_e_counts_E3->GetBinError(bIdx)*pow(h_e_counts_E3->GetXaxis()->GetBinCenter(bIdx), 3));
            }
        }

        // Write output file
        TFile *outfile {TFile::Open(output_file, "RECREATE")};
        if (outfile->IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        h_e_counts->SetName("h_all_e_flux");
        h_e_counts->SetTitle("All-Electron DAMPE flux");
        h_e_counts->GetXaxis()->SetTitle("Energy [GeV]");
        h_e_counts->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

        h_e_counts_E3->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

        h_e_counts->Write();
        h_e_counts_E3->Write();

        outfile->Close();

        if (verbose)
            std::cout << "\n\nOutput file has been written ... [" << output_file << "]\n\n";
    }