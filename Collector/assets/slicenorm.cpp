#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"

#include <iostream>

void HistoSliceNormalization(const char* input_file, const char* histo_path, const char* normalized_histo_path) {
    auto slicenorm = [](TH2D* histo) { 
        for (int bidx=1; bidx<=histo->GetNbinsX(); ++bidx) {
		    auto scale_factor = 1/histo->Integral(bidx, bidx, 1, histo->GetNbinsY());
            if (!(TMath::IsNaN(scale_factor) || !(TMath::Finite(scale_factor)))) {
                std::cout << "\nBin X: " << bidx << " scale factor: " << scale_factor;
                for (int bidy=1; bidy<=histo->GetNbinsY(); ++bidy) {
                    histo->SetBinContent(bidx, bidy, histo->GetBinContent(bidx, bidy)*scale_factor);
                    histo->SetBinError(bidx, bidy, histo->GetBinError(bidx, bidy)*scale_factor);
                }
            }
	    }   
	    return histo; 
    };


    TFile *_input = TFile::Open(input_file, "READ");
    if (_input->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << input_file << "]\n\n";
        exit(100);
    }

    // Simu histos
    auto histo = slicenorm(static_cast<TH2D*>(_input->Get(histo_path))); 
    histo->SetDirectory(0);

    TFile *_output = TFile::Open(normalized_histo_path, "RECREATE");
    if (_output->IsZombie())
    {
        std::cerr << "\n\nError writing output file: [" << normalized_histo_path << "]\n\n";
        exit(100);
    }

    histo->Write();

    _output->Close();

}