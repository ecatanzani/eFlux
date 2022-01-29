#include <iostream>
#include <memory>
#include <tuple>

#include "TH1D.h"
#include "TFile.h"

inline std::shared_ptr<TH1D> GetHistoFromFile(const char* path, const char* histo_name, const bool verbose) {
    TFile* infile {TFile::Open(path, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << path << "]" << std::endl;
        exit(100);
    }

    std::shared_ptr<TH1D> histo = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get(histo_name)));

    if (verbose)
        std::cout << "\nFound TTree in input file [" << histo->GetName() << "] --> [" << path << "]";

    return histo;
}

void buildBackgrounndEstimation(
    const char* rejected_proton_data,
    const char* background_proton_fraction,
    const char* output_file,
    const bool verbose) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        auto h_proton_fraction = GetHistoFromFile(background_proton_fraction, "h_proton_fraction", verbose);
        auto h_proton_rejected = GetHistoFromFile(rejected_proton_data, "h_proton_not_passed", verbose);

        auto h_background = static_cast<TH1D*>(h_proton_rejected->Clone("h_background"));
        h_background->Multiply(h_proton_fraction.get());

        TFile* outfile = TFile::Open(output_file, "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]\n\n";
            exit(100);
        }
        
        h_proton_fraction->Write();
        h_proton_rejected->Write();
        h_background->Write();

        outfile->Close();

    }