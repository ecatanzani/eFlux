#include <iostream>
#include <memory>
#include <tuple>

#include "TH1D.h"
#include "TFile.h"

inline std::tuple<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> extractHistoFromFile(const char* file, const bool verbose) {
    
    TFile* infile {TFile::Open(file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << file << "]" << std::endl;
        exit(100);
    }

    std::shared_ptr<TH1D> h_proton_passed = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get("h_proton_passed")));
    std::shared_ptr<TH1D> h_proton_not_passed = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get("h_proton_not_passed")));
    
    return std::tuple<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>>(h_proton_passed, h_proton_not_passed);
}

void buildBackgrounndFraction(
    const char* input_file, 
    const char* output_file,
    const bool verbose) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        std::shared_ptr<TH1D> h_proton_passed, h_proton_not_passed;
        std::tie(h_proton_passed, h_proton_not_passed) = extractHistoFromFile(input_file, verbose);

        h_proton_passed->Sumw2();
        h_proton_not_passed->Sumw2();

        auto h_proton_fraction = static_cast<TH1D*>(h_proton_passed->Clone("h_proton_fraction"));
        h_proton_fraction->Divide(h_proton_not_passed.get());

        TFile* outfile = TFile::Open(output_file, "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]\n\n";
            exit(100);
        }

        h_proton_passed->Write();
        h_proton_not_passed->Write();
        h_proton_fraction->Write();

        outfile->Close();
    }