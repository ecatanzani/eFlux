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
    const char* bdt_background_fraction,
    const char* xtrl_tight_background_fraction,
    const char* xtrl_loose_background_fraction,
    const char* output_file,
    const bool verbose) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        auto h_bdt_proton_fraction              = GetHistoFromFile(bdt_background_fraction, "h_bdt_proton_fraction", verbose);
        auto h_xtrl_tight_proton_fraction       = GetHistoFromFile(xtrl_tight_background_fraction, "h_xtrl_tight_proton_fraction", verbose);
        auto h_xtrl_loose_proton_fraction       = GetHistoFromFile(xtrl_loose_background_fraction, "h_xtrl_loose_proton_fraction", verbose);
        auto h_bdt_proton_not_passed            = GetHistoFromFile(rejected_proton_data, "h_proton_not_passed", verbose);
        auto h_xtrl_proton_not_passed           = GetHistoFromFile(rejected_proton_data, "h_xtrl_proton_not_passed", verbose);
        
        auto h_bdt_background                   = static_cast<TH1D*>(h_bdt_proton_not_passed->Clone("h_bdt_background"));
        auto h_xtrl_tight_background            = static_cast<TH1D*>(h_xtrl_proton_not_passed->Clone("h_xtrl_tight_background"));
        auto h_xtrl_loose_background            = static_cast<TH1D*>(h_xtrl_proton_not_passed->Clone("h_xtrl_loose_background"));

        h_bdt_background                        ->Multiply(h_bdt_proton_fraction.get());
        h_xtrl_tight_background                 ->Multiply(h_xtrl_tight_proton_fraction.get());
        h_xtrl_loose_background                 ->Multiply(h_xtrl_loose_proton_fraction.get());

        TFile* outfile = TFile::Open(output_file, "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]\n\n";
            exit(100);
        }
        
        h_bdt_background                        ->Write();
        h_xtrl_tight_background                 ->Write();
        h_xtrl_loose_background                 ->Write();

        outfile->Close();

    }