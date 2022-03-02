#include <tuple>
#include <memory>
#include <iostream>

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

void buildBackgrounndFraction(
    const char* input_file, 
    const char* output_file,
    const bool verbose) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        auto h_bdt_background_passed            = GetHistoFromFile(input_file, "h_bdt_proton_passed", verbose);
        auto h_xtrl_tight_background_passed     = GetHistoFromFile(input_file, "h_xtrl_tight_proton_passed", verbose);
        auto h_xtrl_loose_background_passed     = GetHistoFromFile(input_file, "h_xtrl_loose_proton_passed", verbose);
        auto h_bdt_background_not_passed        = GetHistoFromFile(input_file, "h_bdt_proton_not_passed", verbose);
        auto h_xtrl_background_not_passed       = GetHistoFromFile(input_file, "h_xtrl_proton_not_passed", verbose);
        

        h_bdt_background_passed                 ->Sumw2();
        h_xtrl_tight_background_passed          ->Sumw2();
        h_xtrl_loose_background_passed          ->Sumw2();
        h_bdt_background_not_passed             ->Sumw2();
        h_xtrl_background_not_passed            ->Sumw2();
        

        auto h_bdt_proton_fraction              = static_cast<TH1D*>(h_bdt_background_passed->Clone("h_bdt_proton_fraction"));
        auto h_xtrl_tight_proton_fraction       = static_cast<TH1D*>(h_xtrl_tight_background_passed->Clone("h_xtrl_tight_proton_fraction"));
        auto h_xtrl_loose_proton_fraction       = static_cast<TH1D*>(h_xtrl_loose_background_passed->Clone("h_xtrl_loose_proton_fraction"));

        h_bdt_proton_fraction                   ->Divide(h_bdt_background_not_passed.get());
        h_xtrl_tight_proton_fraction            ->Divide(h_xtrl_background_not_passed.get());
        h_xtrl_loose_proton_fraction            ->Divide(h_xtrl_background_not_passed.get());
        
        TFile* outfile = TFile::Open(output_file, "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]\n\n";
            exit(100);
        }

        h_bdt_proton_fraction->Write();
        h_xtrl_tight_proton_fraction->Write();
        h_xtrl_loose_proton_fraction->Write();

        outfile->Close();
    }