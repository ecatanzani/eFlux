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

void buildSignalEfficiency(
    const char* signal_eff_file,
    const char* output_file,
    const bool verbose) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        auto h_bdt_signal_not_passed        = GetHistoFromFile(signal_eff_file, "h_bdt_signal_not_passed", verbose);
        auto h_xtrl_tight_signal_not_passed = GetHistoFromFile(signal_eff_file, "h_xtrl_tight_signal_not_passed", verbose);
        auto h_xtrl_loose_signal_not_passed = GetHistoFromFile(signal_eff_file, "h_xtrl_loose_signal_not_passed", verbose);
        auto h_signal                       = GetHistoFromFile(signal_eff_file, "h_signal", verbose);

        h_bdt_signal_not_passed             ->Sumw2();
        h_xtrl_tight_signal_not_passed      ->Sumw2();
        h_xtrl_loose_signal_not_passed      ->Sumw2();
        h_signal                            ->Sumw2();

        auto h_bdt_signal_efficiency            = static_cast<TH1D*>(h_bdt_signal_not_passed->Clone("h_bdt_signal_efficiency"));
        auto h_xtrl_tight_signal_efficiency     = static_cast<TH1D*>(h_xtrl_tight_signal_not_passed->Clone("h_xtrl_tight_signal_efficiency"));
        auto h_xtrl_loose_signal_efficiency     = static_cast<TH1D*>(h_xtrl_loose_signal_not_passed->Clone("h_xtrl_loose_signal_efficiency"));

        h_bdt_signal_efficiency             ->Divide(h_signal.get());
        h_xtrl_tight_signal_efficiency      ->Divide(h_signal.get());
        h_xtrl_loose_signal_efficiency      ->Divide(h_signal.get());

        for (int bin=1; bin<=h_bdt_signal_efficiency->GetNbinsX(); ++bin) {
            h_bdt_signal_efficiency         ->SetBinContent(bin, 1 - h_bdt_signal_efficiency->GetBinContent(bin));
            h_xtrl_tight_signal_efficiency  ->SetBinContent(bin, 1 - h_xtrl_tight_signal_efficiency->GetBinContent(bin));
            h_xtrl_loose_signal_efficiency  ->SetBinContent(bin, 1 - h_xtrl_loose_signal_efficiency->GetBinContent(bin));
        }

        TFile* outfile = TFile::Open(output_file, "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]\n\n";
            exit(100);
        }
        
        h_bdt_signal_efficiency             ->Write();
        h_xtrl_tight_signal_efficiency      ->Write();
        h_xtrl_loose_signal_efficiency      ->Write();

        outfile->Close();

    }