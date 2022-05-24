#include "signal.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>
#include <vector>

#include "TKey.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <ROOT/RDataFrame.hxx>

inline std::vector<double> extract_bdt_cut_values(const char* ext_tree_path, const unsigned int nenergy_bins) {
    std::vector<double> bdt_cut_values (nenergy_bins, -1);

    auto get_tree_name = [](const char* file) -> const char* {
        TFile* input_file = TFile::Open(file, "READ");
        if (!input_file->IsOpen()) {
            std::cerr << "\n\nError reading input file [" << file << "]\n\n";
            exit(100);
        }
        std::string tree_name;
        for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
            auto key = dynamic_cast<TKey*>(keyAsObject);
            if (!strcmp(key->GetClassName(), "TTree"))
                tree_name = static_cast<std::string>(key->GetName());
        }
        input_file->Close();
        return tree_name.c_str();
    };

    // Read external TFile
    TFile input_file(ext_tree_path, "READ");
    if (!input_file.IsOpen()) {
        std::cerr << "\n\nError reading input ROOT file [" << ext_tree_path << "]\n\n";
        exit(100);
    }

    // Extract the TTree
    std::unique_ptr<TTreeReader> bdt_reader {std::make_unique<TTreeReader>(get_tree_name(ext_tree_path), &input_file)};
    
    /*
    TTreeReaderValue<double> f_ec_X(*bdt_reader, "f_ec_X");
    TTreeReaderValue<double> f_ec_Y(*bdt_reader, "f_ec_Y");
    TTreeReaderValue<double> f_ec_err(*bdt_reader, "f_ec_err");
    */
    TTreeReaderValue<double> f_ec_b_sub_X(*bdt_reader, "f_ec_b_sub_X");
    /*
    TTreeReaderValue<double> f_ec_b_sub_Y(*bdt_reader, "f_ec_b_sub_Y");
    TTreeReaderValue<double> f_ec_b_sub_err(*bdt_reader, "f_ec_b_sub_err");
    */

    unsigned int energy_bin_idx {0};

    while (bdt_reader->Next()) {
        bdt_cut_values[energy_bin_idx] = *f_ec_b_sub_X;
        ++energy_bin_idx;
    }

    return bdt_cut_values;
}

void signal_selection(in_args input_args) {

    std::shared_ptr<energy_config> config = std::make_shared<energy_config>(input_args.energy_config_file);
    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.verbose);

    auto energy_binning = config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Extract external BDT cut values
    auto bdt_cuts = extract_bdt_cut_values(input_args.bdt_cut_values, energy_nbins);

    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evt_parser->GetEvtTree());

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nTotal events: " << *(_data_fr.Count());
    std::cout << "\n\n***************************";

    if (input_args.verbose) std::cout << "\n\nAnlysis running...\n\n";

    auto get_bdt_cut = [&bdt_cuts] (const int energy_bin) -> double {
        return bdt_cuts[energy_bin-1];
    };

    auto xtrl_loose_cut = [] (const double raw_energy_gev) -> double {
        return 3*pow(log10(raw_energy_gev/200), 2) + 12;
    };

    auto h_bdt_selection = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("bdt_cut", get_bdt_cut, {"energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value > tmva_cut; }, {"tmva_classifier", "bdt_cut"})
                                    .Histo1D<double>({"h_bdt_selection", "BDT event selection; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev");

    auto h_xtrl_tight_selection = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Filter("xtrl<8.5")
                                    .Histo1D<double>({"h_xtrl_tight_selection", "XTRL tight event selection; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev");

    auto h_xtrl_loose_selection = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("raw_energy_gev", "energy * 0.001")
                                    .Define("xtrl_loose_cut", xtrl_loose_cut, {"raw_energy_gev"})
                                    .Filter("xtrl<xtrl_loose_cut")
                                    .Histo1D<double>({"h_xtrl_loose_selection", "XTRL loose event selection; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev");

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen()) {
        std::cerr << "\n\nError writing output file... [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_bdt_selection         ->Write();
    h_xtrl_tight_selection  ->Write();
    h_xtrl_loose_selection  ->Write();

    outfile->Close();
}