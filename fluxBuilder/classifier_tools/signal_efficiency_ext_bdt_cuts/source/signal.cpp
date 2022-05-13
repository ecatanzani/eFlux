#include "signal.h"
#include "list_parser.h"
#include "energy_config.h"

#include <tuple>
#include <memory>
#include <vector>

#include "TF1.h"
#include "TKey.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TGraphErrors.h"
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

inline std::tuple<std::tuple<std::shared_ptr<TF1>, std::shared_ptr<TF1>>, std::tuple<std::shared_ptr<TGraphErrors>, std::shared_ptr<TGraphErrors>>> extractCorrectionFunctions(
    const char* file_name, 
    const char* tf1_shift_name="shift_fit_function",
    const char* tf1_sigma_name="sigma_ratio_fit_function",
    const char* gr_shift_name="gr_bin_shift",
    const char* gr_sigma_ratio_name="gr_bin_ratio") {

        TFile *infile = TFile::Open(file_name, "READ");
        if (infile->IsZombie()) {
            std::cout << "Error opening input TF1 file: [" << file_name << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TF1> tf1_shift = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_shift_name)));
        std::shared_ptr<TF1> tf1_sigma = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_sigma_name)));
        std::shared_ptr<TGraphErrors> gr_shift = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(gr_shift_name)));
        std::shared_ptr<TGraphErrors> gr_sigma = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(gr_sigma_ratio_name)));

        return std::make_tuple(std::make_tuple(tf1_shift, tf1_sigma), std::make_tuple(gr_shift, gr_sigma));
    }


void signal_efficiency(in_args input_args) {

    std::shared_ptr<energy_config> config = std::make_shared<energy_config>(input_args.energy_config_file);
    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.verbose);

    auto energy_binning = config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Extract external BDT cut values
    auto bdt_cuts = extract_bdt_cut_values(input_args.bdt_cut_values, energy_nbins);

    // Extract the efficiency correction functions
    auto corrections = extractCorrectionFunctions(input_args.eff_corr_function.c_str());

    auto shift_function = std::get<0>(std::get<0>(corrections));
    auto sigma_function = std::get<1>(std::get<0>(corrections));
    auto gr_shift = std::get<0>(std::get<1>(corrections));
    auto gr_sigma = std::get<1>(std::get<1>(corrections));

    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evt_parser->GetEvtTree());

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nTotal events: " << *(_data_fr.Count());
    std::cout << "\n\n***************************";

    auto get_bdt_cut = [&bdt_cuts, shift_function, sigma_function, gr_shift, gr_sigma] (const double energy_gev, const int energy_bin) -> double {
        double cut {bdt_cuts[energy_bin-1]};
        double cut_corr {0};
        if (energy_bin<=33)
            cut_corr = (cut - gr_shift->GetPointY(energy_bin-1))/gr_sigma->GetPointY(energy_bin-1); 
        else
            cut_corr = (cut - shift_function->Eval(energy_gev))/sigma_function->Eval(energy_gev);
        return cut_corr;
    };

    auto xtrl_loose_cut = [] (const double raw_energy_gev) -> double {
        return 3*pow(log10(raw_energy_gev/200), 2) + 12;
    };

    const double signal_spectral_index = -3;
    auto get_weight = [signal_spectral_index, &energy_binning] (const double energy_gev) -> double {
        return std::pow(energy_gev, signal_spectral_index+1)*std::pow(energy_binning[0], fabs(signal_spectral_index+1));
    };

    if (input_args.verbose) std::cout << "\n\nAnlysis running...\n\n";

    auto h_bdt_signal_not_passed = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("evt_w", get_weight, {"simu_energy_gev"})
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
                                    .Filter([] (const double tmva_value, const double tmva_cut) {return tmva_value < tmva_cut; }, {"tmva_classifier", "bdt_cut"})
                                    .Histo1D<double, double>({"h_bdt_signal_not_passed", "BDT event selection; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "evt_w");

    auto h_xtrl_tight_signal_not_passed = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("evt_w", get_weight, {"simu_energy_gev"})
                                    .Filter("xtrl>8.5")
                                    .Histo1D<double, double>({"h_xtrl_tight_signal_not_passed", "Event Selection - not XTRL tight; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "evt_w");

    auto h_xtrl_loose_signal_not_passed = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("simu_energy_gev", "simu_energy * 0.001")
                                    .Define("evt_w", get_weight, {"simu_energy_gev"})
                                    .Define("raw_energy_gev", "energy * 0.001")
                                    .Define("xtrl_loose_cut", xtrl_loose_cut, {"raw_energy_gev"})
                                    .Filter("xtrl>xtrl_loose_cut")
                                    .Histo1D<double, double>({"h_xtrl_loose_signal_not_passed", "Event Selection - not XTRL loose; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "evt_w"); 

    auto h_signal = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                            .Define("simu_energy_gev", "simu_energy * 0.001")
                            .Define("evt_w", get_weight, {"simu_energy_gev"})
                            .Histo1D<double>({"h_signal", "Signal Total Events; BGO Corr energy [GeV]; entries", energy_nbins, &energy_binning[0]}, "corr_energy_gev", "evt_w");

    TFile* outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen()) {
        std::cerr << "\n\nError writing output file... [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    h_bdt_signal_not_passed             ->Write();
    h_xtrl_tight_signal_not_passed      ->Write();
    h_xtrl_loose_signal_not_passed      ->Write();
    h_signal                            ->Write();
    
    outfile->Close();
}