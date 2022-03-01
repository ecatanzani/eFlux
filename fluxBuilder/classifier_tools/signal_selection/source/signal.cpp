#include "signal.h"
#include "bdt_config.h"
#include "list_parser.h"
#include "energy_config.h"

#include <memory>

#include "TFile.h"
#include <ROOT/RDataFrame.hxx>

inline const char* get_bdt_config_file(const char* energy_config_file) {
    return (std::string(energy_config_file).substr(0, std::string(energy_config_file).find("/eFlux")+6) + std::string("/Classifier/Reader/config/classifier.conf")).c_str();
}

void signal_selection(in_args input_args) {

    std::shared_ptr<energy_config> config = std::make_shared<energy_config>(input_args.energy_config_file);
    std::shared_ptr<parser> evt_parser = std::make_unique<parser>(input_args.input_list, input_args.verbose);
    std::shared_ptr<bdt_config> cl_config = std::make_unique<bdt_config>(get_bdt_config_file(input_args.energy_config_file));

    auto energy_binning = config->GetEnergyBinning();
    auto energy_nbins = (int)energy_binning.size() - 1;

    ROOT::EnableImplicitMT(input_args.threads);
    ROOT::RDataFrame _data_fr(*evt_parser->GetEvtTree());

    std::cout << "\n\n**** Filter statistics ****\n";
    std::cout << "***************************\n";
    std::cout << "\nTotal events: " << *(_data_fr.Count());
    std::cout << "\n\n***************************";

    if (input_args.verbose) std::cout << "\n\nAnlysis running...\n\n";

    auto get_bdt_cut = [cl_config] (const double energy_gev, const int energy_bin) -> double {
        double cut {0};

        if (energy_gev>=10 && energy_gev<100)
            cut = cl_config->GetLowEnergyBDTCut();
        else if (energy_gev>=100 && energy_gev<1000)
            cut = energy_bin != 32 ? cl_config->GetMidEnergyBDTCut() : cl_config->GetMidEnergyBDTCut() - 0.07;
        else if (energy_gev>=1000 && energy_gev<=10000)     
            cut = cl_config->GetHighEnergyBDTCut();
        
        return cut;
    };

    auto xtrl_loose_cut = [] (const double raw_energy_gev) -> double {
        return 3*pow(log10(raw_energy_gev/200), 2) + 12;
    };

    auto h_bdt_selection = _data_fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                    .Define("bdt_cut", get_bdt_cut, {"corr_energy_gev", "energy_bin"})
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