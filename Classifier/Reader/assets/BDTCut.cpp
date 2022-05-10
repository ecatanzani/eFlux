#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TPDF.h"
#include "TPad.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TLegendEntry.h"
#include <ROOT/RDataFrame.hxx>

inline const std::string get_tree_name(const std::string stream) {
    const std::string file = stream.substr(0, stream.find('\n'));
    TFile* input_file = TFile::Open(file.c_str(), "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError reading input file [" << file << "]\n\n";
        exit(100);
    }
    std::string tree_name;
    for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
        auto key = dynamic_cast<TKey*>(keyAsObject);
        if (!strcmp(key->GetClassName(), "TTree")) {
            if (!strcmp(key->GetName(), "total_tree")) {
                tree_name = static_cast<std::string>(key->GetName());
                break;
            }
        }
    }
    input_file->Close();
    return tree_name;
}

std::shared_ptr<TChain> parse_input_list(const std::string input_list, const bool verbose) {

    auto parse_input_file = [](const std::string input_list) {
        std::ifstream input_file(input_list.c_str());
        if (!input_file.is_open())
        {
            std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
            exit(100);
        }
        std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
        input_file.close();
        return input_string;
    };

    std::istringstream input_stream(parse_input_file(input_list));
    std::shared_ptr<TChain> evtch;
    std::string tmp_str;
    bool first_elm {true};
    while (input_stream >> tmp_str) {
        if (first_elm) {
            evtch = std::make_shared<TChain> (get_tree_name(tmp_str).c_str(), "TMVA data set");
            first_elm = false;
        }
        evtch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
    return evtch;
}

struct energy_config {
    std::size_t n_bins;
    double min_event_energy {-999};
    double max_event_energy {-999};
    std::vector<float> energy_binning;
    
    void createLogBinning() {
        energy_binning = std::vector<float>(n_bins + 1, 0);
        double log_interval {(log10(max_event_energy)-log10(min_event_energy))/n_bins};
        for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
            energy_binning[bIdx] = pow(10, log10(min_event_energy) + bIdx * log_interval);
    }
};

inline std::string parse_config_file(const char* config_file) {
	std::ifstream input_file(config_file);
	if (!input_file.is_open()) {
		std::cerr << "\nInput config file not found [" << config_file << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

inline energy_config get_config_info(const std::string parsed_config) {
	energy_config config_pars;
    std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;
	while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> config_pars.n_bins;
		if (!strcmp(tmp_str.c_str(), "min_event_energy")) {
			input_stream >> tmp_str;
			config_pars.min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy")) {
			input_stream >> tmp_str;
			config_pars.max_event_energy = stod(tmp_str, &sz);
		}
	}
    config_pars.createLogBinning();
    return config_pars;
}

inline std::vector<float> parse_energy_config(const char* config_file) {
    return get_config_info(parse_config_file(config_file)).energy_binning;
}

void BDTCut(
    std::string input_list, 
    const char* output_file,
    const unsigned int threads = 1,
    const bool verbose = true) {

    auto chain = parse_input_list(input_list, verbose);

    // Initialize RDF
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame fr(*chain);

    auto h_tmva_classifier_20_100 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                        .Histo1D({"h_tmva_classifier_20_100", "TMVA classifier 20 GeV - 100 GeV", 1000, -1, 1}, "tmva_classifier");

    auto h_tmva_classifier_100_250 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                        .Histo1D({"h_tmva_classifier_100_250", "TMVA classifier 100 GeV - 250 GeV", 1000, -1, 1}, "tmva_classifier");

    auto h_tmva_classifier_250_500 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                        .Histo1D({"h_tmva_classifier_250_500", "TMVA classifier 250 GeV - 500 GeV", 1000, -1, 1}, "tmva_classifier");

    auto h_tmva_classifier_500_1000 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                        .Histo1D({"h_tmva_classifier_500_1000", "TMVA classifier 500 GeV - 1 TeV", 1000, -1, 1}, "tmva_classifier");

    auto h_tmva_classifier_1000_3000 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                        .Histo1D({"h_tmva_classifier_1000_3000", "TMVA classifier 1 TeV - 3 TeV", 1000, -1, 1}, "tmva_classifier");

    auto h_tmva_classifier_3000 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Filter("corr_energy_gev>=3000")
                                        .Histo1D({"h_tmva_classifier_3000", "TMVA classifier 3 TeV - 10 TeV", 1000, -1, 1}, "tmva_classifier");

    TFile *output = TFile::Open(output_file, "RECREATE");
    if (output->IsZombie()) {
        std::cerr << "Error writing output ROOT file [" << output_file << "]\n\n";
        exit(100);
    }

    h_tmva_classifier_20_100->Write();
    h_tmva_classifier_100_250->Write();
    h_tmva_classifier_250_500->Write();
    h_tmva_classifier_500_1000->Write();
    h_tmva_classifier_1000_3000->Write();
    h_tmva_classifier_3000->Write();

    output->Close();
}

void BDTMCCut(
    std::string input_list, 
    const char* output_file,
    const char* energy_config_file,
    const unsigned int threads = 1,
    const bool is_electron = true,
    const bool verbose = true) {

    // Parse input file list
    auto chain = parse_input_list(input_list, verbose);

    // Extract energy binning from config file
    auto energy_binning = parse_energy_config(energy_config_file);
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Initialize RDF
    ROOT::EnableImplicitMT(threads);
    ROOT::RDataFrame fr(*chain);

    const double signal_spectral_index {is_electron ? -3.0 : -2.7};
    auto get_weight = [signal_spectral_index, &energy_binning] (const double energy_gev) -> double {
        return std::pow(energy_gev, signal_spectral_index+1)*std::pow(energy_binning[0], fabs(signal_spectral_index+1));
    };

    auto h_tmva_classifier_20_100 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("evt_w", get_weight, {"corr_energy_gev"})
                                        .Filter("corr_energy_gev>=20 && corr_energy_gev<100")
                                        .Histo1D({"h_tmva_classifier_20_100", "TMVA classifier 20 GeV - 100 GeV", 1000, -1, 1}, "tmva_classifier", "evt_w");
                                        
    auto h_tmva_classifier_100_250 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("evt_w", get_weight, {"corr_energy_gev"})
                                        .Filter("corr_energy_gev>=100 && corr_energy_gev<250")
                                        .Histo1D({"h_tmva_classifier_100_250", "TMVA classifier 100 GeV - 250 GeV", 1000, -1, 1}, "tmva_classifier", "evt_w");
                                        
    auto h_tmva_classifier_250_500 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("evt_w", get_weight, {"corr_energy_gev"})
                                        .Filter("corr_energy_gev>=250 && corr_energy_gev<500")
                                        .Histo1D({"h_tmva_classifier_250_500", "TMVA classifier 250 GeV - 500 GeV", 1000, -1, 1}, "tmva_classifier", "evt_w");
                                        
    auto h_tmva_classifier_500_1000 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("evt_w", get_weight, {"corr_energy_gev"})
                                        .Filter("corr_energy_gev>=500 && corr_energy_gev<1000")
                                        .Histo1D({"h_tmva_classifier_500_1000", "TMVA classifier 500 GeV - 1 TeV", 1000, -1, 1}, "tmva_classifier", "evt_w");
                                        
    auto h_tmva_classifier_1000_3000 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("evt_w", get_weight, {"corr_energy_gev"})
                                        .Filter("corr_energy_gev>=1000 && corr_energy_gev<3000")
                                        .Histo1D({"h_tmva_classifier_1000_3000", "TMVA classifier 1 TeV - 3 TeV", 1000, -1, 1}, "tmva_classifier", "evt_w");
                                        
    auto h_tmva_classifier_3000 = fr.Define("corr_energy_gev", "energy_corr * 0.001")
                                        .Define("evt_w", get_weight, {"corr_energy_gev"})
                                        .Filter("corr_energy_gev>=3000")
                                        .Histo1D({"h_tmva_classifier_3000", "TMVA classifier 3 TeV - 10 TeV", 1000, -1, 1}, "tmva_classifier", "evt_w");
                                        
    TFile *output = TFile::Open(output_file, "RECREATE");
    if (output->IsZombie()) {
        std::cerr << "Error writing output ROOT file [" << output_file << "]\n\n";
        exit(100);
    }

    h_tmva_classifier_20_100->Write();
    h_tmva_classifier_100_250->Write();
    h_tmva_classifier_250_500->Write();
    h_tmva_classifier_500_1000->Write();
    h_tmva_classifier_1000_3000->Write();
    h_tmva_classifier_3000->Write();

    output->Close();
}

void MergeMCDataBDT(
    const char* input_data_file,
    const char* input_mc_electron_file,
    const char* input_mc_proton_file,
    const char* output_file,
    const bool logy = true,
    const bool verbose = true) {

        // Extract data histos
        auto data_file = TFile::Open(input_data_file, "READ");
        if (data_file->IsZombie()) {
            std::cerr << "\n\nError readin input data file [" << input_data_file << "]\n\n";
            exit(100);
        }

        auto data_h_tmva_classifier_20_100      = static_cast<TH1D*>(data_file->Get("h_tmva_classifier_20_100"));
        auto data_h_tmva_classifier_100_250     = static_cast<TH1D*>(data_file->Get("h_tmva_classifier_100_250"));
        auto data_h_tmva_classifier_250_500     = static_cast<TH1D*>(data_file->Get("h_tmva_classifier_250_500"));
        auto data_h_tmva_classifier_500_1000    = static_cast<TH1D*>(data_file->Get("h_tmva_classifier_500_1000"));
        auto data_h_tmva_classifier_1000_3000   = static_cast<TH1D*>(data_file->Get("h_tmva_classifier_1000_3000"));
        auto data_h_tmva_classifier_3000        = static_cast<TH1D*>(data_file->Get("h_tmva_classifier_3000"));

        data_h_tmva_classifier_20_100           ->SetLineWidth(2);
        data_h_tmva_classifier_100_250          ->SetLineWidth(2);
        data_h_tmva_classifier_250_500          ->SetLineWidth(2);
        data_h_tmva_classifier_500_1000         ->SetLineWidth(2);
        data_h_tmva_classifier_1000_3000        ->SetLineWidth(2);
        data_h_tmva_classifier_3000             ->SetLineWidth(2);

        data_h_tmva_classifier_20_100           ->SetLineColor(kBlack);
        data_h_tmva_classifier_100_250          ->SetLineColor(kBlack);
        data_h_tmva_classifier_250_500          ->SetLineColor(kBlack);
        data_h_tmva_classifier_500_1000         ->SetLineColor(kBlack);
        data_h_tmva_classifier_1000_3000        ->SetLineColor(kBlack);
        data_h_tmva_classifier_3000             ->SetLineColor(kBlack);

        data_h_tmva_classifier_20_100           ->SetDirectory(0);
        data_h_tmva_classifier_100_250          ->SetDirectory(0);
        data_h_tmva_classifier_250_500          ->SetDirectory(0);
        data_h_tmva_classifier_500_1000         ->SetDirectory(0);
        data_h_tmva_classifier_1000_3000        ->SetDirectory(0);
        data_h_tmva_classifier_3000             ->SetDirectory(0);

        data_h_tmva_classifier_20_100           ->SetTitle("DATA");
        data_h_tmva_classifier_100_250          ->SetTitle("DATA");
        data_h_tmva_classifier_250_500          ->SetTitle("DATA");
        data_h_tmva_classifier_500_1000         ->SetTitle("DATA");
        data_h_tmva_classifier_1000_3000        ->SetTitle("DATA");
        data_h_tmva_classifier_3000             ->SetTitle("DATA");

        data_file->Close();

        // Extract MC electron histos
        auto mcelectron_file = TFile::Open(input_mc_electron_file, "READ");
        if (mcelectron_file->IsZombie()) {
            std::cerr << "\n\nError readin input MC electron file [" << input_mc_electron_file << "]\n\n";
            exit(100);
        }

        auto mc_electron_h_tmva_classifier_20_100       = static_cast<TH1D*>(mcelectron_file->Get("h_tmva_classifier_20_100"));
        auto mc_electron_h_tmva_classifier_100_250      = static_cast<TH1D*>(mcelectron_file->Get("h_tmva_classifier_100_250"));
        auto mc_electron_h_tmva_classifier_250_500      = static_cast<TH1D*>(mcelectron_file->Get("h_tmva_classifier_250_500"));
        auto mc_electron_h_tmva_classifier_500_1000     = static_cast<TH1D*>(mcelectron_file->Get("h_tmva_classifier_500_1000"));
        auto mc_electron_h_tmva_classifier_1000_3000    = static_cast<TH1D*>(mcelectron_file->Get("h_tmva_classifier_1000_3000"));
        auto mc_electron_h_tmva_classifier_3000         = static_cast<TH1D*>(mcelectron_file->Get("h_tmva_classifier_3000"));

        mc_electron_h_tmva_classifier_20_100            ->SetLineWidth(2);
        mc_electron_h_tmva_classifier_100_250           ->SetLineWidth(2);
        mc_electron_h_tmva_classifier_250_500           ->SetLineWidth(2);
        mc_electron_h_tmva_classifier_500_1000          ->SetLineWidth(2);
        mc_electron_h_tmva_classifier_1000_3000         ->SetLineWidth(2);
        mc_electron_h_tmva_classifier_3000              ->SetLineWidth(2);

        mc_electron_h_tmva_classifier_20_100            ->SetLineColor(kBlue+2);
        mc_electron_h_tmva_classifier_100_250           ->SetLineColor(kBlue+2);
        mc_electron_h_tmva_classifier_250_500           ->SetLineColor(kBlue+2);
        mc_electron_h_tmva_classifier_500_1000          ->SetLineColor(kBlue+2);
        mc_electron_h_tmva_classifier_1000_3000         ->SetLineColor(kBlue+2);
        mc_electron_h_tmva_classifier_3000              ->SetLineColor(kBlue+2);

        mc_electron_h_tmva_classifier_20_100            ->SetMarkerColor(kBlue+2);
        mc_electron_h_tmva_classifier_100_250           ->SetMarkerColor(kBlue+2);
        mc_electron_h_tmva_classifier_250_500           ->SetMarkerColor(kBlue+2);
        mc_electron_h_tmva_classifier_500_1000          ->SetMarkerColor(kBlue+2);
        mc_electron_h_tmva_classifier_1000_3000         ->SetMarkerColor(kBlue+2);
        mc_electron_h_tmva_classifier_3000              ->SetMarkerColor(kBlue+2);

        mc_electron_h_tmva_classifier_20_100            ->SetDirectory(0);
        mc_electron_h_tmva_classifier_100_250           ->SetDirectory(0);
        mc_electron_h_tmva_classifier_250_500           ->SetDirectory(0);
        mc_electron_h_tmva_classifier_500_1000          ->SetDirectory(0);
        mc_electron_h_tmva_classifier_1000_3000         ->SetDirectory(0);
        mc_electron_h_tmva_classifier_3000              ->SetDirectory(0);

        mc_electron_h_tmva_classifier_20_100            ->SetTitle("MC electron");
        mc_electron_h_tmva_classifier_100_250           ->SetTitle("MC electron");
        mc_electron_h_tmva_classifier_250_500           ->SetTitle("MC electron");
        mc_electron_h_tmva_classifier_500_1000          ->SetTitle("MC electron");
        mc_electron_h_tmva_classifier_1000_3000         ->SetTitle("MC electron");
        mc_electron_h_tmva_classifier_3000              ->SetTitle("MC electron");

        mcelectron_file->Close();

        // Extract MC proton histos
        auto mcproton_file = TFile::Open(input_mc_proton_file, "READ");
        if (mcproton_file->IsZombie()) {
            std::cerr << "\n\nError readin input MC proton file [" << input_mc_proton_file << "]\n\n";
            exit(100);
        }

        auto mc_proton_h_tmva_classifier_20_100         = static_cast<TH1D*>(mcproton_file->Get("h_tmva_classifier_20_100"));
        auto mc_proton_h_tmva_classifier_100_250        = static_cast<TH1D*>(mcproton_file->Get("h_tmva_classifier_100_250"));
        auto mc_proton_h_tmva_classifier_250_500        = static_cast<TH1D*>(mcproton_file->Get("h_tmva_classifier_250_500"));
        auto mc_proton_h_tmva_classifier_500_1000       = static_cast<TH1D*>(mcproton_file->Get("h_tmva_classifier_500_1000"));
        auto mc_proton_h_tmva_classifier_1000_3000      = static_cast<TH1D*>(mcproton_file->Get("h_tmva_classifier_1000_3000"));
        auto mc_proton_h_tmva_classifier_3000           = static_cast<TH1D*>(mcproton_file->Get("h_tmva_classifier_3000"));

        mc_proton_h_tmva_classifier_20_100              ->SetLineWidth(2);
        mc_proton_h_tmva_classifier_100_250             ->SetLineWidth(2);
        mc_proton_h_tmva_classifier_250_500             ->SetLineWidth(2);
        mc_proton_h_tmva_classifier_500_1000            ->SetLineWidth(2);
        mc_proton_h_tmva_classifier_1000_3000           ->SetLineWidth(2);
        mc_proton_h_tmva_classifier_3000                ->SetLineWidth(2);

        mc_proton_h_tmva_classifier_20_100              ->SetLineColor(kRed+2);
        mc_proton_h_tmva_classifier_100_250             ->SetLineColor(kRed+2);
        mc_proton_h_tmva_classifier_250_500             ->SetLineColor(kRed+2);
        mc_proton_h_tmva_classifier_500_1000            ->SetLineColor(kRed+2);
        mc_proton_h_tmva_classifier_1000_3000           ->SetLineColor(kRed+2);
        mc_proton_h_tmva_classifier_3000                ->SetLineColor(kRed+2);

        mc_proton_h_tmva_classifier_20_100              ->SetMarkerColor(kRed+2);
        mc_proton_h_tmva_classifier_100_250             ->SetMarkerColor(kRed+2);
        mc_proton_h_tmva_classifier_250_500             ->SetMarkerColor(kRed+2);
        mc_proton_h_tmva_classifier_500_1000            ->SetMarkerColor(kRed+2);
        mc_proton_h_tmva_classifier_1000_3000           ->SetMarkerColor(kRed+2);
        mc_proton_h_tmva_classifier_3000                ->SetMarkerColor(kRed+2);

        mc_proton_h_tmva_classifier_20_100              ->SetDirectory(0);
        mc_proton_h_tmva_classifier_100_250             ->SetDirectory(0);
        mc_proton_h_tmva_classifier_250_500             ->SetDirectory(0);
        mc_proton_h_tmva_classifier_500_1000            ->SetDirectory(0);
        mc_proton_h_tmva_classifier_1000_3000           ->SetDirectory(0);
        mc_proton_h_tmva_classifier_3000                ->SetDirectory(0);

        mc_proton_h_tmva_classifier_20_100              ->SetTitle("MC proton");
        mc_proton_h_tmva_classifier_100_250             ->SetTitle("MC proton");
        mc_proton_h_tmva_classifier_250_500             ->SetTitle("MC proton");
        mc_proton_h_tmva_classifier_500_1000            ->SetTitle("MC proton");
        mc_proton_h_tmva_classifier_1000_3000           ->SetTitle("MC proton");
        mc_proton_h_tmva_classifier_3000                ->SetTitle("MC proton");

        mcproton_file->Close();

        // Extract scale values for electron and proton mc histos
        double h_data_proton_max_20_100         {1};
        double h_data_proton_max_100_250        {1};
        double h_data_proton_max_250_500        {1};
        double h_data_proton_max_500_1000       {1};
        double h_data_proton_max_1000_3000      {1};
        double h_data_proton_max_3000           {1};

        double h_data_electron_max_20_100       {1};
        double h_data_electron_max_100_250      {1};
        double h_data_electron_max_250_500      {1};
        double h_data_electron_max_500_1000     {1};
        double h_data_electron_max_1000_3000    {1};
        double h_data_electron_max_3000         {1};

        data_h_tmva_classifier_20_100       ->GetXaxis()->SetRangeUser(-1, 0);
        data_h_tmva_classifier_100_250      ->GetXaxis()->SetRangeUser(-1, 0);
        data_h_tmva_classifier_250_500      ->GetXaxis()->SetRangeUser(-1, 0);
        data_h_tmva_classifier_500_1000     ->GetXaxis()->SetRangeUser(-1, 0);
        data_h_tmva_classifier_1000_3000    ->GetXaxis()->SetRangeUser(-1, 0);
        data_h_tmva_classifier_3000         ->GetXaxis()->SetRangeUser(-1, 0);

        h_data_proton_max_20_100            = data_h_tmva_classifier_20_100->GetMaximum();
        h_data_proton_max_100_250           = data_h_tmva_classifier_100_250->GetMaximum();
        h_data_proton_max_250_500           = data_h_tmva_classifier_250_500->GetMaximum();
        h_data_proton_max_500_1000          = data_h_tmva_classifier_500_1000->GetMaximum();
        h_data_proton_max_1000_3000         = data_h_tmva_classifier_1000_3000->GetMaximum();
        h_data_proton_max_3000              = data_h_tmva_classifier_3000->GetMaximum();

        data_h_tmva_classifier_20_100       ->GetXaxis()->SetRangeUser(0, 1);
        data_h_tmva_classifier_100_250      ->GetXaxis()->SetRangeUser(0, 1);
        data_h_tmva_classifier_250_500      ->GetXaxis()->SetRangeUser(0, 1);
        data_h_tmva_classifier_500_1000     ->GetXaxis()->SetRangeUser(0, 1);
        data_h_tmva_classifier_1000_3000    ->GetXaxis()->SetRangeUser(0, 1);
        data_h_tmva_classifier_3000         ->GetXaxis()->SetRangeUser(0, 1);

        h_data_electron_max_20_100            = data_h_tmva_classifier_20_100->GetMaximum();
        h_data_electron_max_100_250           = data_h_tmva_classifier_100_250->GetMaximum();
        h_data_electron_max_250_500           = data_h_tmva_classifier_250_500->GetMaximum();
        h_data_electron_max_500_1000          = data_h_tmva_classifier_500_1000->GetMaximum();
        h_data_electron_max_1000_3000         = data_h_tmva_classifier_1000_3000->GetMaximum();
        h_data_electron_max_3000              = data_h_tmva_classifier_3000->GetMaximum();

        // Rebuild original data range
        data_h_tmva_classifier_20_100       ->GetXaxis()->SetRangeUser(-1, 1);
        data_h_tmva_classifier_100_250      ->GetXaxis()->SetRangeUser(-1, 1);
        data_h_tmva_classifier_250_500      ->GetXaxis()->SetRangeUser(-1, 1);
        data_h_tmva_classifier_500_1000     ->GetXaxis()->SetRangeUser(-1, 1);
        data_h_tmva_classifier_1000_3000    ->GetXaxis()->SetRangeUser(-1, 1);
        data_h_tmva_classifier_3000         ->GetXaxis()->SetRangeUser(-1, 1);

        // Scale MC histos
        auto scale_mc_histo = [](TH1D* h_mc, double data_max_val)
        {
            if (data_max_val) {
                if (h_mc->GetMaximum() >= data_max_val)
                    h_mc->Scale(h_mc->GetMaximum()/data_max_val);
                else
                    h_mc->Scale(data_max_val/h_mc->GetMaximum());
            }
        };

        scale_mc_histo(mc_electron_h_tmva_classifier_20_100, h_data_electron_max_20_100);
        scale_mc_histo(mc_electron_h_tmva_classifier_100_250, h_data_electron_max_100_250);
        scale_mc_histo(mc_electron_h_tmva_classifier_250_500, h_data_electron_max_250_500);
        scale_mc_histo(mc_electron_h_tmva_classifier_500_1000, h_data_electron_max_500_1000);
        scale_mc_histo(mc_electron_h_tmva_classifier_1000_3000, h_data_electron_max_1000_3000);
        scale_mc_histo(mc_electron_h_tmva_classifier_3000, h_data_electron_max_3000);

        scale_mc_histo(mc_proton_h_tmva_classifier_20_100, h_data_proton_max_20_100);
        scale_mc_histo(mc_proton_h_tmva_classifier_100_250, h_data_proton_max_100_250);
        scale_mc_histo(mc_proton_h_tmva_classifier_250_500, h_data_proton_max_250_500);
        scale_mc_histo(mc_proton_h_tmva_classifier_500_1000, h_data_proton_max_500_1000);
        scale_mc_histo(mc_proton_h_tmva_classifier_1000_3000, h_data_proton_max_1000_3000);
        scale_mc_histo(mc_proton_h_tmva_classifier_3000, h_data_proton_max_3000);

        // Build canvas
        TCanvas print_canvas("print_canvas", "print_canvas");
        print_canvas.SetTicks();

        // 20 GeV - 100 GeV plot
        data_h_tmva_classifier_20_100->Draw();
        mc_electron_h_tmva_classifier_20_100->Draw("same, histo");
        mc_proton_h_tmva_classifier_20_100->Draw("same, histo");

        if (logy)
            gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        auto legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        TPaveLabel label(0.0, 0.95, 0.3, 1, "BDT classifier - 20 GeV - 100 GeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_classifier_data_mc.pdf(","Title:BDT classifier - 20 GeV - 100 GeV");

        // 100 GeV - 250 GeV plot
        data_h_tmva_classifier_100_250->Draw();
        mc_electron_h_tmva_classifier_100_250->Draw("same, histo");
        mc_proton_h_tmva_classifier_100_250->Draw("same, histo");

        if (logy)
            gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT classifier - 100 GeV - 250 GeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_classifier_data_mc.pdf","Title:BDT classifier - 100 GeV - 250 GeV");

        // 250 GeV - 500 GeV plot
        data_h_tmva_classifier_250_500->Draw();
        mc_electron_h_tmva_classifier_250_500->Draw("same, histo");
        mc_proton_h_tmva_classifier_250_500->Draw("same, histo");

        if (logy)
            gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT classifier - 250 GeV - 500 GeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_classifier_data_mc.pdf","Title:BDT classifier - 250 GeV - 500 GeV");

        // 500 GeV - 1 TeV plot
        data_h_tmva_classifier_500_1000->Draw();
        mc_electron_h_tmva_classifier_500_1000->Draw("same, histo");
        mc_proton_h_tmva_classifier_500_1000->Draw("same, histo");

        if (logy)
            gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT classifier - 500 GeV - 1 TeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_classifier_data_mc.pdf","Title:BDT classifier - 500 GeV - 1 TeV");

        // 1 TeV - 3 TeV plot
        data_h_tmva_classifier_1000_3000->Draw();
        mc_electron_h_tmva_classifier_1000_3000->Draw("same, histo");
        mc_proton_h_tmva_classifier_1000_3000->Draw("same, histo");

        if (logy)
            gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT classifier - 1 TeV - 3 TeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_classifier_data_mc.pdf","Title:BDT classifier - 1 TeV - 3 TeV");

        // 3 TeV plot
        data_h_tmva_classifier_3000->Draw();
        mc_electron_h_tmva_classifier_3000->Draw("same, histo");
        mc_proton_h_tmva_classifier_3000->Draw("same, histo");

        if (logy)
            gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT classifier - 3 TeV - 10 TeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_classifier_data_mc.pdf)","Title:BDT classifier - 3 TeV - 10 TeV");

        // Save plots on ROOT canvas
        TFile out_file(output_file, "RECREATE");
        TCanvas out_canvas("bdt_classifier_data_mc_comparison", "bdt_classifier_data_mc_comparison");

        out_canvas.Divide(3, 2);

        out_canvas.cd(1);

        data_h_tmva_classifier_20_100->Draw();
        mc_electron_h_tmva_classifier_20_100->Draw("same, histo");
        mc_proton_h_tmva_classifier_20_100->Draw("same, histo");

        gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        out_canvas.cd(2);

        data_h_tmva_classifier_100_250->Draw();
        mc_electron_h_tmva_classifier_100_250->Draw("same, histo");
        mc_proton_h_tmva_classifier_100_250->Draw("same, histo");

        gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        out_canvas.cd(3);

        data_h_tmva_classifier_250_500->Draw();
        mc_electron_h_tmva_classifier_250_500->Draw("same, histo");
        mc_proton_h_tmva_classifier_250_500->Draw("same, histo");

        gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        out_canvas.cd(4);

        data_h_tmva_classifier_500_1000->Draw();
        mc_electron_h_tmva_classifier_500_1000->Draw("same, histo");
        mc_proton_h_tmva_classifier_500_1000->Draw("same, histo");

        gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        out_canvas.cd(5);

        data_h_tmva_classifier_1000_3000->Draw();
        mc_electron_h_tmva_classifier_1000_3000->Draw("same, histo");
        mc_proton_h_tmva_classifier_1000_3000->Draw("same, histo");

        gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        out_canvas.cd(6);

        data_h_tmva_classifier_3000->Draw();
        mc_electron_h_tmva_classifier_3000->Draw("same, histo");
        mc_proton_h_tmva_classifier_3000->Draw("same, histo");

        gPad->SetLogy();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);

        legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        out_canvas.Write();
        out_file.Close();
    }