#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TEfficiency.h"
#include <ROOT/RDataFrame.hxx>

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

std::vector<float> createLogBinning(const double eMin, const double eMax, const std::size_t n_bins) {
	std::vector<float> binning(n_bins + 1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);
	return binning;
}

std::vector<float> createLinearBinning(float a, float b, std::size_t N) {
	float h = (b - a) / static_cast<float>(N);
	std::vector<float> xs(N + 1);
	std::vector<float>::iterator x;
	float val;
	for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
		*x = val;
	return xs;
}

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
        if (!strcmp(key->GetClassName(), "TTree"))
            tree_name = static_cast<std::string>(key->GetName());
    }
    input_file->Close();
    return tree_name;
}

void buildEfficiency(
    const char* input_file,
    const char* output_file,
    const char* energy_config_file,
    const bool verbose,
    const unsigned int threads=1) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        // Chain input files
        TFile *infile = TFile::Open(input_file, "READ");
        if (!infile->IsOpen()) {
            std::cerr << "\n\nError reading input file [" << input_file << "]\n\n";
            exit(100);
        }

        // Trigger histos
        auto h_trigger_efficiency_accepted_het_tight_xtrl = static_cast<TH1D*>(infile->Get("h_trigger_efficiency_accepted_het_tight_xtrl"));
        auto h_trigger_efficiency_accepted_het_let_tight_xtrl = static_cast<TH1D*>(infile->Get("h_trigger_efficiency_accepted_het_let_tight_xtrl"));
        auto h_trigger_efficiency_accepted_het_bdt = static_cast<TH1D*>(infile->Get("h_trigger_efficiency_accepted_het_bdt"));
        auto h_trigger_efficiency_accepted_het_let_bdt = static_cast<TH1D*>(infile->Get("h_trigger_efficiency_accepted_het_let_bdt"));
         
        auto h_trigger_efficiency_accepted_het_loose_xtrl = static_cast<TH1D*>(infile->Get("h_trigger_efficiency_accepted_het_loose_xtrl"));
        auto h_trigger_efficiency_accepted_het_let_loose_xtrl = static_cast<TH1D*>(infile->Get("h_trigger_efficiency_accepted_het_let_loose_xtrl"));

        h_trigger_efficiency_accepted_het_tight_xtrl->SetDirectory(0);
        h_trigger_efficiency_accepted_het_let_tight_xtrl->SetDirectory(0);
        h_trigger_efficiency_accepted_het_bdt->SetDirectory(0);
        h_trigger_efficiency_accepted_het_let_bdt->SetDirectory(0);

        h_trigger_efficiency_accepted_het_loose_xtrl->SetDirectory(0);
        h_trigger_efficiency_accepted_het_let_loose_xtrl->SetDirectory(0);

        // MaxRMS histos
        auto h_maxrms_efficiency_accepted_tight_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_efficiency_accepted_tight_xtrl"));
        auto h_maxrms_efficiency_total_tight_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_efficiency_total_tight_xtrl"));
        auto h_maxrms_efficiency_accepted_bdt = static_cast<TH1D*>(infile->Get("h_maxrms_efficiency_accepted_bdt"));
        auto h_maxrms_efficiency_total_bdt = static_cast<TH1D*>(infile->Get("h_maxrms_efficiency_total_bdt"));

        auto h_maxrms_efficiency_accepted_loose_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_efficiency_accepted_loose_xtrl"));
        auto h_maxrms_efficiency_total_loose_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_efficiency_total_loose_xtrl"));

        h_maxrms_efficiency_accepted_tight_xtrl->SetDirectory(0);
        h_maxrms_efficiency_total_tight_xtrl->SetDirectory(0);
        h_maxrms_efficiency_accepted_bdt->SetDirectory(0);
        h_maxrms_efficiency_total_bdt->SetDirectory(0);
        
        h_maxrms_efficiency_accepted_loose_xtrl->SetDirectory(0);
        h_maxrms_efficiency_total_loose_xtrl->SetDirectory(0);

        // nbarlayer13 histos
        auto h_nbarlayer13_efficiency_accepted_tight_xtrl = static_cast<TH1D*>(infile->Get("h_nbarlayer13_efficiency_accepted_tight_xtrl"));
        auto h_nbarlayer13_efficiency_total_tight_xtrl = static_cast<TH1D*>(infile->Get("h_nbarlayer13_efficiency_total_tight_xtrl"));
        auto h_nbarlayer13_efficiency_accepted_bdt = static_cast<TH1D*>(infile->Get("h_nbarlayer13_efficiency_accepted_bdt"));
        auto h_nbarlayer13_efficiency_total_bdt = static_cast<TH1D*>(infile->Get("h_nbarlayer13_efficiency_total_bdt"));

        auto h_nbarlayer13_efficiency_accepted_loose_xtrl = static_cast<TH1D*>(infile->Get("h_nbarlayer13_efficiency_accepted_loose_xtrl"));
        auto h_nbarlayer13_efficiency_total_loose_xtrl = static_cast<TH1D*>(infile->Get("h_nbarlayer13_efficiency_total_loose_xtrl"));

        h_nbarlayer13_efficiency_accepted_tight_xtrl->SetDirectory(0);
        h_nbarlayer13_efficiency_total_tight_xtrl->SetDirectory(0);
        h_nbarlayer13_efficiency_accepted_bdt->SetDirectory(0);
        h_nbarlayer13_efficiency_total_bdt->SetDirectory(0);
        
        h_nbarlayer13_efficiency_accepted_loose_xtrl->SetDirectory(0);
        h_nbarlayer13_efficiency_total_loose_xtrl->SetDirectory(0);

        // MaxRMS and nbarlayer13 histos
        auto h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl"));
        auto h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl"));
        auto h_maxrms_and_nbarlayer13_efficiency_accepted_bdt = static_cast<TH1D*>(infile->Get("h_maxrms_and_nbarlayer13_efficiency_accepted_bdt"));
        auto h_maxrms_and_nbarlayer13_efficiency_total_bdt = static_cast<TH1D*>(infile->Get("h_maxrms_and_nbarlayer13_efficiency_total_bdt"));

        auto h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl"));
        auto h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl = static_cast<TH1D*>(infile->Get("h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl"));

        h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl->SetDirectory(0);
        h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl->SetDirectory(0);
        h_maxrms_and_nbarlayer13_efficiency_accepted_bdt->SetDirectory(0);
        h_maxrms_and_nbarlayer13_efficiency_total_bdt->SetDirectory(0);
        
        h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl->SetDirectory(0);
        h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl->SetDirectory(0);

        // Track Selection histos
        auto h_track_efficiency_accepted_tight_xtrl = static_cast<TH1D*>(infile->Get("h_track_efficiency_accepted_tight_xtrl"));
        auto h_track_efficiency_total_tight_xtrl = static_cast<TH1D*>(infile->Get("h_track_efficiency_total_tight_xtrl"));
        auto h_track_efficiency_accepted_bdt = static_cast<TH1D*>(infile->Get("h_track_efficiency_accepted_bdt"));
        auto h_track_efficiency_total_bdt = static_cast<TH1D*>(infile->Get("h_track_efficiency_total_bdt"));

        auto h_track_efficiency_accepted_loose_xtrl = static_cast<TH1D*>(infile->Get("h_track_efficiency_accepted_loose_xtrl"));
        auto h_track_efficiency_total_loose_xtrl = static_cast<TH1D*>(infile->Get("h_track_efficiency_total_loose_xtrl"));

        h_track_efficiency_accepted_tight_xtrl->SetDirectory(0);
        h_track_efficiency_total_tight_xtrl->SetDirectory(0);
        h_track_efficiency_accepted_bdt->SetDirectory(0);
        h_track_efficiency_total_bdt->SetDirectory(0);
        
        h_track_efficiency_accepted_loose_xtrl->SetDirectory(0);
        h_track_efficiency_total_loose_xtrl->SetDirectory(0);

        // PSD-STK match histos
        auto h_psdstkmatch_efficiency_accepted_tight_xtrl = static_cast<TH1D*>(infile->Get("h_psdstkmatch_efficiency_accepted_tight_xtrl"));
        auto h_psdstkmatch_efficiency_total_tight_xtrl = static_cast<TH1D*>(infile->Get("h_psdstkmatch_efficiency_total_tight_xtrl"));
        auto h_psdstkmatch_efficiency_accepted_bdt = static_cast<TH1D*>(infile->Get("h_psdstkmatch_efficiency_accepted_bdt"));
        auto h_psdstkmatch_efficiency_total_bdt = static_cast<TH1D*>(infile->Get("h_psdstkmatch_efficiency_total_bdt"));

        auto h_psdstkmatch_efficiency_accepted_loose_xtrl = static_cast<TH1D*>(infile->Get("h_psdstkmatch_efficiency_accepted_loose_xtrl"));
        auto h_psdstkmatch_efficiency_total_loose_xtrl = static_cast<TH1D*>(infile->Get("h_psdstkmatch_efficiency_total_loose_xtrl"));

        h_psdstkmatch_efficiency_accepted_tight_xtrl->SetDirectory(0);
        h_psdstkmatch_efficiency_total_tight_xtrl->SetDirectory(0);
        h_psdstkmatch_efficiency_accepted_bdt->SetDirectory(0);
        h_psdstkmatch_efficiency_total_bdt->SetDirectory(0);
        
        h_psdstkmatch_efficiency_accepted_loose_xtrl->SetDirectory(0);
        h_psdstkmatch_efficiency_total_loose_xtrl->SetDirectory(0);

        // PSD charge histos
        auto h_psdcharge_efficiency_accepted_tight_xtrl = static_cast<TH1D*>(infile->Get("h_psdcharge_efficiency_accepted_tight_xtrl"));
        auto h_psdcharge_efficiency_total_tight_xtrl = static_cast<TH1D*>(infile->Get("h_psdcharge_efficiency_total_tight_xtrl"));
        auto h_psdcharge_efficiency_accepted_bdt = static_cast<TH1D*>(infile->Get("h_psdcharge_efficiency_accepted_bdt"));
        auto h_psdcharge_efficiency_total_bdt = static_cast<TH1D*>(infile->Get("h_psdcharge_efficiency_total_bdt"));

        auto h_psdcharge_efficiency_accepted_loose_xtrl = static_cast<TH1D*>(infile->Get("h_psdcharge_efficiency_accepted_loose_xtrl"));
        auto h_psdcharge_efficiency_total_loose_xtrl = static_cast<TH1D*>(infile->Get("h_psdcharge_efficiency_total_loose_xtrl"));

        h_psdcharge_efficiency_accepted_tight_xtrl->SetDirectory(0);
        h_psdcharge_efficiency_total_tight_xtrl->SetDirectory(0);
        h_psdcharge_efficiency_accepted_bdt->SetDirectory(0);
        h_psdcharge_efficiency_total_bdt->SetDirectory(0);
        
        h_psdcharge_efficiency_accepted_loose_xtrl->SetDirectory(0);
        h_psdcharge_efficiency_total_loose_xtrl->SetDirectory(0);

        infile->Close();

        // Build efficiencies
        std::unique_ptr<TEfficiency> trigger_eff_het_xtrl_tight;
        std::unique_ptr<TEfficiency> maxrms_eff_xtrl_tight;
        std::unique_ptr<TEfficiency> nbarlayer13_eff_xtrl_tight;
        std::unique_ptr<TEfficiency> maxrms_and_nbarlayer13_eff_xtrl_tight;
        std::unique_ptr<TEfficiency> track_selection_eff_xtrl_tight;
        std::unique_ptr<TEfficiency> psd_stk_match_eff_xtrl_tight;
        std::unique_ptr<TEfficiency> psd_charge_eff_xtrl_tight;

        std::unique_ptr<TEfficiency> trigger_eff_het_xtrl_loose;
        std::unique_ptr<TEfficiency> maxrms_eff_xtrl_loose;
        std::unique_ptr<TEfficiency> nbarlayer13_eff_xtrl_loose;
        std::unique_ptr<TEfficiency> maxrms_and_nbarlayer13_eff_xtrl_loose;
        std::unique_ptr<TEfficiency> track_selection_eff_xtrl_loose;
        std::unique_ptr<TEfficiency> psd_stk_match_eff_xtrl_loose;
        std::unique_ptr<TEfficiency> psd_charge_eff_xtrl_loose;

        std::unique_ptr<TEfficiency> trigger_eff_het_bdt;
        std::unique_ptr<TEfficiency> maxrms_eff_bdt;
        std::unique_ptr<TEfficiency> nbarlayer13_eff_bdt;
        std::unique_ptr<TEfficiency> maxrms_and_nbarlayer13_eff_bdt;
        std::unique_ptr<TEfficiency> track_selection_eff_bdt;
        std::unique_ptr<TEfficiency> psd_stk_match_eff_bdt;
        std::unique_ptr<TEfficiency> psd_charge_eff_bdt;

        if (TEfficiency::CheckConsistency(*h_trigger_efficiency_accepted_het_tight_xtrl, *h_trigger_efficiency_accepted_het_let_tight_xtrl))
            trigger_eff_het_xtrl_tight = std::make_unique<TEfficiency>(*h_trigger_efficiency_accepted_het_tight_xtrl, *h_trigger_efficiency_accepted_het_let_tight_xtrl);

        if (TEfficiency::CheckConsistency(*h_maxrms_efficiency_accepted_tight_xtrl, *h_maxrms_efficiency_total_tight_xtrl))
            maxrms_eff_xtrl_tight = std::make_unique<TEfficiency>(*h_maxrms_efficiency_accepted_tight_xtrl, *h_maxrms_efficiency_total_tight_xtrl);

        if (TEfficiency::CheckConsistency(*h_nbarlayer13_efficiency_accepted_tight_xtrl, *h_nbarlayer13_efficiency_total_tight_xtrl))
            nbarlayer13_eff_xtrl_tight = std::make_unique<TEfficiency>(*h_nbarlayer13_efficiency_accepted_tight_xtrl, *h_nbarlayer13_efficiency_total_tight_xtrl);

        if (TEfficiency::CheckConsistency(*h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl, *h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl))
            maxrms_and_nbarlayer13_eff_xtrl_tight = std::make_unique<TEfficiency>(*h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl, *h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl);

        if (TEfficiency::CheckConsistency(*h_track_efficiency_accepted_tight_xtrl, *h_track_efficiency_total_tight_xtrl))
            track_selection_eff_xtrl_tight = std::make_unique<TEfficiency>(*h_track_efficiency_accepted_tight_xtrl, *h_track_efficiency_total_tight_xtrl);

        if (TEfficiency::CheckConsistency(*h_psdstkmatch_efficiency_accepted_tight_xtrl, *h_psdstkmatch_efficiency_total_tight_xtrl))
            psd_stk_match_eff_xtrl_tight = std::make_unique<TEfficiency>(*h_psdstkmatch_efficiency_accepted_tight_xtrl, *h_psdstkmatch_efficiency_total_tight_xtrl);
        
        if (TEfficiency::CheckConsistency(*h_psdcharge_efficiency_accepted_tight_xtrl, *h_psdcharge_efficiency_total_tight_xtrl))
            psd_charge_eff_xtrl_tight = std::make_unique<TEfficiency>(*h_psdcharge_efficiency_accepted_tight_xtrl, *h_psdcharge_efficiency_total_tight_xtrl);

        if (TEfficiency::CheckConsistency(*h_trigger_efficiency_accepted_het_loose_xtrl, *h_trigger_efficiency_accepted_het_let_loose_xtrl))
            trigger_eff_het_xtrl_loose = std::make_unique<TEfficiency>(*h_trigger_efficiency_accepted_het_loose_xtrl, *h_trigger_efficiency_accepted_het_let_loose_xtrl);

        if (TEfficiency::CheckConsistency(*h_maxrms_efficiency_accepted_loose_xtrl, *h_maxrms_efficiency_total_loose_xtrl))
            maxrms_eff_xtrl_loose = std::make_unique<TEfficiency>(*h_maxrms_efficiency_accepted_loose_xtrl, *h_maxrms_efficiency_total_loose_xtrl);

        if (TEfficiency::CheckConsistency(*h_nbarlayer13_efficiency_accepted_loose_xtrl, *h_nbarlayer13_efficiency_total_loose_xtrl))
            nbarlayer13_eff_xtrl_loose = std::make_unique<TEfficiency>(*h_nbarlayer13_efficiency_accepted_loose_xtrl, *h_nbarlayer13_efficiency_total_loose_xtrl);

        if (TEfficiency::CheckConsistency(*h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl, *h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl))
            maxrms_and_nbarlayer13_eff_xtrl_loose = std::make_unique<TEfficiency>(*h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl, *h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl);

        if (TEfficiency::CheckConsistency(*h_track_efficiency_accepted_loose_xtrl, *h_track_efficiency_total_loose_xtrl))
            track_selection_eff_xtrl_loose = std::make_unique<TEfficiency>(*h_track_efficiency_accepted_loose_xtrl, *h_track_efficiency_total_loose_xtrl);

        if (TEfficiency::CheckConsistency(*h_psdstkmatch_efficiency_accepted_loose_xtrl, *h_psdstkmatch_efficiency_total_loose_xtrl))
            psd_stk_match_eff_xtrl_loose = std::make_unique<TEfficiency>(*h_psdstkmatch_efficiency_accepted_loose_xtrl, *h_psdstkmatch_efficiency_total_loose_xtrl);
        
        if (TEfficiency::CheckConsistency(*h_psdcharge_efficiency_accepted_loose_xtrl, *h_psdcharge_efficiency_total_loose_xtrl))
            psd_charge_eff_xtrl_loose = std::make_unique<TEfficiency>(*h_psdcharge_efficiency_accepted_loose_xtrl, *h_psdcharge_efficiency_total_loose_xtrl);

        if (TEfficiency::CheckConsistency(*h_trigger_efficiency_accepted_het_bdt, *h_trigger_efficiency_accepted_het_let_bdt))
            trigger_eff_het_bdt = std::make_unique<TEfficiency>(*h_trigger_efficiency_accepted_het_bdt, *h_trigger_efficiency_accepted_het_let_bdt);
        
        if (TEfficiency::CheckConsistency(*h_maxrms_efficiency_accepted_bdt, *h_maxrms_efficiency_total_bdt))
            maxrms_eff_bdt = std::make_unique<TEfficiency>(*h_maxrms_efficiency_accepted_bdt, *h_maxrms_efficiency_total_bdt);

        if (TEfficiency::CheckConsistency(*h_nbarlayer13_efficiency_accepted_bdt, *h_nbarlayer13_efficiency_total_bdt))
            nbarlayer13_eff_bdt = std::make_unique<TEfficiency>(*h_nbarlayer13_efficiency_accepted_bdt, *h_nbarlayer13_efficiency_total_bdt);

        if (TEfficiency::CheckConsistency(*h_maxrms_and_nbarlayer13_efficiency_accepted_bdt, *h_maxrms_and_nbarlayer13_efficiency_total_bdt))
            maxrms_and_nbarlayer13_eff_bdt = std::make_unique<TEfficiency>(*h_maxrms_and_nbarlayer13_efficiency_accepted_bdt, *h_maxrms_and_nbarlayer13_efficiency_total_bdt);

        if (TEfficiency::CheckConsistency(*h_track_efficiency_accepted_bdt, *h_track_efficiency_total_bdt))
            track_selection_eff_bdt = std::make_unique<TEfficiency>(*h_track_efficiency_accepted_bdt, *h_track_efficiency_total_bdt);

        if (TEfficiency::CheckConsistency(*h_psdstkmatch_efficiency_accepted_bdt, *h_psdstkmatch_efficiency_total_bdt))
            psd_stk_match_eff_bdt = std::make_unique<TEfficiency>(*h_psdstkmatch_efficiency_accepted_bdt, *h_psdstkmatch_efficiency_total_bdt);

        if (TEfficiency::CheckConsistency(*h_psdcharge_efficiency_accepted_bdt, *h_psdcharge_efficiency_total_bdt))
            psd_charge_eff_bdt = std::make_unique<TEfficiency>(*h_psdcharge_efficiency_accepted_bdt, *h_psdcharge_efficiency_total_bdt);

        trigger_eff_het_xtrl_tight                  ->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_eff_xtrl_tight                       ->SetStatisticOption(TEfficiency::kBUniform);
        nbarlayer13_eff_xtrl_tight                  ->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_and_nbarlayer13_eff_xtrl_tight       ->SetStatisticOption(TEfficiency::kBUniform);
        track_selection_eff_xtrl_tight              ->SetStatisticOption(TEfficiency::kBUniform);
        psd_stk_match_eff_xtrl_tight                ->SetStatisticOption(TEfficiency::kBUniform);
        psd_charge_eff_xtrl_tight                   ->SetStatisticOption(TEfficiency::kBUniform);
        trigger_eff_het_xtrl_loose                  ->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_eff_xtrl_loose                       ->SetStatisticOption(TEfficiency::kBUniform);
        nbarlayer13_eff_xtrl_loose                  ->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_and_nbarlayer13_eff_xtrl_loose       ->SetStatisticOption(TEfficiency::kBUniform);
        track_selection_eff_xtrl_loose              ->SetStatisticOption(TEfficiency::kBUniform);
        psd_stk_match_eff_xtrl_loose                ->SetStatisticOption(TEfficiency::kBUniform);
        psd_charge_eff_xtrl_loose                   ->SetStatisticOption(TEfficiency::kBUniform);

        trigger_eff_het_bdt                         ->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_eff_bdt                              ->SetStatisticOption(TEfficiency::kBUniform);
        nbarlayer13_eff_bdt                         ->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_and_nbarlayer13_eff_bdt              ->SetStatisticOption(TEfficiency::kBUniform);
        track_selection_eff_bdt                     ->SetStatisticOption(TEfficiency::kBUniform);
        psd_stk_match_eff_bdt                       ->SetStatisticOption(TEfficiency::kBUniform);
        psd_charge_eff_bdt                          ->SetStatisticOption(TEfficiency::kBUniform);

        trigger_eff_het_xtrl_tight                  ->SetName("trigger_eff_het_xtrl_tight");
        maxrms_eff_xtrl_tight                       ->SetName("maxrms_eff_xtrl_tight");
        nbarlayer13_eff_xtrl_tight                  ->SetName("nbarlayer13_eff_xtrl_tight");
        maxrms_and_nbarlayer13_eff_xtrl_tight       ->SetName("maxrms_and_nbarlayer13_eff_xtrl_tight");
        track_selection_eff_xtrl_tight              ->SetName("track_selection_eff_xtrl_tight");
        psd_stk_match_eff_xtrl_tight                ->SetName("psd_stk_match_eff_xtrl_tight");
        psd_charge_eff_xtrl_tight                   ->SetName("psd_charge_eff_xtrl_tight");
        trigger_eff_het_xtrl_loose                  ->SetName("trigger_eff_het_xtrl_loose");
        maxrms_eff_xtrl_loose                       ->SetName("maxrms_eff_xtrl_loose");
        nbarlayer13_eff_xtrl_loose                  ->SetName("nbarlayer13_eff_xtrl_loose");
        maxrms_and_nbarlayer13_eff_xtrl_loose       ->SetName("maxrms_and_nbarlayer13_eff_xtrl_loose");
        track_selection_eff_xtrl_loose              ->SetName("track_selection_eff_xtrl_loose");
        psd_stk_match_eff_xtrl_loose                ->SetName("psd_stk_match_eff_xtrl_loose");
        psd_charge_eff_xtrl_loose                   ->SetName("psd_charge_eff_xtrl_loose");

        trigger_eff_het_xtrl_tight                  ->SetTitle("trigger_eff_het_xtrl_tight");
        maxrms_eff_xtrl_tight                       ->SetTitle("maxrms_eff_xtrl_tight");
        nbarlayer13_eff_xtrl_tight                  ->SetTitle("nbarlayer13_eff_xtrl_tight");
        maxrms_and_nbarlayer13_eff_xtrl_tight       ->SetTitle("maxrms_and_nbarlayer13_eff_xtrl_tight");
        track_selection_eff_xtrl_tight              ->SetTitle("track_selection_eff_xtrl_tight");
        psd_stk_match_eff_xtrl_tight                ->SetTitle("psd_stk_match_eff_xtrl_tight");
        psd_charge_eff_xtrl_tight                   ->SetTitle("psd_charge_eff_xtrl_tight");
        trigger_eff_het_xtrl_loose                  ->SetTitle("trigger_eff_het_xtrl_loose");
        maxrms_eff_xtrl_loose                       ->SetTitle("maxrms_eff_xtrl_loose");
        nbarlayer13_eff_xtrl_loose                  ->SetTitle("nbarlayer13_eff_xtrl_loose");
        maxrms_and_nbarlayer13_eff_xtrl_loose       ->SetTitle("maxrms_and_nbarlayer13_eff_xtrl_loose");
        track_selection_eff_xtrl_loose              ->SetTitle("track_selection_eff_xtrl_loose");
        psd_stk_match_eff_xtrl_loose                ->SetTitle("psd_stk_match_eff_xtrl_loose");
        psd_charge_eff_xtrl_loose                   ->SetTitle("psd_charge_eff_xtrl_loose");

        trigger_eff_het_bdt                         ->SetName("trigger_eff_het_bdt");
        maxrms_eff_bdt                              ->SetName("maxrms_eff_bdt");
        nbarlayer13_eff_bdt                         ->SetName("nbarlayer13_eff_bdt");
        maxrms_and_nbarlayer13_eff_bdt              ->SetName("maxrms_and_nbarlayer13_eff_bdt");
        track_selection_eff_bdt                     ->SetName("track_selection_eff_bdt");
        psd_stk_match_eff_bdt                       ->SetName("psd_stk_match_eff_bdt");
        psd_charge_eff_bdt                          ->SetName("psd_charge_eff_bdt");

        trigger_eff_het_bdt                         ->SetTitle("trigger_eff_het_bdt");
        maxrms_eff_bdt                              ->SetTitle("maxrms_eff_bdt");
        nbarlayer13_eff_bdt                         ->SetTitle("nbarlayer13_eff_bdt");
        maxrms_and_nbarlayer13_eff_bdt              ->SetTitle("maxrms_and_nbarlayer13_eff_bdt");
        track_selection_eff_bdt                     ->SetTitle("track_selection_eff_bdt");
        psd_stk_match_eff_bdt                       ->SetTitle("psd_stk_match_eff_bdt");
        psd_charge_eff_bdt                          ->SetTitle("psd_charge_eff_bdt");

        // Write output TFile
        TFile* outfile = TFile::Open(output_file, "RECREATE");
        if (outfile->IsZombie()) {
            std::cerr << "Error writing output ROOT file [" << output_file << "] \n\n";
            exit(100);
        }

        outfile->mkdir("efficiencies");
        outfile->cd("efficiencies");

        trigger_eff_het_xtrl_tight                  ->Write();
        maxrms_eff_xtrl_tight                       ->Write();
        nbarlayer13_eff_xtrl_tight                  ->Write();
        maxrms_and_nbarlayer13_eff_xtrl_tight       ->Write();
        track_selection_eff_xtrl_tight              ->Write();
        psd_stk_match_eff_xtrl_tight                ->Write();
        psd_charge_eff_xtrl_tight                   ->Write();
        trigger_eff_het_xtrl_loose                  ->Write();
        maxrms_eff_xtrl_loose                       ->Write();
        nbarlayer13_eff_xtrl_loose                  ->Write();
        maxrms_and_nbarlayer13_eff_xtrl_loose       ->Write();
        track_selection_eff_xtrl_loose              ->Write();
        psd_stk_match_eff_xtrl_loose                ->Write();
        psd_charge_eff_xtrl_loose                   ->Write();

        trigger_eff_het_bdt                         ->Write();
        maxrms_eff_bdt                              ->Write();
        nbarlayer13_eff_bdt                         ->Write();
        maxrms_and_nbarlayer13_eff_bdt              ->Write();
        track_selection_eff_bdt                     ->Write();
        psd_stk_match_eff_bdt                       ->Write();
        psd_charge_eff_bdt                          ->Write();

        outfile->mkdir("histos");
        outfile->cd("histos");

        h_trigger_efficiency_accepted_het_tight_xtrl                ->Write();
        h_trigger_efficiency_accepted_het_let_tight_xtrl            ->Write();
        h_maxrms_efficiency_accepted_tight_xtrl                     ->Write();
        h_maxrms_efficiency_total_tight_xtrl                        ->Write();
        h_nbarlayer13_efficiency_accepted_tight_xtrl                ->Write();
        h_nbarlayer13_efficiency_total_tight_xtrl                   ->Write();
        h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl     ->Write();
        h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl        ->Write();
        h_track_efficiency_accepted_tight_xtrl                      ->Write();
        h_track_efficiency_total_tight_xtrl                         ->Write();
        h_psdstkmatch_efficiency_accepted_tight_xtrl                ->Write();
        h_psdstkmatch_efficiency_total_tight_xtrl                   ->Write();
        h_psdcharge_efficiency_accepted_tight_xtrl                  ->Write();
        h_psdcharge_efficiency_total_tight_xtrl                     ->Write();
        h_trigger_efficiency_accepted_het_loose_xtrl                ->Write();
        h_trigger_efficiency_accepted_het_let_loose_xtrl            ->Write();
        h_maxrms_efficiency_accepted_loose_xtrl                     ->Write();
        h_maxrms_efficiency_total_loose_xtrl                        ->Write();
        h_nbarlayer13_efficiency_accepted_loose_xtrl                ->Write();
        h_nbarlayer13_efficiency_total_loose_xtrl                   ->Write();
        h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl     ->Write();
        h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl        ->Write();
        h_track_efficiency_accepted_loose_xtrl                      ->Write();
        h_track_efficiency_total_loose_xtrl                         ->Write();
        h_psdstkmatch_efficiency_accepted_loose_xtrl                ->Write();
        h_psdstkmatch_efficiency_total_loose_xtrl                   ->Write();
        h_psdcharge_efficiency_accepted_loose_xtrl                  ->Write();
        h_psdcharge_efficiency_total_loose_xtrl                     ->Write();

        h_trigger_efficiency_accepted_het_bdt                       ->Write();
        h_trigger_efficiency_accepted_het_let_bdt                   ->Write();
        h_maxrms_efficiency_accepted_bdt                            ->Write();
        h_maxrms_efficiency_total_bdt                               ->Write();
        h_nbarlayer13_efficiency_accepted_bdt                       ->Write();
        h_nbarlayer13_efficiency_total_bdt                          ->Write();
        h_maxrms_and_nbarlayer13_efficiency_accepted_bdt            ->Write();
        h_maxrms_and_nbarlayer13_efficiency_total_bdt               ->Write();
        h_track_efficiency_accepted_bdt                             ->Write();
        h_track_efficiency_total_bdt                                ->Write();
        h_psdstkmatch_efficiency_accepted_bdt                       ->Write();
        h_psdstkmatch_efficiency_total_bdt                          ->Write();
        h_psdcharge_efficiency_accepted_bdt                         ->Write();    
        h_psdcharge_efficiency_total_bdt                            ->Write();

        outfile->Close();

        if (verbose)
            std::cout << "\n\nOutput file has been written: [" << output_file << "]\n\n";
    }