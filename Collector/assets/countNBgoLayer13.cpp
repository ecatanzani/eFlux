#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TH2D.h"

inline std::vector<float> createLogBinning(const double eMin, const double eMax, const std::size_t n_bins) {
	std::vector<float> binning(n_bins + 1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);
	return binning;
}

struct energy_config {
    std::size_t n_bins;
    double min_event_energy;
    double max_event_energy;
};

inline std::string parse_config_file(const char* config_file_path) {
	std::ifstream input_file(config_file_path);
	if (!input_file.is_open())
	{
		std::cerr << "\nInput config file not found [" << config_file_path << "]\n\n";
		exit(100);
	}
	std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

inline void get_config_info(const std::string parsed_config, energy_config &config) {
	std::string tmp_str;
	std::istringstream input_stream(parsed_config);
	std::string::size_type sz;

	while (input_stream >> tmp_str) {
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> config.n_bins;
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			config.min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			config.max_event_energy = stod(tmp_str, &sz);
		}
	}
}

void countNBgoLayer13(
    const char* input_path, 
    const char* histo_name,
    const char* energy_config_file) {
    // Get energy config to obtain binning info
    energy_config config;
    get_config_info(parse_config_file(energy_config_file), config);
    auto energy_binning = createLogBinning(
		config.min_event_energy,
		config.max_event_energy,
		config.n_bins);
    
    // Read input file
    TFile *input_file = TFile::Open(input_path, "READ");
    if (!input_file->IsOpen()) {
        std::cerr << "\n\nError opening input file [" << input_path << "]\n\n";
        exit(100);
    }

    std::unique_ptr<TH2D> histo = std::unique_ptr<TH2D>(static_cast<TH2D*>(input_file->Get(histo_name)));

    double up_events {0};
    double tot_events {0};
    double bar_th {0};
    auto tmp_energy = [&energy_binning](int low, int high) -> double {return 0.5*(energy_binning[high]-energy_binning[low]);};

    for (int bX=1; bX<=histo->GetNbinsX(); ++bX) {
        bar_th = 8*log10(energy_binning[bX-1]) - 5;
        if (bar_th > 22)
            bar_th = 22;
        std::cout << "\nEnergy: " << energy_binning[bX-1] << " bars: " << (int)bar_th;
        auto proj = histo->ProjectionY("proj", bX, bX);
        up_events += proj->Integral(proj->FindBin(bar_th), proj->GetNbinsX()+1);
        tot_events += proj->Integral();
    }

    std::cout << "\n\nEvents above threshold: " << up_events;
    std::cout << "\nEvents: " << tot_events;
    std::cout << "\n\nFraction of rejected events: " << up_events/tot_events << std::endl;
}