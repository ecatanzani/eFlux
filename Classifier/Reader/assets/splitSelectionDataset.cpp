#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
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

inline void split_energy_bins(
    std::shared_ptr<TTree> cumulative_tree, 
    const std::string wd, 
    const int energy_nbins,
    const unsigned int threads,
    const bool verbose) {

        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame _data_fr(*cumulative_tree);

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
            auto bin_filter = [bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };

            std::string fname;
            if (wd.empty()) fname = cumulative_tree->GetName() + std::string("_energybin_") + std::to_string(bin_idx) + std::string(".root");
            else fname = wd + std::string("/") + cumulative_tree->GetName() + std::string("_energybin_") + std::to_string(bin_idx) + std::string(".root");
            auto tmp_entries = *(_data_fr.Filter(bin_filter, {"energy_bin"}).Count());
            if (tmp_entries) {
                if (verbose) std::cout << "\nOutput TFile has been written [" << fname << "]\t entries: " << tmp_entries;
                _data_fr.Filter(bin_filter, {"energy_bin"}).Snapshot(cumulative_tree->GetName(), fname.c_str());
            }
        }
}

inline std::shared_ptr<TTree> getTreeFromFile(const char* input_file, const char* tree_name, const bool verbose) {
    std::shared_ptr<TFile> input = std::make_shared<TFile>(input_file, "READ");
    if (!input->IsOpen()) {
        std::cerr << "\n\nInput file not found [" << input_file << "]\n\n";
        exit(100);
    }

    std::shared_ptr<TTree> tree(static_cast<TTree*>(input->Get(tree_name)));
    if (verbose) std::cout << "\n\nExtracting from input file ... [" << tree_name << "]\n\n";
    return tree;
}

void splitSelectionDataset(
    const char* input_file,
    const char* output_directory,
    const char* energy_config_file,
    const bool verbose=true,
    const char* tree_name="total_tree",
    const unsigned int threads=1) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto n_energy_bins = static_cast<int>(energy_binning.size()) - 1;

        // Extract tree from file
        //auto tree = getTreeFromFile(input_file, tree_name, verbose);
        std::shared_ptr<TFile> input = std::make_shared<TFile>(input_file, "READ");
        if (!input->IsOpen()) {
            std::cerr << "\n\nInput file not found [" << input_file << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TTree> tree(static_cast<TTree*>(input->Get(tree_name)));
        if (verbose) std::cout << "\n\nExtracting from input file ... [" << tree_name << "]\n\n";

        // Split tree into energy bins
        //split_energy_bins(tree, output_directory, n_energy_bins, threads, verbose);
        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame _data_fr(*tree);

        for (int bin_idx = 1; bin_idx <= n_energy_bins; ++bin_idx) {
            auto bin_filter = [bin_idx](int energy_bin) -> bool { return energy_bin == bin_idx; };

            std::string fname = output_directory + std::string("/") + tree->GetName() + std::string("_energybin_") + std::to_string(bin_idx) + std::string(".root");
            auto tmp_entries = *(_data_fr.Filter(bin_filter, {"energy_bin"}).Count());
            if (tmp_entries) {
                if (verbose) std::cout << "\nOutput TFile has been written [" << fname << "]\t entries: " << tmp_entries;
                _data_fr.Filter(bin_filter, {"energy_bin"}).Snapshot(tree->GetName(), fname.c_str());
            }
        }
    }