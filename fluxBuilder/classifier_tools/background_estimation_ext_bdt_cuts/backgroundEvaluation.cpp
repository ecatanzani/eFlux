#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TKey.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"


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

void backgroundEvaluation(
    const char* bdt_best_cut,
    const char* data_bdt_profiles,
    const char* energy_config_file) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        // Extract external BDT cut values
        auto bdt_cuts = extract_bdt_cut_values(input_args.bdt_cut_values, energy_nbins);

        // Read proton fit file
        TFile input_proton_fit_file(data_bdt_profiles, "READ");
        if (!input_proton_fit_file.IsOpen()) {
            std::cout << "ERROR opening input proton fit file: [" << input_proton_fit_file << "]\n\n";
            exit(100);
        }

        // Evaluate the background contamination in each energy bin
        for (unsigned int energy_bin_idx = 0; energy_bin_idx < energy_nbins; ++energy_bin_idx) {
            
            // Extract TMVA classifier distribution and proton linear fit for each energy bin
            std::string h_name = "energybin_" + std::to_string(energy_bin_idx+1) + "/h_classifier_bin_" + std::to_string(energy_bin_idx+1);
            std::string f_name = "energybin_" + std::to_string(energy_bin_idx+1) + "/proton_linear_fit_" + std::to_string(energy_bin_idx+);;
            auto h_data = static_cast<TH1D*>(input_proton_fit_file.Get(h_name.c_str()));
            auto background_fit = static_cast<TF1*>(input_proton_fit_file.Get(f_name.c_str()));

            double h_data_bin_width = h_data->GetBinWidth(1);   // The histos have uniform binning. The binning is used to rescale the integral of the TF1
            
            double bacground_contamination = background_fit->Integral(bdt_cuts[energy_bin_idx], 1.0)/h_data_bin_width;
        }
        

    }