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

void extractElectronsCountHisto(
    const char* input_file, 
    const char* output_file,
    const char* energy_config_file,
    const bool verbose,
    const unsigned int threads=1) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);

        // Read input file
        TFile* infile {TFile::Open(input_file, "READ")};
        if (infile->IsZombie()) {
            std::cerr << "\n\nError opening input file [" << input_file << "]" << std::endl;
            exit(100);
        }

        TIter nextkey(infile->GetListOfKeys());
        TKey *key {nullptr};
        std::unique_ptr<TTree> mytree;
        while ((key=static_cast<TKey*>(nextkey()))) 
        {
            TObject *obj {key->ReadObj()};
            if (obj->IsA()->InheritsFrom(TTree::Class())) {
                mytree = std::unique_ptr<TTree>(static_cast<TTree*>(obj));
                break;
            }
        }

        if (verbose)
            std::cout << "\nFound TTree in input file [" << mytree->GetName() << "]";

        // Initialize RDF
        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame fr(*mytree);

        // Build histo
        auto h_electron = fr.Define("energy_corr_gev", "energy_corr*0.001").Histo1D<double>(
            {"h_electron", "Electrons Energy Distribution", (int)energy_binning.size()-1, &energy_binning[0]},
            "energy_corr_gev");

        TFile *outfile {TFile::Open(output_file, "RECREATE")};
        if (infile->IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        h_electron->Write();
        outfile->Close();

        if (verbose)
            std::cout << "\nOutput file has been written [" << output_file << "]\n\n";
    }