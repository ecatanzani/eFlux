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

void extractClassifier(
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

        // Build histos
        auto bdt_binning = createLinearBinning(-10, 10, 100);
        auto xtrl_binning = createLinearBinning(0, 1000, 1000);
        auto h_classifier = fr.Histo2D<double, double>({"h_classifier", "XTRL/BDT classifiers; BDT; xtrl", (int)bdt_binning.size()-1, &bdt_binning[0], (int)xtrl_binning.size()-1, &xtrl_binning[0]}, "tmva_classifier", "xtrl");
        auto h_classifier_energy = fr.Histo3D<double, double, double>({"h_classifier_energy", "XTRL/BDT classifiers; BDT; xtrl; Corrected Energy [GeV]", (int)bdt_binning.size()-1, &bdt_binning[0], (int)xtrl_binning.size()-1, &xtrl_binning[0], (int)energy_binning.size()-1, &energy_binning[0]}, "tmva_classifier", "xtrl");

        std::vector<ROOT::RDF::RResultPtr<TH2D>> h_classifier_energy_bin((int)energy_binning.size()-1);
        for (unsigned int bIdx = 0; bIdx < energy_binning.size()-1; ++bIdx)
            h_classifier_energy_bin[bIdx] = fr.Histo2D<double, double>({
                (std::string("h_classifier_energy_bin_") + std::to_string(bIdx)).c_str(), 
                ("XTRL/BDT classifiers - energy bin " + std::to_string(bIdx+1) + std::string("; BDT; xtrl")).c_str(), 
                (int)bdt_binning.size()-1, &bdt_binning[0], (int)xtrl_binning.size()-1, &xtrl_binning[0]}, "tmva_classifier", "xtrl");

        TFile* outfile {TFile::Open(output_file, "RECREATE")};
        if (infile->IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        h_classifier->Write();
        h_classifier_energy->Write();
        for (unsigned int bIdx = 0; bIdx < energy_binning.size()-1; ++bIdx)
            h_classifier_energy_bin[bIdx]->Write();

        outfile->Close();
    }