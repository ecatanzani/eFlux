#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TKey.h"
#include "TH1D.h"
#include "TFile.h"
#include "TGraphErrors.h"

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

inline std::shared_ptr<TH1D> extractHistoFromFile(const char* file, const bool verbose) {
    TFile* infile {TFile::Open(file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key {nullptr};
    std::shared_ptr<TH1D> myhisto;
    while ((key=static_cast<TKey*>(nextkey())))  {
        TObject *obj {key->ReadObj()};
        if (obj->IsA()->InheritsFrom(TH1D::Class())) {
            myhisto = std::shared_ptr<TH1D>(static_cast<TH1D*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound TTree in input file [" << myhisto->GetName() << "] --> [" << file << "]";
    
    return myhisto;
}

void fromTHtoGR(
    const char* infile, 
    const char* outfile,
    const double energy_th_gev,
    //const char* energy_config_file, 
    const bool verbose) {

        /*
        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        
        // Build new restricted binning
        std::vector<double> restricted_energy_binning;

        for (auto&& elm : energy_binning)
            if (elm > energy_th_gev)
                restricted_energy_binning.push_back(elm);
        */
        
        auto h_acceptance = extractHistoFromFile(infile, verbose);

        std::vector<double> energies, acceptance;
        std::vector<double> err_energies, err_acceptance;

        for (int b_idx=1; b_idx<=h_acceptance->GetNbinsX(); ++b_idx) {
            auto bin_center = h_acceptance->GetBinCenter(b_idx);
            if (bin_center > energy_th_gev) {
                energies.push_back(bin_center);
                acceptance.push_back(h_acceptance->GetBinContent(b_idx));
                err_energies.push_back(0);
                err_acceptance.push_back(h_acceptance->GetBinError(b_idx));
            }
        }
        
        TGraphErrors gr(energies.size(), &energies[0], &acceptance[0], &err_energies[0], &err_acceptance[0]);
        gr.SetName("gr_acceptance");
        gr.GetXaxis()->SetTitle("Energy [GeV]");
        gr.GetYaxis()->SetTitle("Acceptance [m^{2} sr]");

        TFile* outfile_ptr {TFile::Open(outfile, "RECREATE")};
        if (outfile_ptr->IsZombie()) {
            std::cerr << "\n\nError opening output file [" << outfile << "]" << std::endl;
            exit(100);
        }

        gr.Write();

        outfile_ptr->Close();

    }