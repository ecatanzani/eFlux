#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TF1.h"
#include "TKey.h"
#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"
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

void mcShift(
    const char* input_file,
    const char* output_file,
    const char* energy_config_file,
    const bool verbose,
    const bool mc = false,
    const unsigned int threads=1) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

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

        std::vector<ROOT::RDF::RResultPtr<TH1D>> h_classifier_bin(energy_nbins);

        const double signal_spectral_index = -3;
        auto get_weight = [signal_spectral_index, &energy_binning] (const double energy_gev) -> double {
            return std::pow(energy_gev, signal_spectral_index -1)*std::pow(energy_binning[0], 2);
        };

        // Build the tmva classifier distribution for each bin
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
            auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };
            h_classifier_bin[bin_idx-1] = !mc ? fr.Filter(bin_filter, {"energy_bin"}).Histo1D<double>("tmva_classifier") : 
                                                fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("simu_energy_gev", [](const double energy) -> double {return energy*0.001;}, {"simu_energy"})
                                                    .Define("evt_w", get_weight, {"simu_energy_gev"})
                                                    .Histo1D<double, double>("tmva_classifier", "evt_w");
        }

        std::vector<TF1> gaus_fit_1(energy_nbins), gaus_fit_2(energy_nbins);
        
        // First fit
        double mean {0}, rms {0};
        for (auto it=std::begin(h_classifier_bin); it != std::end(h_classifier_bin); ++it) {
            mean = it->GetPtr()->GetMean();
            rms = it->GetPtr()->GetRMS();
            TF1 gaus_fit("gaus_fit", "gaus", mean-rms, mean+rms);
            it->GetPtr()->Fit(&gaus_fit, "Q");
            gaus_fit_1[std::distance(std::begin(h_classifier_bin), it)] = gaus_fit;
        }

        // Second fit
        double tf1_mean {0}, tf1_rms {0};
        for (auto it=std::begin(h_classifier_bin); it != std::end(h_classifier_bin); ++it) {
            tf1_mean = gaus_fit_1[std::distance(std::begin(h_classifier_bin), it)].GetParameter(1);
            tf1_rms = gaus_fit_1[std::distance(std::begin(h_classifier_bin), it)].GetParameter(2);
            TF1 gaus_fit("gaus_fit", "gaus", tf1_mean-tf1_rms, tf1_mean+tf1_rms);
            it->GetPtr()->Fit(&gaus_fit, "Q");
            gaus_fit_2[std::distance(std::begin(h_classifier_bin), it)] = gaus_fit;
        }

        std::vector<double> mean_shift(energy_nbins, 0), rms_shift(energy_nbins, 0);
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
            mean_shift[bin_idx-1] = gaus_fit_2[bin_idx-1].GetParameter(1);
            rms_shift[bin_idx-1] = gaus_fit_2[bin_idx-1].GetParameter(2);
        }

        std::vector<double> bins(energy_nbins);
        std::iota(std::begin(bins), std::end(bins), 1);
        TGraph gr_shift_mean(energy_nbins, &bins[0], &mean_shift[0]);
        gr_shift_mean.SetName("gr_shift_mean");
        gr_shift_mean.GetXaxis()->SetTitle("Energy Bin");
        gr_shift_mean.GetYaxis()->SetTitle("Mean Shift");

        TGraph gr_shift_rms(energy_nbins, &bins[0], &rms_shift[0]);
        gr_shift_mean.SetName("gr_shift_rms");
        gr_shift_mean.GetXaxis()->SetTitle("Energy Bin");
        gr_shift_mean.GetYaxis()->SetTitle("RMS Shift");

        TFile outfile(output_file, "RECREATE");
        if (outfile.IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        for (int bidx = 0; bidx < energy_nbins; ++bidx) {
            auto tmp_dir_name = std::string("energybin_") + std::to_string(bidx + 1);
            outfile.mkdir(tmp_dir_name.c_str());
            outfile.cd(tmp_dir_name.c_str());
            h_classifier_bin[bidx]->Write();
            gaus_fit_1[bidx].Write();
            gaus_fit_2[bidx].Write();
        }

        gr_shift_mean.Write();
        gr_shift_rms.Write();

        outfile.Close();
    }