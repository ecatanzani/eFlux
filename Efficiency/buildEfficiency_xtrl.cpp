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

std::string parse_input_file(const std::string input_list)
{
	std::ifstream input_file(input_list.c_str());
	if (!input_file.is_open())
	{
		std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
		exit(100);
	}
	std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
	input_file.close();
	return input_string;
}

std::shared_ptr<TChain> read_input_list(const char* input_list, const bool verbose) {
    std::istringstream input_stream(parse_input_file(std::string(input_list)));
    auto evtch = std::make_shared<TChain> (get_tree_name(input_stream.str()).c_str(), "DAMPE event tree");
    std::string tmp_str;
    while (input_stream >> tmp_str) {
        evtch->Add(tmp_str.c_str());
        if (verbose) std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
    return evtch;
}

void buildEfficiency(
    const char* input_file_list,
    const char* output_file,
    const char* energy_config_file,
    const bool verbose,
    const unsigned int threads=1) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        // Chain input files
        auto mytree = read_input_list(input_file_list, verbose);

        // Initialize RDF
        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame fr(*mytree);

        auto compute_xtrl = [](const double last_layer_energy_fraction, const double sumrms) -> double {
            return last_layer_energy_fraction != -999 ? 0.125e-6 * pow(sumrms, 4) * last_layer_energy_fraction : -999;
        };

        auto xtrl_loose_cut = [](const double bgo_total_energy_gev, const double xtrl) -> bool {
            auto xtrl_cut_value {3*log10(bgo_total_energy_gev/200) + 12};
            if (xtrl<xtrl_cut_value && xtrl!=-999)
                return true;
            else
                return false;
        };

        // Trigger histos
        auto h_trigger_efficiency_accepted_het_tight_xtrl = fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_trigger_efficiency_accepted_het_tight_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_trigger_efficiency_accepted_het_let_tight_xtrl = fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let)")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_trigger_efficiency_accepted_het_let_tight_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // MaxRMS histos
        auto h_maxrms_efficiency_accepted_tight_xtrl = fr.Filter("maxrms_efficiency_preselection==1 && maxrms_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_maxrms_efficiency_accepted_tight_xtrl", "HET Trigger + MaxRMS Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_maxrms_efficiency_total_tight_xtrl = fr.Filter("maxrms_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_maxrms_efficiency_total_tight_xtrl", "HET Trigger + MaxRMS", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // nbarlayer13 histos
        auto h_nbarlayer13_efficiency_accepted_tight_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1 && nbarlayer13_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_nbarlayer13_efficiency_accepted_tight_xtrl", "HET Trigger + nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_nbarlayer13_efficiency_total_tight_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_nbarlayer13_efficiency_total_tight_xtrl", "HET Trigger + nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        
        // MaxRms && nbarlayer13 histos
        auto h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && maxrms_and_nbarlayer13_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_accepted_tight_xtrl", "HET Trigger + maxrms & nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_total_tight_xtrl", "HET Trigger + maxrms & nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
       
        // Track Selection histos
        auto h_track_efficiency_accepted_tight_xtrl = fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_track_efficiency_accepted_tight_xtrl", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_track_efficiency_total_tight_xtrl = fr.Filter("track_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_track_efficiency_total_tight_xtrl", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // PSD-STK match histos
        auto h_psdstkmatch_efficiency_accepted_tight_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && psdstkmatch_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_psdstkmatch_efficiency_accepted_tight_xtrl", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_psdstkmatch_efficiency_total_tight_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_psdstkmatch_efficiency_total_tight_xtrl", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // PSD charge histos
        auto h_psdcharge_efficiency_accepted_tight_xtrl = fr.Filter("psdcharge_efficiency_preselection==1 && psdcharge_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_psdcharge_efficiency_accepted_tight_xtrl", "HET Trigger + PSD Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_psdcharge_efficiency_total_tight_xtrl = fr.Filter("psdcharge_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter("xtrl_evt<8.5 && xtrl_evt!= -999")
                                                .Histo1D({"h_psdcharge_efficiency_total_tight_xtrl", "HET Trigger + PSD Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});


        // Trigger histos
        auto h_trigger_efficiency_accepted_het_loose_xtrl = fr.Filter("trigger_efficiency_preselection==1 && trigger_efficiency_preselection_is_het==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_trigger_efficiency_accepted_het_loose_xtrl", "HET Trigger", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_trigger_efficiency_accepted_het_let_loose_xtrl = fr.Filter("trigger_efficiency_preselection==1 && (trigger_efficiency_preselection_is_het==1 || trigger_efficiency_preselection_is_let)")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_trigger_efficiency_accepted_het_let_loose_xtrl", "LET + HET Trigger", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // MaxRMS histos
        auto h_maxrms_efficiency_accepted_loose_xtrl = fr.Filter("maxrms_efficiency_preselection==1 && maxrms_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_maxrms_efficiency_accepted_loose_xtrl", "HET Trigger + MaxRMS Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_maxrms_efficiency_total_loose_xtrl = fr.Filter("maxrms_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_maxrms_efficiency_total_loose_xtrl", "HET Trigger + MaxRMS", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // nbarlayer13 histos
        auto h_nbarlayer13_efficiency_accepted_loose_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1 && nbarlayer13_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_nbarlayer13_efficiency_accepted_loose_xtrl", "HET Trigger + nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_nbarlayer13_efficiency_total_loose_xtrl = fr.Filter("nbarlayer13_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_nbarlayer13_efficiency_total_loose_xtrl", "HET Trigger + nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
                                            
        // MaxRms && nbarlayer13 histos
        auto h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && maxrms_and_nbarlayer13_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_accepted_loose_xtrl", "HET Trigger + maxrms & nbarlayer13 Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl = fr.Filter("maxrms_and_nbarlayer13_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_maxrms_and_nbarlayer13_efficiency_total_loose_xtrl", "HET Trigger + maxrms & nbarlayer13", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // Track Selection histos
        auto h_track_efficiency_accepted_loose_xtrl = fr.Filter("track_efficiency_preselection==1 && track_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_track_efficiency_accepted_loose_xtrl", "HET Trigger + Track Selection Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_track_efficiency_total_loose_xtrl = fr.Filter("track_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_track_efficiency_total_loose_xtrl", "HET Trigger + Track Selection", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // PSD-STK match histos
        auto h_psdstkmatch_efficiency_accepted_loose_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && psdstkmatch_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_psdstkmatch_efficiency_accepted_loose_xtrl", "HET Trigger + PSD-STK Match Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_psdstkmatch_efficiency_total_loose_xtrl = fr.Filter("psdstkmatch_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_psdstkmatch_efficiency_total_loose_xtrl", "HET Trigger + PSD-STK Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // PSD charge histos
         auto h_psdcharge_efficiency_accepted_loose_xtrl = fr.Filter("psdcharge_efficiency_preselection==1 && psdcharge_efficiency_preselection_accepted==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_psdcharge_efficiency_accepted_loose_xtrl", "HET Trigger + PSD Charge Accepted", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});
        auto h_psdcharge_efficiency_total_loose_xtrl = fr.Filter("psdcharge_efficiency_preselection==1 && HET_trigger==1")
                                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                                .Filter(xtrl_loose_cut, {"energy", "xtrl_evt"})
                                                .Histo1D({"h_psdcharge_efficiency_total_loose_xtrl", "HET Trigger + PSD Charge Match", energy_nbins, &energy_binning[0]}, {"corr_energy_gev"});

        // XTRL vs STK cosine histo
        auto h_xtrl_stk_cosine = fr.Filter("HET_trigger==1 && evtfilter_correct_bgo_reco==1")
                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                .Histo2D({"h_xtrl_stk_cosine", "XTRL vs STK direction; cosine STK direction #cos(#theta); xtrl", 100, 0, 1, 300, 0, 300}, "STK_bestTrack_costheta", "xtrl_evt");
        
        auto h_xtrl_bgo_cosine = fr.Filter("HET_trigger==1 && evtfilter_correct_bgo_reco==1")
                                .Define("corr_energy_gev", "energy_corr * 0.001")
                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                .Define("xtrl_evt", compute_xtrl, {"fracLast_13", "sumRms"})
                                .Histo2D({"h_xtrl_bgo_cosine", "XTRL vs BGO direction; cosine BGO direction #cos(#theta); xtrl", 100, 0, 1, 300, 0, 300}, "bgorec_cosine", "xtrl_evt");

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

        h_xtrl_stk_cosine                                           ->Write();
        h_xtrl_bgo_cosine                                           ->Write();

        outfile->Close();

        if (verbose)
            std::cout << "\n\nOutput file has been written: [" << output_file << "]\n\n";
    }