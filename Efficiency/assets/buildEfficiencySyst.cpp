#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>
#include <string>
#include <tuple>

#include "TH1D.h"
#include "TAxis.h"
#include "TKey.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

struct energy_config {
    std::size_t n_bins;
    double min_event_energy {-999};
    double max_event_energy {-999};
    std::vector<double> energy_binning;
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
	
    int introduction {4};
	int elm_in_row {5};
    std::vector<std::string> row;

	int column_counter {0};
    int line_counter {0};

    while (input_stream >> tmp_str) {

        if (!line_counter) {
			// This is the first line... we are not interested in it
			++column_counter;
			if (column_counter==introduction) {
				++line_counter;
				column_counter = 0;
			}
		}
		else {
			// This is a general line...
			row.push_back(tmp_str);
			++column_counter;
			
			if (column_counter == elm_in_row) {

				// The row of the binning has been completed... let's extract the info
				if (line_counter==1) 
					config_pars.energy_binning.push_back(stod(row[2], &sz));
				config_pars.energy_binning.push_back(stod(row.back(), &sz));

				// Reset
				column_counter = 0;
				++line_counter;
				row.clear();
			}
		}
    }
    
    return config_pars;
}

inline std::vector<double> parse_energy_config(const char* config_file) {
    return get_config_info(parse_config_file(config_file)).energy_binning;
}

inline std::shared_ptr<TGraphErrors> extractGraphFromFile(const char* file, const bool verbose) {
    TFile* infile {TFile::Open(file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key {nullptr};
    std::shared_ptr<TGraphErrors> mygraph;
    while ((key=static_cast<TKey*>(nextkey())))  {
        TObject *obj {key->ReadObj()};
        if (obj->IsA()->InheritsFrom(TGraphErrors::Class())) {
            mygraph = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound TTree in input file [" << mygraph->GetName() << "] --> [" << file << "]";
    
    return mygraph;
}

inline std::tuple<std::shared_ptr<TGraph>, std::shared_ptr<TH1D>> get_systematics_from_eff(
    const char* efficiency_ratios,
    const char* energy_config_file,
    std::shared_ptr<TGraphErrors> flux) {
        
        TFile input_eff_ratios(efficiency_ratios, "READ");
        if (!input_eff_ratios.IsOpen()) {
            std::cerr << "\n\nError opening input file [" << efficiency_ratios << "]\n\n";
            exit(100);
        }

        auto f_ratio_trigger_eff_het_over_let_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_trigger_eff_het_over_let_bdt"));
        auto f_ratio_maxrms_eff_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_maxrms_eff_bdt"));
        auto f_ratio_nbarlayer13_eff_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_nbarlayer13_eff_bdt"));
        auto f_ratio_track_selection_eff_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_track_selection_eff_bdt"));
        auto f_ratio_psd_stk_match_eff_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_psd_stk_match_eff_bdt"));
        auto f_ratio_psd_charge_eff_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_psd_charge_eff_bdt"));
        auto f_ratio_stk_charge_eff_bdt = static_cast<TF1*>(input_eff_ratios.Get("f_ratio_stk_charge_eff_bdt"));

        auto get_sys = [&](const double energy_gev) -> double {
            double sys_trigger          {pow(1 - f_ratio_trigger_eff_het_over_let_bdt->Eval(energy_gev) , 2)};
            double sys_maxrms           {pow(1 - f_ratio_maxrms_eff_bdt->Eval(energy_gev), 2)};
            double sys_nbarlayer13      {pow(1 - f_ratio_nbarlayer13_eff_bdt->Eval(energy_gev), 2)};
            double sys_track_selection  {pow(1 - f_ratio_track_selection_eff_bdt->Eval(energy_gev), 2)};
            double sys_psd_stk_match    {pow(1 - f_ratio_psd_stk_match_eff_bdt->Eval(energy_gev), 2)};
            double sys_psd_charge       {pow(1 - f_ratio_psd_charge_eff_bdt->Eval(energy_gev) , 2)};
            double sys_stk_charge       {pow(1 - f_ratio_stk_charge_eff_bdt->Eval(energy_gev), 2)};
            return sqrt(sys_trigger + sys_maxrms + sys_nbarlayer13 + sys_track_selection + sys_psd_stk_match + sys_psd_charge + sys_stk_charge);
        };

        std::vector<double> energy (flux->GetN(), 0);
        std::vector<double> syst (flux->GetN(), 0);

        for (int pidx=0; pidx<flux->GetN(); ++pidx) {
            energy[pidx] = flux->GetPointX(pidx);
            syst[pidx] = get_sys(energy[pidx])*100;
        }

        std::shared_ptr<TGraph> syst_graph = std::make_shared<TGraph>(static_cast<int>(energy.size()), energy.data(), syst.data());
        syst_graph->SetName("syst_graph");
        syst_graph->GetXaxis()->SetTitle("Energy [GeV]");
        syst_graph->GetYaxis()->SetTitle("Systematic [%]");

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        std::shared_ptr<TH1D> histo_err = std::make_shared<TH1D>("efficiency_syst_err", "efficiency_syst_err", energy_nbins, energy_binning.data());

        for (int pidx=0; pidx<syst_graph->GetN(); ++pidx) {
            int bin = histo_err->FindBin(syst_graph->GetPointX(pidx));
            histo_err->SetBinContent(bin, 0);
            histo_err->SetBinError(bin, syst_graph->GetPointY(pidx));
        }

        
        return std::make_tuple(syst_graph, histo_err);
    }

void buildEfficiencySyst(
    const char* eff_ratio_file, 
    const char* flux_file,
    const char* energy_config_file,
    const char* output_file,
    const bool verbose = true) {

        auto flux = extractGraphFromFile(flux_file, verbose);

        std::shared_ptr<TGraph> syst_graph;
        std::shared_ptr<TH1D> histo_err;
        std::tie(syst_graph, histo_err) = get_systematics_from_eff(eff_ratio_file, energy_config_file, flux);

        std::unique_ptr<TFile> file = std::make_unique<TFile>(output_file, "RECREATE");
        if (file->IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        syst_graph->Write();
        histo_err->Write();

        file->Close();
    }