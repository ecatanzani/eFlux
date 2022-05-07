#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TPDF.h"
#include "TPad.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TPaveLabel.h"
#include "TEfficiency.h"
#include "TLegendEntry.h"
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

void produceFluxPlots(
    const char* input_flux_file, 
    const char* output_flux_file, 
    const char* energy_config_file,
    const bool logy,
    const bool verbose) {
    
    // Extract energy binning from config file
    auto energy_binning = parse_energy_config(energy_config_file);
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Extract fluxes from file 
    TFile input_file(input_flux_file, "READ");
    if (input_file.IsZombie()) {
        std::cerr << "\n\nError reading input ROOT file [" << input_flux_file << "]\n\n";
        exit(100);
    }

    std::vector<std::shared_ptr<TGraphErrors>> flux_graphs_eff_corrected_we (energy_nbins);
    std::vector<std::shared_ptr<TGraphErrors>> flux_graphs (energy_nbins);

    for (int bidx=0; bidx<energy_nbins; ++bidx) {
        std::string gr_name = "energybin_" + std::to_string(bidx+1) + "/gr_flux_E3_bdt_eff_corr_energybin_" + std::to_string(bidx+1);
        flux_graphs_eff_corrected_we[bidx] = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(input_file.Get(gr_name.c_str())));
        gr_name = "energybin_" + std::to_string(bidx+1) + "/gr_flux_E3_bdt_energybin_" + std::to_string(bidx+1);
        flux_graphs[bidx] = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(input_file.Get(gr_name.c_str())));

        flux_graphs[bidx]->SetLineWidth(2);
        flux_graphs[bidx]->SetLineColor(kBlack);
        flux_graphs[bidx]->SetMarkerColor(kBlack);
        flux_graphs_eff_corrected_we[bidx]->SetLineWidth(2);
        flux_graphs_eff_corrected_we[bidx]->SetLineColor(kRed+2);
        flux_graphs_eff_corrected_we[bidx]->SetMarkerColor(kRed+2);
    }
    
    input_file.Close();
    
    // Build canvas
    TCanvas print_canvas("print_canvas", "print_canvas");
    print_canvas.SetTicks();

    TPaveLabel label(0.0, 0.95, 0.3, 1, "BDT classifier", "tlNDC");

    for (int bidx = 0; bidx < energy_nbins; ++bidx) {
        if (bidx)
            print_canvas.Clear();

        flux_graphs[bidx]->Draw();
        flux_graphs_eff_corrected_we[bidx]->Draw("same");
        
        gPad->SetGrid(1,1);
        if (logy) {
            gPad->SetLogy();
            gPad->Update(); 
            flux_graphs[bidx]->SetMinimum(1);
            flux_graphs[bidx]->SetMaximum(1e+5); 
            gPad->Update();
        }
        
        std::string label_name = "Flux samples - energy bin " + std::to_string(bidx+1) + " - [" + std::to_string(energy_binning[bidx]) + ", " + std::to_string(energy_binning[bidx+1]) + "] GeV";
        label = TPaveLabel(0.0, 0.95, 0.3, 1, label_name.c_str(), "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);

        if (!bidx)
            print_canvas.Print("flux_samples_summary.pdf(","Title:BDT classifier");
        else if (bidx==(energy_nbins-1))
            print_canvas.Print("flux_samples_summary.pdf)","Title:BDT classifier");
        else
            print_canvas.Print("flux_samples_summary.pdf","Title:BDT classifier");
    }

}