#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

// Evaluate the background contamination using fitting technique with a pol1 in background linear descent
const std::vector<std::tuple<double, double>> estimate_background(
    std::string background_linear_fit_file,
    const std::vector<std::tuple<double, unsigned int>> bdt_cuts,
    const unsigned int energy_bin,
    const bool verbose,
    const unsigned int threads) {

        if (verbose)
            std::cout << "\n\n*** background Estimation Facility***\n\n";

        // Read proton fit file
        TFile input_proton_fit_file(background_linear_fit_file.c_str(), "READ");
        if (!input_proton_fit_file.IsOpen()) {
            std::cout << "ERROR opening input proton fit file: [" << background_linear_fit_file << "]\n\n";
            exit(100);
        }

        // Extract the histogram and fit function
        std::string h_name = "energybin_" + std::to_string(energy_bin) + "/h_classifier_bin_" + std::to_string(energy_bin);
        std::string f_name = "energybin_" + std::to_string(energy_bin) + "/proton_linear_fit_" + std::to_string(energy_bin);;
        auto h_data = static_cast<TH1D*>(input_proton_fit_file.Get(h_name.c_str()));
        auto background_fit = static_cast<TF1*>(input_proton_fit_file.Get(f_name.c_str()));

        // Evaluate the background contamination
        double h_data_bin_width = h_data->GetBinWidth(1);   // The histos have uniform binning. The binning is used to rescale the integral of the TF1
        
        // Build the background tuple
        std::vector<std::tuple<double, double>> background_values (bdt_cuts.size());

        for (size_t tidx=0; tidx<background_values.size(); ++tidx) {
            double bdt_cut = std::get<0>(bdt_cuts[tidx]);
            double bacground_contamination = background_fit->Integral(bdt_cut, 1.0)/h_data_bin_width;
            background_values[tidx] = std::make_tuple(bdt_cut, bacground_contamination);
        }

        return background_values;
    }