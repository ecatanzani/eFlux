#include <iostream>
#include <memory>
#include <vector>
#include<numeric>

#include "TAxis.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"

// Sliding window technique using 3 points
void spectralAnalysis(
    const char* input_flux_file,
    const char* output_file,
    const bool verbose = true) {

        // Extract flux graph
        TFile flux_file(input_flux_file, "READ");
        if (!flux_file.IsOpen()) {
            std::cout << "Error opening input ROOT file [" << input_flux_file << "]\n\n";
            exit(100);
        }

        const char* flux_file_name = "gr_flux_E3";
        TGraphErrors* flux_graph = static_cast<TGraphErrors*>(flux_file.Get(flux_file_name));

        std::vector<double> energy;
        std::vector<double> spectral_index;
        std::vector<std::shared_ptr<TGraphErrors>> spectral_index_graphs;
        std::vector<std::shared_ptr<TF1>> spectral_index_functions;

        // Loop over the points to fit the spectrum
        for (int pidx=1; pidx<flux_graph->GetN()-1; ++pidx) {
            std::vector<double> energy_tmp {flux_graph->GetPointX(pidx-1), flux_graph->GetPointX(pidx), flux_graph->GetPointX(pidx+1)};
            std::vector<double> flux_tmp {flux_graph->GetPointY(pidx-1), flux_graph->GetPointY(pidx), flux_graph->GetPointY(pidx+1)};
            std::vector<double> energy_err_tmp (3, 0);
            std::vector<double> flux_err_tmp {flux_graph->GetErrorY(pidx-1), flux_graph->GetErrorY(pidx), flux_graph->GetErrorY(pidx+1)};

            std::shared_ptr<TGraphErrors> gr_tmp = std::make_shared<TGraphErrors>(3, energy_tmp.data(), flux_tmp.data(), energy_err_tmp.data(), flux_err_tmp.data());
            gr_tmp->SetName((std::string("gr_flux_E3_") + std::to_string(pidx)).c_str());
            gr_tmp->GetXaxis()->SetTitle("Energy [GeV]");
            std::shared_ptr<TF1> fitfunc_tmp = std::make_shared<TF1>("fitfunc_tmp", "pow(x,3) * pow(10, ([0] + [1]*log10(x)))", energy_tmp.front(), energy_tmp.back());
            gr_tmp->Fit(fitfunc_tmp.get(), "QIR");

            energy.push_back(flux_graph->GetPointX(pidx));
            spectral_index.push_back(fitfunc_tmp->GetParameter(1));

            spectral_index_graphs.push_back(gr_tmp);
            spectral_index_functions.push_back(fitfunc_tmp);
        }
        
        // Plot the results
        std::unique_ptr<TGraphErrors> gr_spectral_index = std::make_unique<TGraphErrors>(energy.size(), energy.data(), spectral_index.data());
        gr_spectral_index->SetName("gr_spectral_index");
        gr_spectral_index->GetXaxis()->SetTitle("Energy [GeV]");
        gr_spectral_index->GetYaxis()->SetTitle("Spectral index");

        TFile ofile(output_file, "RECREATE");
        if (!ofile.IsOpen()) {
            std::cout << "Error writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        gr_spectral_index->Write();


        for (size_t idxp=0; idxp<spectral_index_graphs.size(); ++idxp) {
            ofile.mkdir(std::to_string(idxp+1).c_str());
            ofile.cd(std::to_string(idxp+1).c_str());
            spectral_index_graphs[idxp]->Write();
            spectral_index_functions[idxp]->Write();
        }

        ofile.Close();

    }

// Sliding window technique using N points
void spectralAnalysisNP(
    const char* input_flux_file,
    const char* output_file,
    const int npoints = 5,
    const bool verbose = true) {

        // Extract flux graph
        TFile flux_file(input_flux_file, "READ");
        if (!flux_file.IsOpen()) {
            std::cout << "Error opening input ROOT file [" << input_flux_file << "]\n\n";
            exit(100);
        }

        const char* flux_file_name = "gr_flux_E3";
        TGraphErrors* flux_graph = static_cast<TGraphErrors*>(flux_file.Get(flux_file_name));

        std::vector<double> energy;
        std::vector<double> spectral_index;
        std::vector<std::shared_ptr<TGraphErrors>> spectral_index_graphs;
        std::vector<std::shared_ptr<TF1>> spectral_index_functions;

        int collected_points {0};
        int iteration {0};
        std::vector<double> energy_tmp;
        std::vector<double> energy_err_tmp;
        std::vector<double> flux_tmp;
        std::vector<double> flux_err_tmp;

        auto fit_points = [&energy_tmp, &energy_err_tmp, &flux_tmp, &flux_err_tmp, &npoints, &iteration] (
            std::vector<double> &energy, 
            std::vector<double> &spectral_index,
            std::vector<std::shared_ptr<TGraphErrors>> &spectral_index_graphs,
            std::vector<std::shared_ptr<TF1>> &spectral_index_functions) {

                std::shared_ptr<TGraphErrors> gr_tmp = std::make_shared<TGraphErrors>(npoints, energy_tmp.data(), flux_tmp.data(), energy_err_tmp.data(), flux_err_tmp.data());
                gr_tmp->SetName((std::string("gr_flux_E3_") + std::to_string(iteration)).c_str());
                gr_tmp->GetXaxis()->SetTitle("Energy [GeV]");
                std::shared_ptr<TF1> fitfunc_tmp = std::make_shared<TF1>("fitfunc_tmp", "pow(x,3) * pow(10, ([0] + [1]*log10(x)))", energy_tmp.front(), energy_tmp.back());
                gr_tmp->Fit(fitfunc_tmp.get(), "QIR");

                double mean_energy = std::accumulate(energy_tmp.begin(), energy_tmp.end(), 0.) / energy_tmp.size();

                energy.push_back(mean_energy);
                spectral_index.push_back(fitfunc_tmp->GetParameter(1));

                spectral_index_graphs.push_back(gr_tmp);
                spectral_index_functions.push_back(fitfunc_tmp);
            };

        // Loop over the points to fit the spectrum
        for (int pidx=0; pidx<flux_graph->GetN(); ++pidx) {

            if (collected_points==npoints) {
                ++iteration;
                
                // Fit
                fit_points(energy, spectral_index, spectral_index_graphs, spectral_index_functions);

                // Reset
                collected_points = 0;
                energy_tmp.clear();
                energy_err_tmp.clear();
                flux_tmp.clear();
                flux_err_tmp.clear();
            }
            else {
                
                energy_tmp.push_back(flux_graph->GetPointX(pidx));
                energy_err_tmp.push_back(0);
                flux_tmp.push_back(flux_graph->GetPointY(pidx));
                flux_err_tmp.push_back(flux_graph->GetErrorY(pidx));
                
                ++collected_points;
            }
        }

        if (energy_tmp.size()) {
            ++iteration;
            fit_points(energy, spectral_index, spectral_index_graphs, spectral_index_functions);
        }
        
        // Plot the results
        std::unique_ptr<TGraphErrors> gr_spectral_index = std::make_unique<TGraphErrors>(energy.size(), energy.data(), spectral_index.data());
        gr_spectral_index->SetName("gr_spectral_index");
        gr_spectral_index->GetXaxis()->SetTitle("Energy [GeV]");
        gr_spectral_index->GetYaxis()->SetTitle("Spectral index");

        TFile ofile(output_file, "RECREATE");
        if (!ofile.IsOpen()) {
            std::cout << "Error writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        gr_spectral_index->Write();


        for (size_t idxp=0; idxp<spectral_index_graphs.size(); ++idxp) {
            ofile.mkdir(std::to_string(idxp+1).c_str());
            ofile.cd(std::to_string(idxp+1).c_str());
            spectral_index_graphs[idxp]->Write();
            spectral_index_functions[idxp]->Write();
        }

        ofile.Close();

    }