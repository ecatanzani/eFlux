#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include <tuple>

#include "TF1.h"
#include "TAxis.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

// Sliding window technique using 3 (not-independent) points
void spectralAnalysis_3PNI(
    const char* input_flux_file,
    const char* output_file,
    const bool verbose = true,
    const bool remove_first_points = true) {

        // Extract flux graph
        TFile flux_file(input_flux_file, "READ");
        if (!flux_file.IsOpen()) {
            std::cout << "Error opening input ROOT file [" << input_flux_file << "]\n\n";
            exit(100);
        }

        const char* flux_file_name = "gr_flux_E3";
        TGraphErrors* flux_graph = static_cast<TGraphErrors*>(flux_file.Get(flux_file_name));
        if (remove_first_points) {
            flux_graph->RemovePoint(0);
            flux_graph->RemovePoint(1);
        }

        std::vector<double> energy;
        std::vector<double> energy_err;
        std::vector<double> spectral_index;
        std::vector<double> spectral_index_err;

        std::vector<std::shared_ptr<TGraphErrors>> spectral_index_graphs;
        std::vector<std::shared_ptr<TF1>> spectral_index_functions;

        auto wtsydp = [](const double minene, const double maxene, const double index) -> double {
            auto dene = maxene - minene;
            if (index != -1)
                return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
            else
                return dene / log(maxene / minene);
        };

        // Loop over the points to fit the spectrum
        for (int pidx=1; pidx<flux_graph->GetN()-1; ++pidx) {
            std::vector<double> energy_tmp {flux_graph->GetPointX(pidx-1), flux_graph->GetPointX(pidx), flux_graph->GetPointX(pidx+1)};
            std::vector<double> flux_tmp {flux_graph->GetPointY(pidx-1), flux_graph->GetPointY(pidx), flux_graph->GetPointY(pidx+1)};
            std::vector<double> energy_err_tmp (3, 0);
            std::vector<double> flux_err_tmp {flux_graph->GetErrorY(pidx-1), flux_graph->GetErrorY(pidx), flux_graph->GetErrorY(pidx+1)};

            std::shared_ptr<TGraphErrors> gr_tmp = std::make_shared<TGraphErrors>(static_cast<int>(energy_tmp.size()), energy_tmp.data(), flux_tmp.data(), energy_err_tmp.data(), flux_err_tmp.data());
            gr_tmp->SetName((std::string("gr_flux_E3_") + std::to_string(pidx)).c_str());
            gr_tmp->GetXaxis()->SetTitle("Energy [GeV]");
            std::shared_ptr<TF1> fitfunc_tmp = std::make_shared<TF1>("fitfunc_tmp", "pow(x,3) * pow(10, ([0] + [1]*log10(x)))", energy_tmp.front(), energy_tmp.back());
            gr_tmp->Fit(fitfunc_tmp.get(), "QIR");

            energy.push_back(wtsydp(energy_tmp.front(), energy_tmp.back(), fitfunc_tmp->GetParameter(1)));
            energy_err.push_back(0);
            spectral_index.push_back(fitfunc_tmp->GetParameter(1));
            spectral_index_err.push_back(fitfunc_tmp->GetParError(1));

            spectral_index_graphs.push_back(gr_tmp);
            spectral_index_functions.push_back(fitfunc_tmp);
        }
        
        // Plot the results
        std::unique_ptr<TGraphErrors> gr_spectral_index = std::make_unique<TGraphErrors>(energy.size(), energy.data(), spectral_index.data(), energy_err.data(), spectral_index_err.data());
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

// Sliding window technique using N independent points
void spectralAnalysisNP(
    const char* input_flux_file,
    const char* output_file,
    const int npoints = 5,
    const bool verbose = true,
    const bool remove_first_points = true) {

        // Extract flux graph
        TFile flux_file(input_flux_file, "READ");
        if (!flux_file.IsOpen()) {
            std::cout << "Error opening input ROOT file [" << input_flux_file << "]\n\n";
            exit(100);
        }

        const char* flux_file_name = "gr_flux_E3";
        TGraphErrors* flux_graph = static_cast<TGraphErrors*>(flux_file.Get(flux_file_name));
        if (remove_first_points) {
            flux_graph->RemovePoint(0);
            flux_graph->RemovePoint(1);
        }

        std::vector<double> energy;
        std::vector<double> spectral_index;
        std::vector<double> energy_err;
        std::vector<double> spectral_index_err;

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
            std::vector<double> &energy_err,
            std::vector<double> &spectral_index,
            std::vector<double> &spectral_index_err,
            std::vector<std::shared_ptr<TGraphErrors>> &spectral_index_graphs,
            std::vector<std::shared_ptr<TF1>> &spectral_index_functions) {
                
                auto wtsydp = [](const double minene, const double maxene, const double index) -> double {
                    auto dene = maxene - minene;
                    if (index != -1)
                        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
                    else
                        return dene / log(maxene / minene);
                };

                std::shared_ptr<TGraphErrors> gr_tmp = std::make_shared<TGraphErrors>(static_cast<int>(energy_tmp.size()), energy_tmp.data(), flux_tmp.data(), energy_err_tmp.data(), flux_err_tmp.data());
                gr_tmp->SetName((std::string("gr_flux_E3_") + std::to_string(iteration)).c_str());
                gr_tmp->GetXaxis()->SetTitle("Energy [GeV]");
                std::shared_ptr<TF1> fitfunc_tmp = std::make_shared<TF1>("fitfunc_tmp", "pow(x,3) * pow(10, ([0] + [1]*log10(x)))", energy_tmp.front(), energy_tmp.back());
                gr_tmp->Fit(fitfunc_tmp.get(), "QIR");

                //double mean_energy = std::accumulate(energy_tmp.begin(), energy_tmp.end(), 0.) / energy_tmp.size();
                double wtsydp_energy = wtsydp(energy_tmp.front(), energy_tmp.back(), fitfunc_tmp->GetParameter(1));

                //energy.push_back(mean_energy);
                energy.push_back(wtsydp_energy);
                energy_err.push_back(0);
                spectral_index.push_back(fitfunc_tmp->GetParameter(1));
                spectral_index_err.push_back(fitfunc_tmp->GetParError(1));

                spectral_index_graphs.push_back(gr_tmp);
                spectral_index_functions.push_back(fitfunc_tmp);

                std::cout << "\nEmin: " << energy_tmp.front() << " (GeV) - Emax: " << energy_tmp.back() << " (GeV) - WTSYDP: " << wtsydp_energy << " (GeV) - " << fitfunc_tmp->GetParameter(1) << " +- " << fitfunc_tmp->GetParError(1);
            };

        // Loop over the points to fit the spectrum
        for (int pidx=0; pidx<flux_graph->GetN(); ++pidx) {

            if (collected_points==npoints) {
                ++iteration;
                
                // Fit
                fit_points(energy, energy_err, spectral_index, spectral_index_err, spectral_index_graphs, spectral_index_functions);

                // Reset
                collected_points = 0;
                energy_tmp.clear();
                energy_err_tmp.clear();
                flux_tmp.clear();
                flux_err_tmp.clear();
            }
                
            energy_tmp.push_back(flux_graph->GetPointX(pidx));
            energy_err_tmp.push_back(0);
            flux_tmp.push_back(flux_graph->GetPointY(pidx));
            flux_err_tmp.push_back(flux_graph->GetErrorY(pidx));
            
            ++collected_points;
            
        }

        if (energy_tmp.size()>1) {
            ++iteration;
            fit_points(energy, energy_err, spectral_index, spectral_index_err, spectral_index_graphs, spectral_index_functions);
        }
        
        // Plot the results
        std::unique_ptr<TGraphErrors> gr_spectral_index = std::make_unique<TGraphErrors>(static_cast<int>(energy.size()), energy.data(), spectral_index.data(), energy_err.data(), spectral_index_err.data());
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

void spectralAnalysisCoupleOfPoints(
    const char* input_flux_file,
    const char* output_file) {

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
        std::vector<double> energy_err;
        std::vector<double> spectral_index_err;

        std::cout << "\nThe flux has " << flux_graph->GetN() << " points\n";

        auto wtsydp = [](const double minene, const double maxene, const double index) -> double {
            auto dene = maxene - minene;
            if (index != -1)
                return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
            else
                return dene / log(maxene / minene);
        };

        for (int pidx=1; pidx<flux_graph->GetN()-1; pidx+=2) {

            // Here I can use point pidx and pidx+1

            double tmp_min_energy {flux_graph->GetPointX(pidx)};
            double tmp_max_energy {flux_graph->GetPointX(pidx+1)};

            double tmp_flux_first_point {flux_graph->GetPointY(pidx)};
            double tmp_flux_second_point {flux_graph->GetPointY(pidx+1)};

            double tmp_flux_err_first_point {flux_graph->GetErrorY(pidx)};
            double tmp_flux_err_second_point {flux_graph->GetErrorY(pidx+1)};

            // Calculate the spectral index
            double tmp_spectral_index = (log10(tmp_flux_second_point/pow(tmp_max_energy, 3)) - log10(tmp_flux_first_point/pow(tmp_min_energy, 3)))/(log10(tmp_max_energy) - log10(tmp_min_energy));

            // Calculate the error on the spectral index
            double tmp_spectral_index_err = sqrt(
                
                pow( (1./(log10(tmp_max_energy)-log10(tmp_min_energy))) * (1./(tmp_flux_first_point*log(10))) * tmp_flux_err_first_point, 2) + 
                pow( (1./(log10(tmp_max_energy)-log10(tmp_min_energy))) * (1./(tmp_flux_second_point*log(10))) * tmp_flux_err_second_point, 2)
            );

            // Calculate the WTSYDP energy
            double wtsydp_energy = wtsydp(tmp_min_energy, tmp_max_energy, tmp_spectral_index);

            //std::cout << "\nEnergy: " << wtsydp_energy << "\tSpectral index: " << spectral_index << "\tpoints: " << pidx << "-" << pidx+1;
            std::cout << "\nFlux: " << tmp_flux_first_point << "-" << tmp_flux_second_point << "\tWTSYDP energy: " << wtsydp_energy << "\tSpectral index: " << tmp_spectral_index << "\t error: " << tmp_spectral_index_err << "\tpoints: " << pidx << "-" << pidx+1;

            // Fill the vectors
            energy.push_back(wtsydp_energy);
            energy_err.push_back(0);
            spectral_index.push_back(tmp_spectral_index);
            spectral_index_err.push_back(tmp_spectral_index_err);

        }

        // Plot the results
        std::unique_ptr<TGraphErrors> gr_spectral_index = std::make_unique<TGraphErrors>(static_cast<int>(energy.size()), energy.data(), spectral_index.data(), energy_err.data(), spectral_index_err.data());
        gr_spectral_index->SetName("gr_spectral_index");
        gr_spectral_index->GetXaxis()->SetTitle("Energy [GeV]");
        gr_spectral_index->GetYaxis()->SetTitle("Spectral index");

        TFile ofile(output_file, "RECREATE");
        if (!ofile.IsOpen()) {
            std::cout << "Error writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        gr_spectral_index->Write();

        ofile.Close();

    }

// GammaDiff

// Working broken power-law fit
double fit_func(double *x, double *par) {
    double energy               {x[0]};
    double function_value       {0};
    const double break_energy   {2800};

    if (energy<break_energy)
        function_value = pow(energy, 3) * pow(10, par[0] + par[1]*log10(energy));
    else
        function_value = pow(energy, 3) * pow(10, + (par[0]+(par[1]-par[2])*log10(break_energy)) + par[2]*log10(energy));
    return function_value;
}

double fit_func_fixed_EB(double *x, double *par) {
    double energy               {x[0]};
    double function_value       {0};
    const double break_energy   {2505.07};

    if (energy<break_energy)
        function_value = pow(energy, 3) * pow(10, par[0] + par[1]*log10(energy));
    else
        function_value = pow(energy, 3) * pow(10, + (par[0]+(par[1]-par[2])*log10(break_energy)) + par[2]*log10(energy));
    return function_value;
}

double fit_func_floating_break(double *x, double *par) {
    double energy               {x[0]};
    double function_value       {0};

    if (energy<par[3])
        function_value = pow(energy, 3) * pow(10, par[0] + par[1]*log10(energy));
    else
        function_value = pow(energy, 3) * pow(10, + (par[0]+(par[1]-par[2])*log10(par[3])) + par[2]*log10(energy));
    return function_value;
}

double fit_func_floating_break_gammadiff(double *x, double *par) {
    double energy               {x[0]};
    double function_value       {0};

    if (energy<par[3])
        function_value = pow(energy, 3) * pow(10, par[0] + par[1]*log10(energy));
    else 
        //par[2] = gamma diff
        function_value = pow(energy, 3) * pow(10, + (par[0]+par[2]*log10(par[3])) + (par[1]-par[2])*log10(energy));
    
    return function_value;
}

double fit_func_floating_break_gammadiff_fixed_EB(double *x, double *par) {
    double energy               {x[0]};
    double function_value       {0};
    const double break_energy   {2505.07};

    if (energy<break_energy)
        function_value = pow(energy, 3) * pow(10, par[0] + par[1]*log10(energy));
    else 
        //par[2] = gamma diff
        function_value = pow(energy, 3) * pow(10, + (par[0]+par[2]*log10(break_energy)) + (par[1]-par[2])*log10(energy));
    
    return function_value;
}

/*
double fit_func(double *x, double *par) {
    double energy               {x[0]};
    double function_value       {0};
    const double break_energy   {2800};

    if (energy<break_energy)
        function_value = pow(energy, 3) * pow(10, par[0] + par[1]*log10(energy));
    else {
        par[3] = par[1]-par[2];
        function_value = pow(energy, 3) * pow(10, + (par[0]+(par[1]+par[3])*log10(break_energy)) + (par[1]-par[3])*log10(energy));
    }
    return function_value;
}
*/

void GammaDiff(
        const char* input_flux_file,
        const char* output_flux_file = "gamma_diff.root",
        const double fit_start_energy = 900) {
            
            // Extract flux graph
            TFile flux_file(input_flux_file, "READ");
            if (!flux_file.IsOpen()) {
                std::cout << "Error opening input ROOT file [" << input_flux_file << "]\n\n";
                exit(100);
            }

            const char* flux_file_name = "gr_flux_E3";
            TGraphErrors* flux_graph = static_cast<TGraphErrors*>(flux_file.Get(flux_file_name));

            #if 0
            std::vector<double> pre_break_energy_tmp;
            std::vector<double> pre_break_energy_err_tmp;
            std::vector<double> pre_break_flux_tmp;
            std::vector<double> pre_break_flux_err_tmp;

            std::vector<double> post_break_energy_tmp;
            std::vector<double> post_break_energy_err_tmp;
            std::vector<double> post_break_flux_tmp;
            std::vector<double> post_break_flux_err_tmp;

            for (int pidx=0; pidx<flux_graph->GetN(); ++pidx) {
                
                auto energy_tmp = flux_graph->GetPointX(pidx);

                if (energy_tmp>=pre_break_power_law_start && energy_tmp<pre_break_power_law_end) {
                    pre_break_energy_tmp.push_back(energy_tmp);
                    pre_break_energy_err_tmp.push_back(0);
                    pre_break_flux_tmp.push_back(flux_graph->GetPointY(pidx));
                    pre_break_flux_err_tmp.push_back(flux_graph->GetErrorY(pidx));
                }
                else if (energy_tmp>=post_break_power_law_start && energy_tmp<post_break_power_law_end) {
                    post_break_energy_tmp.push_back(energy_tmp);
                    post_break_energy_err_tmp.push_back(0);
                    post_break_flux_tmp.push_back(flux_graph->GetPointY(pidx));
                    post_break_flux_err_tmp.push_back(flux_graph->GetErrorY(pidx));
                }
            }

            std::vector<std::shared_ptr<TGraphErrors>> spectral_index_graphs;
            std::vector<std::shared_ptr<TF1>> spectral_index_functions;

            auto fit_points = [&spectral_index_graphs, &spectral_index_functions](
                    std::vector<double> &energy, 
                    std::vector<double> &energy_err,
                    std::vector<double> &flux,
                    std::vector<double> &flux_err) -> tuple<double, double> {
                
                        std::shared_ptr<TGraphErrors> gr_tmp = std::make_shared<TGraphErrors>(static_cast<int>(energy.size()), energy.data(), flux.data(), energy_err.data(), flux_err.data());
                        std::string gr_name {std::string("gr_flux_E3_") + std::to_string(energy.front()) + std::string("_") + std::to_string(energy.back())};
                        gr_tmp->SetName(gr_name.c_str());
                        gr_tmp->GetXaxis()->SetTitle("Energy [GeV]");
                        std::shared_ptr<TF1> fitfunc_tmp = std::make_shared<TF1>("fitfunc_tmp", "pow(x,3) * pow(10, ([0] + [1]*log10(x)))", energy.front(), energy.back());
                        gr_tmp->Fit(fitfunc_tmp.get(), "QIR");

                        spectral_index_graphs.push_back(gr_tmp);
                        spectral_index_functions.push_back(fitfunc_tmp);

                        return std::make_tuple(fitfunc_tmp->GetParameter(1), fitfunc_tmp->GetParError(1));
                    };

            
            double pre_break_gamma          {0};
            double pre_break_gamma_err      {0};
            double post_break_gamma         {0};
            double post_break_gamma_err     {0};

            std::tie(pre_break_gamma, pre_break_gamma_err) = fit_points(pre_break_energy_tmp, pre_break_energy_err_tmp, pre_break_flux_tmp, pre_break_flux_err_tmp);
            std::tie(post_break_gamma, post_break_gamma_err) = fit_points(post_break_energy_tmp, post_break_energy_err_tmp, post_break_flux_tmp, post_break_flux_err_tmp);

            std::cout << "\n\nAnalysis results...\n";
            std::cout << "\nPre-Break spectral index: " << pre_break_gamma << " +/- " << pre_break_gamma_err;
            std::cout << "\nPost-Break spectral index: " << post_break_gamma << " +/- " << post_break_gamma_err;
            std::cout << "\n\nDeltaGamma: " << (post_break_gamma - pre_break_gamma) << " +/- " << sqrt(pow(post_break_gamma_err, 2) + pow(pre_break_gamma_err, 2));

            std::cout << "\n\nAnalysis completed\n\n";

            TFile output_file(output_flux_file, "RECREATE");
            if (!output_file.IsOpen()) {
                std::cout << "Error writing output ROOT file [" << output_flux_file << "]\n\n";
                exit(100);
            }

            for (size_t idx=0; idx<spectral_index_graphs.size(); ++idx) {

                if (!idx) {
                    output_file.mkdir("pre-break");
                    output_file.cd("pre-break");
                }
                else {
                    output_file.mkdir("post-break");
                    output_file.cd("post-break");
                }

                spectral_index_graphs[idx]->Write();
                spectral_index_functions[idx]->Write();
            }

            output_file.cd();
            
            TCanvas he_fit("he_fit", "he_fit", 500, 500);
            
            flux_graph->SetMarkerStyle(20);
            flux_graph->Draw("AP");

            gPad->Update(); 
            flux_graph->SetMinimum(0);
            flux_graph->SetMaximum(300); 
            gPad->Update();
            
            spectral_index_functions[0]->SetLineWidth(5);
            spectral_index_functions[0]->SetLineStyle(3);
            spectral_index_functions[0]->SetLineColor(kRed);

            spectral_index_functions[1]->SetLineWidth(5);
            spectral_index_functions[1]->SetLineStyle(3);
            spectral_index_functions[1]->SetLineColor(kBlue);

            spectral_index_functions[0]->Draw("same");
            spectral_index_functions[1]->Draw("same"); 

            gPad->SetGrid(1,1);
            gPad->SetLogx();
            gPad->SetLogz();
            gStyle->SetOptStat(0);

            he_fit.Write();

            output_file.Close();

            #endif

            // Define the fitting function

            //std::shared_ptr<TF1> broken_pl = std::make_shared<TF1>("broken_pl", fit_func, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 3); // Working broken powe-law
            std::shared_ptr<TF1> broken_pl = std::make_shared<TF1>("broken_pl", fit_func, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 3);
            
            flux_graph->Fit(broken_pl.get(), "QIR");

            std::shared_ptr<TF1> broken_pl_floating_break = std::make_shared<TF1>("broken_pl_floating_break", fit_func_floating_break, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 4);
            broken_pl_floating_break->SetParameter(0, broken_pl->GetParameter(0));
            broken_pl_floating_break->SetParameter(1, broken_pl->GetParameter(1));
            broken_pl_floating_break->SetParameter(2, broken_pl->GetParameter(2));
            broken_pl_floating_break->SetParameter(3, 2800);
            broken_pl_floating_break->SetParLimits(3, 2000, 4000);

            flux_graph->Fit(broken_pl_floating_break.get(), "QIR");

            std::cout << "\n\nAnalysis results...\n";
            std::cout << "\nEnergy break has been fixed at: " << broken_pl_floating_break->GetParameter(3) << " +/- " << broken_pl_floating_break->GetParError(3) << std::endl;
            std::cout << "\nPre-Break spectral index: " << broken_pl_floating_break->GetParameter(1) << " +/- " << broken_pl_floating_break->GetParError(1);
            std::cout << "\nPost-Break spectral index: " <<  broken_pl_floating_break->GetParameter(2) << " +/- " << broken_pl_floating_break->GetParError(2);

            double deltagamma = broken_pl_floating_break->GetParameter(1)-broken_pl_floating_break->GetParameter(2);
            double deltagamma_err = sqrt(pow(broken_pl_floating_break->GetParError(1), 2) + pow(broken_pl_floating_break->GetParError(2), 2));
            std::cout << "\n\nDeltaGamma: " <<  deltagamma << " +/- " << deltagamma_err;

            std::cout << "\n\nAnalysis completed\n\n";

            
            std::shared_ptr<TF1> broken_pl_floating_break_gdiff = std::make_shared<TF1>("broken_pl_floating_break_gdiff", fit_func_floating_break_gammadiff, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 4);

            broken_pl_floating_break_gdiff->SetParameter(0, broken_pl_floating_break->GetParameter(0));
            broken_pl_floating_break_gdiff->SetParameter(1, broken_pl_floating_break->GetParameter(1));
            broken_pl_floating_break_gdiff->SetParameter(3, broken_pl_floating_break->GetParameter(3));
            
            flux_graph->Fit(broken_pl_floating_break_gdiff.get(), "QIR");
            
            std::cout << "\n\nAnalysis results...\n";
            std::cout << "\nEnergy break has been fixed at: " << broken_pl_floating_break_gdiff->GetParameter(3) << " +/- " << broken_pl_floating_break_gdiff->GetParError(3) << std::endl;
            std::cout << "\nPre-Break spectral index: " << broken_pl_floating_break_gdiff->GetParameter(1) << " +/- " << broken_pl_floating_break_gdiff->GetParError(1);
            //std::cout << "\nPost-Break spectral index: " <<  broken_pl_floating_break_gdiff->GetParameter(2) << " +/- " << broken_pl_floating_break_gdiff->GetParError(2);
            std::cout << "\n\nDeltaGamma: " << broken_pl_floating_break_gdiff->GetParameter(2) << " +/- " << broken_pl_floating_break_gdiff->GetParError(2);

            std::cout << "\n\nAnalysis completed\n\n";

            TFile output_file(output_flux_file, "RECREATE");
            if (!output_file.IsOpen()) {
                std::cout << "Error writing output ROOT file [" << output_flux_file << "]\n\n";
                exit(100);
            }

            TCanvas he_fit("he_fit", "he_fit", 500, 500);
            
            flux_graph->SetMarkerStyle(20);
            flux_graph->GetYaxis()->SetTitle("E^{3} #times #Phi (GeV [m^{2} sr s]^{-1})");
            flux_graph->Draw("AP");

            gPad->Update(); 
            flux_graph->SetMinimum(10);
            flux_graph->SetMaximum(1000); 
            gPad->Update();
            
            broken_pl_floating_break->Draw("same");

            //broken_pl_floating_break_gdiff->SetLineColor(kBlue);
            //broken_pl_floating_break_gdiff->Draw("same");

            gPad->SetGrid(1,1);
            gPad->SetLogx();
            gPad->SetLogy();
            gStyle->SetOptStat(0);

            he_fit.Write();

            flux_graph->Write();
            broken_pl->Write();

            output_file.Close();            

        }


void GammaDiffNoEB(
        const char* input_flux_file,
        const char* output_flux_file = "gamma_diff_nofixedeb.root",
        const double break_energy = 2505.07,
        const double fit_start_energy = 900) {
            
            // Extract flux graph
            TFile flux_file(input_flux_file, "READ");
            if (!flux_file.IsOpen()) {
                std::cout << "Error opening input ROOT file [" << input_flux_file << "]\n\n";
                exit(100);
            }

            const char* flux_file_name = "gr_flux_E3";
            TGraphErrors* flux_graph = static_cast<TGraphErrors*>(flux_file.Get(flux_file_name));

            // Define the fitting function

            //std::shared_ptr<TF1> broken_pl = std::make_shared<TF1>("broken_pl", fit_func, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 3); // Working broken powe-law
            std::shared_ptr<TF1> broken_pl = std::make_shared<TF1>("broken_pl", fit_func_fixed_EB, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 3);
            
            flux_graph->Fit(broken_pl.get(), "QIR");
            
            std::cout << "\n\nAnalysis results...\n";
            std::cout << "\nPre-Break spectral index: " << broken_pl->GetParameter(1) << " +/- " << broken_pl->GetParError(1);
            std::cout << "\nPost-Break spectral index: " <<  broken_pl->GetParameter(2) << " +/- " << broken_pl->GetParError(2);

            double pre_break_si = broken_pl->GetParameter(1);
            double pre_break_si_err = broken_pl->GetParError(1);

            double post_break_si = broken_pl->GetParameter(2);
            double post_break_si_err = broken_pl->GetParError(2);

            std::shared_ptr<TF1> broken_pl_gdiff = std::make_shared<TF1>("broken_pl_gdiff", fit_func_floating_break_gammadiff_fixed_EB, fit_start_energy, flux_graph->GetPointX(flux_graph->GetN()-1), 3);

            broken_pl_gdiff->SetParameter(0, broken_pl->GetParameter(0));
            broken_pl_gdiff->SetParameter(1, broken_pl->GetParameter(1));
            
            flux_graph->Fit(broken_pl_gdiff.get(), "QIR");
            
            std::cout << "\n\nAnalysis results...\n";
            std::cout << "\nPre-Break spectral index: " << broken_pl_gdiff->GetParameter(1) << " +/- " << broken_pl_gdiff->GetParError(1);
            //std::cout << "\nPost-Break spectral index: " <<  broken_pl_floating_break_gdiff->GetParameter(2) << " +/- " << broken_pl_floating_break_gdiff->GetParError(2);
            std::cout << "\n\nDeltaGamma: " << broken_pl_gdiff->GetParameter(2) << " +/- " << broken_pl_gdiff->GetParError(2);

            std::cout << "\n\nAnalysis completed\n\n";

            std::vector<double> preb_energy;
            std::vector<double> preb_energy_err;
            std::vector<double> preb_si;
            std::vector<double> preb_si_err;

            std::vector<double> postb_energy;
            std::vector<double> postb_energy_err;
            std::vector<double> postb_si;
            std::vector<double> postb_si_err;

            for (int idx=0; idx<flux_graph->GetN(); idx++) {
                if (flux_graph->GetPointX(idx) >= fit_start_energy && flux_graph->GetPointX(idx) < break_energy) {
                    preb_energy.push_back(flux_graph->GetPointX(idx));
                    preb_energy_err.push_back(0);
                    preb_si.push_back(pre_break_si);
                    preb_si_err.push_back(post_break_si_err);
                }
                else if (flux_graph->GetPointX(idx) >= break_energy) {
                    if (!postb_energy.size()) {
                        postb_energy.push_back(preb_energy.back());
                        postb_energy_err.push_back(0);
                        postb_si.push_back(post_break_si);
                        postb_si_err.push_back(post_break_si_err);
                    }
                    postb_energy.push_back(flux_graph->GetPointX(idx));
                    postb_energy_err.push_back(0);
                    postb_si.push_back(post_break_si);
                    postb_si_err.push_back(post_break_si_err);
                }
            }

            std::shared_ptr<TGraphErrors> pre_break_si_shape = std::make_shared<TGraphErrors>(static_cast<int>(preb_energy.size()), preb_energy.data(), preb_si.data(), preb_energy_err.data(), preb_si_err.data());
            pre_break_si_shape->SetName("pre_break_si_shape");
            pre_break_si_shape->SetFillColor(kGreen);
            pre_break_si_shape->SetFillStyle(3005);

            std::shared_ptr<TGraphErrors> post_break_si_shape = std::make_shared<TGraphErrors>(static_cast<int>(postb_energy.size()), postb_energy.data(), postb_si.data(), postb_energy_err.data(), postb_si_err.data());
            post_break_si_shape->SetName("post_break_si_shape");
            post_break_si_shape->SetFillColor(kOrange);
            post_break_si_shape->SetFillStyle(3005);

            std::shared_ptr<TMultiGraph> multigr = std::make_shared<TMultiGraph>();
            multigr->SetName("multigr");
            multigr->Add(pre_break_si_shape.get());
            multigr->Add(post_break_si_shape.get());

            std::cout << "\n\nPre-Break SI shape: " << preb_energy.front() << " - " << preb_energy.back() << "\n";
            std::cout << "\n\nPost-Break SI shape: " << postb_energy.front() << " - " << postb_energy.back() << "\n";

            // Create deparate TF1s

            std::shared_ptr<TF1> pre_break_fit_func = std::make_shared<TF1>("pre_break_fit_func", "pow(x, 3) * pow(10, [0] + [1]*log10(x))", fit_start_energy, break_energy);
            pre_break_fit_func->SetParameter(0, broken_pl_gdiff->GetParameter(0));
            pre_break_fit_func->SetParameter(1, broken_pl_gdiff->GetParameter(1));

            pre_break_fit_func->SetLineColor(kGreen-3);
            pre_break_fit_func->SetLineWidth(3);

            TLine pre_break_line(900, pre_break_si, break_energy, pre_break_si);
            pre_break_line.SetLineWidth(3);
            pre_break_line.SetLineColor(kGreen+3);

            std::shared_ptr<TF1> post_break_fit_func = std::make_shared<TF1>("post_break_fit_func", "pow(x, 3) * pow(10, + ([0] + [2]*log10(2505.07)) + ([1]-[2])*log10(x))", break_energy, flux_graph->GetPointX(flux_graph->GetN()-1));
            post_break_fit_func->SetParameter(0, broken_pl_gdiff->GetParameter(0));
            post_break_fit_func->SetParameter(1, broken_pl_gdiff->GetParameter(1));
            post_break_fit_func->SetParameter(2, broken_pl_gdiff->GetParameter(2));

            post_break_fit_func->SetLineColor(kOrange-3);
            post_break_fit_func->SetLineWidth(3);

            TLine post_break_line (break_energy, post_break_si, flux_graph->GetPointX(flux_graph->GetN()-1), post_break_si);
            pre_break_line.SetLineWidth(3);
            pre_break_line.SetLineColor(kOrange+3);

            TFile output_file(output_flux_file, "RECREATE");
            if (!output_file.IsOpen()) {
                std::cout << "Error writing output ROOT file [" << output_flux_file << "]\n\n";
                exit(100);
            }

            TCanvas he_fit("he_fit", "he_fit", 500, 500);
            
            flux_graph->SetMarkerStyle(20);
            flux_graph->GetYaxis()->SetTitle("E^{3} #times #Phi (GeV [m^{2} sr s]^{-1})");
            flux_graph->Draw("AP");

            gPad->Update(); 
            flux_graph->SetMinimum(10);
            flux_graph->SetMaximum(1000); 
            gPad->Update();
            
            broken_pl_gdiff->Draw("same");

            //broken_pl_floating_break_gdiff->SetLineColor(kBlue);
            //broken_pl_floating_break_gdiff->Draw("same");

            gPad->SetGrid(1,1);
            gPad->SetLogx();
            gPad->SetLogy();
            gStyle->SetOptStat(0);

            he_fit.Write();

            flux_graph->Write();
            broken_pl->Write();

            pre_break_si_shape->Write();
            post_break_si_shape->Write();

            multigr->Write();

            pre_break_fit_func->Write();
            post_break_fit_func->Write();

            pre_break_line.Write();
            post_break_line.Write();

            output_file.Close();            

        }