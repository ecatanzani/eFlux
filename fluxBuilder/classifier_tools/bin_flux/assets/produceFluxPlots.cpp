#include <utility>
#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <numeric>
#include <limits>
#include <tuple>

#include "TAxis.h"
#include "TTree.h"
#include "TPDF.h"
#include "TPad.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TAttMarker.h"
#include "TPaveLabel.h"
#include "TEfficiency.h"
#include "TLegendEntry.h"
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

void produceFluxPlots(
    const char* input_flux_file, 
    const char* output_flux_bestcuts_file, 
    const char* energy_config_file,
    const int number_of_steps = 50,
    const bool logy = true,
    const bool draw_bestpoint_errors = false,
    const bool verbose = true,
    const bool justroot = false) {
    
    // Extract energy binning from config file
    auto energy_binning = parse_energy_config(energy_config_file);
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Find best point in energy bin
    auto getBestPointCoordinates = [](std::shared_ptr<TGraph> gr, int steps) -> std::tuple<double, double, double> {
        
        std::vector<std::pair<double, double>> pstack;      // Collection of tmp points
        std::vector<std::pair<double, double>> pmean;       // Collection of mean values, each one in a 2*steps points
        std::vector<std::pair<double, double>> pstd;        // Collection of std values, each one in a 2*steps points
        std::pair<double, double> reference_point;          // Tmp reference point

        // Evaluate mean on std::vector<std::pair<double, double>>
        auto vmean = [&reference_point](const std::vector<std::pair<double, double>> &v) -> double {
            std::vector<double> diff;
            for (auto& elm : v)
                diff.push_back(std::abs(elm.second - reference_point.second));
            return std::accumulate(diff.begin(), diff.end(), 0.)/static_cast<double>(diff.size());
        };

        // Evaluate std on std::vector<std::pair<double, double>>
        auto vstd = [&reference_point](const std::vector<std::pair<double, double>> &v) -> double {
            double reference_flux = reference_point.second;
            double accum {0};
            std::for_each (std::begin(v), std::end(v), [&](const std::pair<double, double> tmp_flux) {
                accum += std::pow((tmp_flux.second - reference_flux), 2);
            });
            double std_flux {sqrt(accum/static_cast<double>(v.size()))};
            //double std_flux_mean = std_flux/sqrt(v.size());
            return std_flux;
        };

        // Loop over graph points
        for (int gr_idx=0; gr_idx<gr->GetN(); ++gr_idx) {

            if(!gr_idx || gr_idx%steps)
                pstack.push_back(std::make_pair(gr->GetPointX(gr_idx), gr->GetPointY(gr_idx)));
            else
                reference_point = std::make_pair(gr->GetPointX(gr_idx), gr->GetPointY(gr_idx));
            
            if (static_cast<int>(pstack.size()) == 2*steps) {
                // Update mean vector
                pmean.push_back(std::make_pair(reference_point.first, vmean(pstack)));
                // Update std vector
                pstd.push_back(std::make_pair(reference_point.first, vstd(pstack)));
                // Clear stack of points
                pstack.clear();
                // Go back of steps
                gr_idx -= steps;
            }
        
        }

        // Use remaining points, if any
        if (pstack.size()) {
            // Update mean vector
            pmean.push_back(std::make_pair(reference_point.first, vmean(pstack)));
            // Update std vector
            pstd.push_back(std::make_pair(reference_point.first, vstd(pstack)));   
        }

        // Get best BDT point in range
        auto GetBestBDT = [&gr](std::vector<std::pair<double, double>> &pmean, std::vector<std::pair<double, double>> &pstd) -> std::tuple<double, double, double> {
            double best_bdt {0.25};
            double min_diff {std::numeric_limits<double>::max()};
            double low_interval_edge {-0.4};
            double high_interval_edge {0.4};
            double err {0};

            for (size_t idx=0; idx<pmean.size(); ++idx) {
                if (pmean[idx].first > low_interval_edge && pmean[idx].first < high_interval_edge) {
                    if (pmean[idx].second < min_diff) {
                        min_diff = pmean[idx].second;
                        best_bdt = pmean[idx].first;
                        err = pstd[idx].second;
                    }
                }
            }
            
            return std::make_tuple(best_bdt, gr->Eval(best_bdt), err);
        };
        
        return GetBestBDT(pmean, pstd);
    };

    // Extract fluxes from file 
    TFile input_file(input_flux_file, "READ");
    if (input_file.IsZombie()) {
        std::cerr << "\n\nError reading input ROOT file [" << input_flux_file << "]\n\n";
        exit(100);
    }

    std::vector<std::shared_ptr<TGraph>> flux_graphs_eff_corrected (energy_nbins);
    std::vector<std::shared_ptr<TGraph>> flux_graphs_eff_corrected_b_sub (energy_nbins);
    std::vector<std::shared_ptr<TGraph>> flux_graphs (energy_nbins);

    for (int bidx=0; bidx<energy_nbins; ++bidx) {
        std::string gr_name = "energybin_" + std::to_string(bidx+1) + "/gr_flux_E3_bdt_eff_corr_energybin_" + std::to_string(bidx+1);
        flux_graphs_eff_corrected[bidx] = std::shared_ptr<TGraph>(static_cast<TGraph*>(input_file.Get(gr_name.c_str())));
        gr_name = "energybin_" + std::to_string(bidx+1) + "/gr_flux_E3_bdt_energybin_" + std::to_string(bidx+1);
        flux_graphs[bidx] = std::shared_ptr<TGraph>(static_cast<TGraph*>(input_file.Get(gr_name.c_str())));
        gr_name = "energybin_" + std::to_string(bidx+1) + "/gr_flux_E3_bdt_eff_corr_b_subtracted_energybin_" + std::to_string(bidx+1);
        flux_graphs_eff_corrected_b_sub[bidx] = std::shared_ptr<TGraph>(static_cast<TGraph*>(input_file.Get(gr_name.c_str())));

        flux_graphs[bidx]->SetLineWidth(2);
        flux_graphs[bidx]->SetLineColor(kBlack);
        flux_graphs[bidx]->SetMarkerColor(kBlack);
        flux_graphs_eff_corrected[bidx]->SetLineWidth(2);
        flux_graphs_eff_corrected[bidx]->SetLineColor(kRed+2);
        flux_graphs_eff_corrected[bidx]->SetMarkerColor(kRed+2);
        flux_graphs_eff_corrected_b_sub[bidx]->SetLineWidth(2);
        flux_graphs_eff_corrected_b_sub[bidx]->SetLineColor(kBlue+2);
        flux_graphs_eff_corrected_b_sub[bidx]->SetMarkerColor(kBlue+2);
    }
    
    input_file.Close();

    // Vectors to store best points for each energy bin
    std::vector<double> best_p_X_f_ec (energy_nbins, -1);
    std::vector<double> best_p_Y_f_ec (energy_nbins, -1);
    std::vector<double> best_p_err_f_ec (energy_nbins, -1);
    std::vector<double> best_p_X_f_ec_b_sub (energy_nbins, -1);
    std::vector<double> best_p_Y_f_ec_b_sub (energy_nbins, -1);
    std::vector<double> best_p_err_f_ec_b_sub (energy_nbins, -1);

    // Build canvas
    std::vector<TCanvas*> canvas_collection (energy_nbins);
    TCanvas print_canvas("print_canvas", "print_canvas");
    print_canvas.SetTicks();

    TPaveLabel label(0.0, 0.95, 0.3, 1, "BDT classifier", "tlNDC");

    for (int bidx = 0; bidx < energy_nbins; ++bidx) {
        if (bidx) {
            print_canvas.cd();
            print_canvas.Clear();
        }

        flux_graphs[bidx]->SetTitle("flux (not corrected)");
        flux_graphs_eff_corrected[bidx]->SetTitle("flux (se)");
        flux_graphs_eff_corrected_b_sub[bidx]->SetTitle("flux (se + bs)");

        flux_graphs[bidx]->Draw();
        flux_graphs_eff_corrected[bidx]->Draw("same");
        flux_graphs_eff_corrected_b_sub[bidx]->Draw("same");

        auto best_point_coordinates = getBestPointCoordinates(flux_graphs_eff_corrected[bidx], number_of_steps);
        auto best_point_coordinates_b_sub = getBestPointCoordinates(flux_graphs_eff_corrected_b_sub[bidx], number_of_steps);

        best_p_X_f_ec[bidx] = std::get<0>(best_point_coordinates);
        best_p_Y_f_ec[bidx] = std::get<1>(best_point_coordinates);
        best_p_err_f_ec[bidx] = std::get<2>(best_point_coordinates);
        best_p_X_f_ec_b_sub[bidx] = std::get<0>(best_point_coordinates_b_sub);
        best_p_Y_f_ec_b_sub[bidx] = std::get<1>(best_point_coordinates_b_sub);
        best_p_err_f_ec_b_sub[bidx] = std::get<2>(best_point_coordinates_b_sub);

        /*
        TMarker marker(best_point_coordinates.first, best_point_coordinates.second, 41);
        TMarker marker_b_sub(best_point_coordinates_b_sub.first, best_point_coordinates_b_sub.second, 41);
        marker.SetMarkerColor(kGreen+2);
        marker.SetMarkerSize(2);
        marker.Draw("same");
        marker_b_sub.SetMarkerColor(kOrange+2);
        marker_b_sub.SetMarkerSize(2);
        marker_b_sub.Draw("same");
        */

        // Instead of TMarker I use TGraph with a single point. This permits to add errors along X axis
        std::vector<double> bp_x_f_ec {std::get<0>(best_point_coordinates)};
        std::vector<double> bp_y_f_ec {std::get<1>(best_point_coordinates)};
        std::vector<double> bp_ex_f_ec {0.};
        std::vector<double> bp_ey_f_ec {std::get<2>(best_point_coordinates)};

        std::vector<double> bp_x_f_ec_b_sub {std::get<0>(best_point_coordinates_b_sub)};
        std::vector<double> bp_y_f_ec_b_sub {std::get<1>(best_point_coordinates_b_sub)};
        std::vector<double> bp_ex_f_ec_b_sub {0.};
        std::vector<double> bp_ey_f_ec_b_sub {std::get<2>(best_point_coordinates_b_sub)};

        if (!draw_bestpoint_errors) {
            bp_ey_f_ec[0] = 0;
            bp_ey_f_ec_b_sub[0] = 0;
        }

        TGraphErrors gr_best_point_fec(bp_ex_f_ec.size(), bp_x_f_ec.data(), bp_y_f_ec.data(), bp_ex_f_ec.data(), bp_ey_f_ec.data());
        TGraphErrors gr_best_point_fec_b_sub(bp_ex_f_ec_b_sub.size(), bp_x_f_ec_b_sub.data(), bp_y_f_ec_b_sub.data(), bp_ex_f_ec_b_sub.data(), bp_ey_f_ec_b_sub.data());

        gr_best_point_fec.SetTitle("BDT cut (se)");
        gr_best_point_fec_b_sub.SetTitle("BDT cut (se + bs)");

        gr_best_point_fec.SetLineWidth(2);
        gr_best_point_fec.SetMarkerStyle(20);
        gr_best_point_fec.SetMarkerColor(kGreen+2);
        gr_best_point_fec.SetLineColor(kGreen+2);

        gr_best_point_fec_b_sub.SetLineWidth(2);
        gr_best_point_fec_b_sub.SetMarkerStyle(20);
        gr_best_point_fec_b_sub.SetMarkerColor(kOrange+2);
        gr_best_point_fec_b_sub.SetLineColor(kOrange+2);

        gr_best_point_fec.Draw("P");
        gr_best_point_fec_b_sub.Draw("P");

        gPad->SetGrid(1,1);
        if (!justroot) {
            if (logy) {
                gPad->SetLogy();
                gPad->Update(); 
                flux_graphs[bidx]->SetMinimum(1);
                flux_graphs[bidx]->SetMaximum(1e+5); 
                gPad->Update();
            }
            else {
                gPad->Update(); 
                flux_graphs[bidx]->SetMinimum(0);
                flux_graphs[bidx]->SetMaximum(300); 
                gPad->Update();
            }
        }
        
        std::string label_name = "Flux samples - energy bin " + std::to_string(bidx+1) + " - [" + std::to_string(energy_binning[bidx]) + ", " + std::to_string(energy_binning[bidx+1]) + "] GeV";
        label = TPaveLabel(0.0, 0.95, 0.3, 1, label_name.c_str(), "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);

        auto legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        if (!bidx)
            print_canvas.Print("flux_samples_summary.pdf(","Title:BDT classifier");
        else if (bidx==(energy_nbins-1))
            print_canvas.Print("flux_samples_summary.pdf)","Title:BDT classifier");
        else
            print_canvas.Print("flux_samples_summary.pdf","Title:BDT classifier");


        // Custom TCanvas
        TCanvas custom_canvas((std::string("canvas_")+std::to_string(bidx)).c_str(), (std::string("canvas_")+std::to_string(bidx)).c_str(), 500, 500);
        custom_canvas.cd();

        flux_graphs[bidx]->GetYaxis()->SetTitle("E^{3} #times Raw flux (GeV^{2} [m^{2} sr s]^{-1})");

        flux_graphs[bidx]->SetTitle("raw flux (not corrected)");
        flux_graphs_eff_corrected_b_sub[bidx]->SetTitle("flux (background + eff correction)");

        flux_graphs[bidx]->SetLineColor(kBlack);
        flux_graphs_eff_corrected_b_sub[bidx]->SetLineColor(kRed);

        flux_graphs[bidx]->SetLineWidth(2);
        flux_graphs_eff_corrected_b_sub[bidx]->SetLineWidth(2);

        flux_graphs[bidx]->Draw();
        flux_graphs_eff_corrected_b_sub[bidx]->Draw("same");
        
        gr_best_point_fec_b_sub.SetTitle("BDT - Best Cut Value");
        gr_best_point_fec_b_sub.SetMarkerColor(kBlue);
        gr_best_point_fec_b_sub.Draw("P");

        gPad->SetGrid(1,1);
        gPad->SetLogy();

        legend = custom_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            if (!strncmp(primitive->GetLabel(), "BDT - Best Cut Value", strlen("BDT - Best Cut Value")))
                primitive->SetOption("p");
            else
                primitive->SetOption("l");
        }

        canvas_collection[bidx] = static_cast<TCanvas*>(custom_canvas.DrawClone());
        canvas_collection[bidx]->SetName((std::string("canvas_")+std::to_string(bidx)).c_str());
    }

    // Build final TTree with best point coordinates
    TFile best_cuts(output_flux_bestcuts_file, "RECREATE");
    if (!best_cuts.IsOpen()) {
        std::cout << "Error writing output ROOT file [" << output_flux_bestcuts_file << "]\n\n";
        exit(100);
    }

    TTree best_cuts_tree("best_cuts_tree", "best_cuts_tree");
    
    double f_ec_X, f_ec_Y, f_ec_err;
    double f_ec_b_sub_X, f_ec_b_sub_Y, f_ec_b_sub_err;
    int energy_bin;

    best_cuts_tree.Branch("f_ec_X", &f_ec_X, "f_ec_X/D");
    best_cuts_tree.Branch("f_ec_Y", &f_ec_Y, "f_ec_Y/D");
    best_cuts_tree.Branch("f_ec_err", &f_ec_err, "f_ec_err/D");

    best_cuts_tree.Branch("f_ec_b_sub_X", &f_ec_b_sub_X, "f_ec_b_sub_X/D");
    best_cuts_tree.Branch("f_ec_b_sub_Y", &f_ec_b_sub_Y, "f_ec_b_sub_Y/D");
    best_cuts_tree.Branch("f_ec_b_sub_err", &f_ec_b_sub_err, "f_ec_b_sub_err/D");

    best_cuts_tree.Branch("energy_bin", &energy_bin, "energy_bin/I");

    if (verbose)
        std::cout << "\n\nSystematics connedted to the BDT best point\n\n";

    for (int bidx = 0; bidx < energy_nbins; ++bidx) {
        f_ec_X = best_p_X_f_ec[bidx];
        f_ec_Y = best_p_Y_f_ec[bidx];
        f_ec_err = best_p_err_f_ec[bidx];
        f_ec_b_sub_X = best_p_X_f_ec_b_sub[bidx];
        f_ec_b_sub_Y = best_p_Y_f_ec_b_sub[bidx];
        f_ec_b_sub_err = best_p_err_f_ec_b_sub[bidx];

        if (verbose)
            std::cout << "\nEnergy bin " << bidx+1 << ": syst: " << f_ec_b_sub_err << " -> [%]: " << 100*(f_ec_b_sub_err/f_ec_b_sub_Y);

        if (bidx==24)       // bin number 25
            f_ec_b_sub_X = 0;
        else if (bidx==27)  // bin number 28
            f_ec_b_sub_X = 0.2;

        # if 0
        // Exceptions
        if (bidx==23)       // bin number 24
            f_ec_b_sub_X = 0.05;
        else if (bidx==26)  // bin number 27
            f_ec_b_sub_X = 0.17;

        else if (bidx==30)  // bin number 31
            f_ec_b_sub_X = 0;
        else if (bidx==31)  // bin number 32
            f_ec_b_sub_X = 0.15;

        /*
        else if (bidx==32)  // bin number 33
            f_ec_b_sub_X = 0.2;
        */
        else if (bidx==33)  // bin number 34
            f_ec_b_sub_X = 0.2;
        else if (bidx==34)  // bin number 35
            f_ec_b_sub_X = 0.2;

        else if (bidx==39)  // bin number 40
            f_ec_b_sub_X = 0.15;
        else if (bidx==40)  // bin number 41
            f_ec_b_sub_X = 0.2;
        #endif

        energy_bin = bidx + 1;
    
        best_cuts_tree.Fill();
    }

    best_cuts_tree.Write();
    
    best_cuts.Close();

    TFile canvas_out_file("canvas.root", "RECREATE");
    if (!canvas_out_file.IsOpen()) {
        std::cout << "Error writing output ROOT file [canvas.root] \n\n";
        exit(100);
    }

    for (auto canvas : canvas_collection) {
        canvas->Write();
    }

    canvas_out_file.Close();

}