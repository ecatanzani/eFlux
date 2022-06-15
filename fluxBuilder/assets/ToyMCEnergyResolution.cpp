#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <tuple>

#include "TKey.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
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

std::shared_ptr<TGraph> buildEnergyResolution(const char* energy_resolution_csv_file) {
    std::ifstream input_file(energy_resolution_csv_file);
    if (!input_file.is_open()) {
        std::cerr << "\nInput energ file not found [" << energy_resolution_csv_file << "]\n\n";
        exit(100);
    }
    std::string input_string(
        (std::istreambuf_iterator<char>(input_file)),
        (std::istreambuf_iterator<char>()));
    input_file.close();

    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;

    std::vector<double> energy;
    std::vector<double> energy_resolution;
    int column_counter {0};
    while (input_stream >> tmp_str) {
        if (!column_counter) {
            energy.push_back(stod(tmp_str, &sz));
            ++column_counter;
        }
        else {
            energy_resolution.push_back(stod(tmp_str, &sz));
            column_counter = 0;
        }
    }

    std::shared_ptr<TGraph> energy_resolution_graph = std::make_shared<TGraph>(energy.size(), energy.data(), energy_resolution.data());
    energy_resolution_graph->SetName("energy_resolution_graph");
    energy_resolution_graph->GetXaxis()->SetTitle("Energy [GeV]");
    energy_resolution_graph->GetYaxis()->SetTitle("Energy resolution");

    return energy_resolution_graph;
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
        std::cout << "\nFound TH1D in input file [" << myhisto->GetName() << "] --> [" << file << "]";
    
    return myhisto;
}

void ToyMCEnergyResolution(
    const char* energy_config_file,
    const char* energy_resolution_csv_file,
    const char* acceptance_file,
    const double spectral_index = -3,
    const unsigned long int nevents = 1e+10,
    const unsigned int start_seed = 7,
    const double energy_shift_value = 0.0175,
    const char* output_file = "ToyMCEnergyResolution.root",
    const bool verbose = true) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        // Extract the acceptance
        auto h_acceptance = extractHistoFromFile(acceptance_file, verbose);
        h_acceptance->Sumw2();

        // Create output file
        std::shared_ptr<TFile> output = std::make_shared<TFile>(output_file, "RECREATE");

        // Get energy resolution
        auto gr_energy_resolution = buildEnergyResolution(energy_resolution_csv_file);

        // Create histograms
        std::shared_ptr<TH1D> energy_hist = std::make_shared<TH1D>("energy_hist", "energy_hist", energy_nbins, energy_binning.data());
        std::shared_ptr<TH1D> energy_resolution_hist = std::make_shared<TH1D>("energy_resolution_hist", "energy_resolution_hist", energy_nbins, energy_binning.data());
        std::shared_ptr<TH1D> energy_resolution_p_shift_hist = std::make_shared<TH1D>("energy_resolution_p_shift_hist", "energy_resolution_p_shift_hist", energy_nbins, energy_binning.data());
        std::shared_ptr<TH1D> energy_resolution_n_shift_hist = std::make_shared<TH1D>("energy_resolution_n_shift_hist", "energy_resolution_n_shift_hist", energy_nbins, energy_binning.data());

        energy_hist->GetXaxis()->SetTitle("Energy [GeV]");
        energy_hist->GetYaxis()->SetTitle("counts GeV^{-1} sr^{-1}");

        energy_resolution_hist->GetXaxis()->SetTitle("Energy [GeV]");
        energy_resolution_hist->GetYaxis()->SetTitle("counts GeV^{-1} sr^{-1}");

        energy_resolution_p_shift_hist->GetXaxis()->SetTitle("Energy [GeV]");
        energy_resolution_p_shift_hist->GetYaxis()->SetTitle("counts GeV^{-1} sr^{-1}");
        
        energy_resolution_n_shift_hist->GetXaxis()->SetTitle("Energy [GeV]");
        energy_resolution_n_shift_hist->GetYaxis()->SetTitle("counts GeV^{-1} sr^{-1}");

        energy_hist->Sumw2();
        energy_resolution_hist->Sumw2();
        energy_resolution_p_shift_hist->Sumw2();
        energy_resolution_n_shift_hist->Sumw2();

        // Get random generator
        std::shared_ptr<TRandom3> random_generator = std::make_shared<TRandom3>(start_seed);

        /* 
            Set the max and min energy windows

            double energy_gen_min = energy_binning.front();
            double energy_gen_max = energy_binning.back();

            If the exact same binning is used for the simulation, the edge bins will contain less statistics. A wider window is used to avoid this.
        */

        double energy_gen_min = 10;
        double energy_gen_max = 2e+4;

        auto generate_energy = [&energy_gen_min, &energy_gen_max, &spectral_index] (std::shared_ptr<TRandom3> engine) -> std::tuple<double, double> {
            auto generate_weight = [&spectral_index] (const double generated_energy) -> double {
                return pow(generated_energy, spectral_index)/pow(generated_energy, -1.);
            };

            // Generate energy log-flat / E^-1
            double energy = exp(engine->Uniform(log(energy_gen_min), log(energy_gen_max)));
            double weight = generate_weight(energy);
            
            return std::make_tuple(energy, weight);
        };

        auto smear_energy = [gr_energy_resolution, &energy_shift_value] (
            const double energy, std::shared_ptr<TRandom3> engine, 
            const bool shift = false, 
            const bool positive_shift = true) -> double {
                double energy_resolution = static_cast<double>(gr_energy_resolution->Eval(energy))*energy;
                double random_energy_from_gaussian {0};
                if (!shift) {
                    random_energy_from_gaussian = engine->Gaus(energy, energy_resolution);
                }
                else {
                    if (positive_shift) {
                        random_energy_from_gaussian = engine->Gaus(energy + energy_shift_value*energy, energy_resolution);
                    }
                    else {
                        random_energy_from_gaussian = engine->Gaus(energy - energy_shift_value*energy, energy_resolution);
                    }
                }
                return random_energy_from_gaussian;
            };

        if (verbose)
            std::cout << "\n\nGenerating " << nevents << " events in [" << energy_gen_min << ", " << energy_gen_max << "] GeV\n";

        // Start loop
        unsigned long int step {static_cast<unsigned long int>(nevents/10)};
        for (unsigned long int idx=0; idx<nevents; ++idx) {
            if (!((idx+1)%step))
                std::cout << "\n\rGenerating events ... [" << (idx+1)*100./nevents << " %]";

            double genergy {0};
            double gweight {0};

            // Generate the energy
            std::tie(genergy, gweight) = generate_energy(random_generator);
            
            // Smear the energy
            auto genergy_smeared = smear_energy(genergy, random_generator);
            // Smear the energy with positive shift
            auto genergy_smeared_p_shift = smear_energy(genergy, random_generator, true, true);
            // Smear the energy with negative shift
            auto genergy_smeared_n_shift = smear_energy(genergy, random_generator, true, false);
            
            // Fill the histos
            energy_hist->Fill(genergy, gweight);
            energy_resolution_hist->Fill(genergy_smeared, gweight);
            energy_resolution_p_shift_hist->Fill(genergy_smeared_p_shift, gweight);
            energy_resolution_n_shift_hist->Fill(genergy_smeared_n_shift, gweight);
        }
            

        if (verbose)
            std::cout << "\n\nSimulation has been completed\n\n";

        // Divide histos by bin width
        auto divide_by_bin_width = [](std::shared_ptr<TH1D> hist) {
            for (int idx=1; idx<=hist->GetNbinsX(); ++idx) {
                hist->SetBinContent(idx, hist->GetBinContent(idx)/hist->GetBinWidth(idx));
                hist->SetBinError(idx, hist->GetBinError(idx)/hist->GetBinWidth(idx));
            }
        };

        // Divide by the acceptance
        energy_hist->Divide(h_acceptance.get());
        energy_resolution_hist->Divide(h_acceptance.get());
        energy_resolution_p_shift_hist->Divide(h_acceptance.get());
        energy_resolution_n_shift_hist->Divide(h_acceptance.get());

        // Divide by the energy bin width
        divide_by_bin_width(energy_hist);
        divide_by_bin_width(energy_resolution_hist);
        divide_by_bin_width(energy_resolution_p_shift_hist);
        divide_by_bin_width(energy_resolution_n_shift_hist);
        
        // Build the ratio histos
        auto h_ratio_no_shift = static_cast<TH1D*>(energy_hist->Clone("h_ratio_no_shift"));
        auto h_ratio_p_shift = static_cast<TH1D*>(energy_hist->Clone("h_ratio_p_shift"));
        auto h_ratio_n_shift = static_cast<TH1D*>(energy_hist->Clone("h_ratio_n_shift"));
        
        h_ratio_no_shift->GetYaxis()->SetTitle("Injected/Reconstructed Flux");
        h_ratio_p_shift->GetYaxis()->SetTitle("Injected/Reconstructed Flux");
        h_ratio_n_shift->GetYaxis()->SetTitle("Injected/Reconstructed Flux");

        h_ratio_no_shift->Sumw2();
        h_ratio_p_shift->Sumw2();
        h_ratio_n_shift->Sumw2();

        h_ratio_no_shift->Divide(energy_resolution_hist.get());
        h_ratio_p_shift->Divide(energy_resolution_p_shift_hist.get());
        h_ratio_n_shift->Divide(energy_resolution_n_shift_hist.get());

        // Create grs
        auto create_gr = [&energy_binning, &energy_nbins, &spectral_index] (TH1D* histo) -> std::shared_ptr<TGraphErrors> {
            
            auto wtsydp = [&spectral_index] (const double minene, const double maxene) -> double {
                double dene = maxene - minene;
                if (spectral_index != -1)
                    return pow(fabs((pow(maxene, spectral_index + 1) - pow(minene, spectral_index + 1)) / ((spectral_index + 1) * dene)), 1. / spectral_index);
                else
                    return dene / log(maxene / minene);
            };

            std::vector<double> energy (energy_nbins, 0);
            std::vector<double> y (energy_nbins, 0);
            std::vector<double> err_energy (energy_nbins, 0);
            std::vector<double> err_y (energy_nbins, 0);

            for (int bidx=1; bidx<=energy_nbins; ++bidx) {
                energy[bidx-1] = wtsydp(energy_binning[bidx-1], energy_binning[bidx]);
                y[bidx-1] = histo->GetBinContent(bidx);
                err_y[bidx-1] = histo->GetBinError(bidx);
            }

            std::shared_ptr<TGraphErrors> gr = std::make_shared<TGraphErrors>(energy_nbins, energy.data(), y.data(), err_energy.data(), err_y.data());
            gr->SetName((std::string("gr_") + std::string(histo->GetName())).c_str());
            gr->SetTitle((std::string("gr_") + std::string(histo->GetTitle())).c_str());
            gr->GetXaxis()->SetTitle("Energy [GeV]");
            gr->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());

            return gr;
        };

        // Build the Graphs
        auto gr_energy_hist = create_gr(energy_hist.get());
        auto gr_energy_resolution_hist = create_gr(energy_resolution_hist.get());
        
        auto gr_ratio_no_shift = create_gr(h_ratio_no_shift);
        auto gr_ratio_p_shift = create_gr(h_ratio_p_shift);
        auto gr_ratio_n_shift = create_gr(h_ratio_n_shift);

        // Write output
        gr_energy_resolution->Write();

        energy_hist->Write();
        energy_resolution_hist->Write();
        
        h_ratio_no_shift->Write();
        h_ratio_p_shift->Write();
        h_ratio_n_shift->Write();

        gr_ratio_no_shift->Write();
        gr_ratio_p_shift->Write();
        gr_ratio_n_shift->Write();

        // Build final canvas
        TCanvas c_fluxes("c_fluxes", "c_fluxes", 500, 500);

        gr_energy_hist->SetLineColor(kBlue);
        gr_energy_hist->SetLineWidth(2);
        gr_energy_hist->SetMarkerStyle(20);
        gr_energy_hist->SetMarkerColor(kBlue);

        gr_energy_resolution_hist->SetLineColor(kRed);
        gr_energy_resolution_hist->SetLineWidth(2);
        gr_energy_resolution_hist->SetMarkerStyle(20);
        gr_energy_resolution_hist->SetMarkerColor(kRed);

        gr_energy_hist->Draw();
        gr_energy_resolution_hist->Draw("same");

        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGrid(1,1);

        TCanvas c_ratio_no_shift("c_ratio_no_shift", "c_ratio_no_shift", 500, 500);
        
        gr_ratio_no_shift->SetLineColor(kBlack);
        gr_ratio_no_shift->SetLineWidth(2);
        gr_ratio_no_shift->SetMarkerStyle(20);
        gr_ratio_no_shift->SetMarkerColor(kBlack);

        gr_ratio_no_shift->Draw();
        
        gPad->SetLogx();
        gPad->SetGrid(1,1);

        TCanvas c_ratio_p_shift("c_ratio_p_shift", "c_ratio_p_shift", 500, 500);
        
        gr_ratio_p_shift->SetLineColor(kBlack);
        gr_ratio_p_shift->SetLineWidth(2);
        gr_ratio_p_shift->SetMarkerStyle(20);
        gr_ratio_p_shift->SetMarkerColor(kBlack);

        gr_ratio_p_shift->Draw();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        TCanvas c_ratio_n_shift("c_ratio_n_shift", "c_ratio_n_shift", 500, 500);
        
        gr_ratio_n_shift->SetLineColor(kBlack);
        gr_ratio_n_shift->SetLineWidth(2);
        gr_ratio_n_shift->SetMarkerStyle(20);
        gr_ratio_n_shift->SetMarkerColor(kBlack);

        gr_ratio_n_shift->Draw();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        TCanvas c_ratio("c_ratio", "c_ratio", 500, 500);
        
        gr_ratio_no_shift->SetLineColor(kMagenta);
        gr_ratio_no_shift->SetLineWidth(2);
        gr_ratio_no_shift->SetMarkerStyle(20);
        gr_ratio_no_shift->SetMarkerColor(kMagenta);

        gr_ratio_p_shift->SetLineColor(kBlue);
        gr_ratio_p_shift->SetLineWidth(2);
        gr_ratio_p_shift->SetMarkerStyle(20);
        gr_ratio_p_shift->SetMarkerColor(kBlue);

        gr_ratio_n_shift->SetLineColor(kRed);
        gr_ratio_n_shift->SetLineWidth(2);
        gr_ratio_n_shift->SetMarkerStyle(20);
        gr_ratio_n_shift->SetMarkerColor(kRed);

        gr_ratio_no_shift->Draw("ALP");
        gr_ratio_p_shift->Draw("LP");
        gr_ratio_n_shift->Draw("LP");

        gPad->SetLogx();
        gPad->SetGrid(1,1);


        // Save canvas
        c_fluxes.Write();
        c_ratio_no_shift.Write();
        c_ratio_p_shift.Write();
        c_ratio_n_shift.Write();
        c_ratio.Write();
        
        output->Close();
    }