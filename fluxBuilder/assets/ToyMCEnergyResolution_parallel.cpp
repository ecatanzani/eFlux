#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <tuple>
#include <thread>
#include <pthread.h>

#include "TH1D.h"
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

struct material {
    unsigned int idx                                {0};
    unsigned int step                               {0};
    unsigned int nevents                            {0};
    unsigned int nevents_per_thread                 {0};
    double energy_gen_min                           {0};
    double energy_gen_max                           {0};
    double spectral_index                           {0};
    std::shared_ptr<TGraph> gr_energy_resolution    {nullptr};
    std::shared_ptr<TRandom3> random_generator      {nullptr};
    std::shared_ptr<TH1D> energy_hist               {nullptr};
    std::shared_ptr<TH1D> energy_resolution_hist    {nullptr};
};

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

void* buildEvent(void *th_args) {
    
    auto simulation_material = static_cast<material*>(th_args);

    std::cout << "\nSpowning thread with " << simulation_material->nevents_per_thread << " events ..." << std::endl;

    double genergy {0};
    double gweight {0};
    double energy_resolution {0};
    double genergy_smeared {0};

    for (unsigned int idx=0; idx < simulation_material->nevents_per_thread; ++idx) {
        
        double genergy = exp(simulation_material->random_generator->Uniform(log(simulation_material->energy_gen_min), log(simulation_material->energy_gen_max)));
        double gweight = pow(genergy, simulation_material->spectral_index)/pow(genergy, -1.);
        
        double energy_resolution = static_cast<double>(simulation_material->gr_energy_resolution->Eval(genergy))*genergy;
        double genergy_smeared = simulation_material->random_generator->Gaus(genergy, energy_resolution);
        
        // Fill the histos
        simulation_material->energy_hist->Fill(genergy, gweight);
        simulation_material->energy_resolution_hist->Fill(genergy_smeared, gweight);
    }

    pthread_exit(NULL);
}

void ToyMCEnergyResolution(
    const char* energy_config_file,
    const char* energy_resolution_csv_file,
    const double spectral_index = -3,
    const unsigned int nevents = 1e+9,
    const unsigned int start_seed = 7,
    const char* output_file = "ToyMCEnergyResolution.root",
    const bool verbose = true) {

        // Initialize the args
        std::shared_ptr<material> simulation_material = std::make_shared<material>();

        simulation_material->nevents = nevents;

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        // Create output file
        std::shared_ptr<TFile> output = std::make_shared<TFile>(output_file, "RECREATE");

        // Get energy resolution
        simulation_material->gr_energy_resolution = buildEnergyResolution(energy_resolution_csv_file);

        // Create histograms
        simulation_material->energy_hist = std::make_shared<TH1D>("energy_hist", "energy_hist", energy_nbins, energy_binning.data());
        simulation_material->energy_resolution_hist = std::make_shared<TH1D>("energy_resolution_hist", "energy_resolution_hist", energy_nbins, energy_binning.data());

        simulation_material->energy_hist->GetXaxis()->SetTitle("Energy [GeV]");
        simulation_material->energy_hist->GetYaxis()->SetTitle("Counts");

        simulation_material->energy_resolution_hist->GetXaxis()->SetTitle("Energy [GeV]");
        simulation_material->energy_resolution_hist->GetYaxis()->SetTitle("Counts");

        // Get random generator
        simulation_material->random_generator = std::make_shared<TRandom3>(start_seed);

        /* 
            Set the max and min energy windows

            double energy_gen_min = energy_binning.front();
            double energy_gen_max = energy_binning.back();

            If the exact same binning is used for the simulation, the edge bins will contain less statistics. A wider window is used to avoid this.
        */

        simulation_material->energy_gen_min = 10;
        simulation_material->energy_gen_max = 20000;

        if (verbose)
            std::cout << "\n\nGenerating " << simulation_material->nevents << " events in [" << simulation_material->energy_gen_min << ", " << simulation_material->energy_gen_max << "] GeV\n";

        // Start loop
        auto nthreads = std::thread::hardware_concurrency();

        std::vector<pthread_t> threads (nthreads);
        simulation_material->nevents_per_thread = nevents/threads.size();
        simulation_material->step = static_cast<unsigned int>(simulation_material->nevents_per_thread/10);

        for(size_t th_idx=0; th_idx<threads.size(); ++th_idx) {
            simulation_material->idx = th_idx;
            auto th_result = pthread_create(&threads[th_idx], NULL, buildEvent, (void *) &simulation_material);
            if (th_result) {
                std::cout << "\n\nERROR! return code from pthread_create() is " << th_result << "\n\n";
                exit(-1);
            }
        }

        if (verbose)
            std::cout << "\n\nSimulation has been completed\n\n";

        simulation_material->energy_hist->Sumw2();
        simulation_material->energy_resolution_hist->Sumw2();

        auto energy_ratio_hist = static_cast<TH1D*>(simulation_material->energy_hist->Clone("energy_ratio_hist"));
        energy_ratio_hist->Sumw2();

        energy_ratio_hist->Divide(simulation_material->energy_resolution_hist.get());

        // Create grs
        auto create_gr = [&energy_binning, &energy_nbins, &spectral_index] (TH1D* histo, const char* titleY) -> std::shared_ptr<TGraphErrors> {
            
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
            gr->GetYaxis()->SetTitle(titleY);

            return gr;
        };

        auto gr_energy_hist = create_gr(simulation_material->energy_hist.get(), "counts");
        auto gr_energy_resolution_hist = create_gr(simulation_material->energy_resolution_hist.get(), "counts");
        auto gr_energy_ratio_hist = create_gr(energy_ratio_hist, "counts ratio");

        // Write output
        simulation_material->gr_energy_resolution->Write();
        simulation_material->energy_hist->Write();
        simulation_material->energy_resolution_hist->Write();
        energy_ratio_hist->Write();

        gr_energy_hist->Write();
        gr_energy_resolution_hist->Write();
        gr_energy_ratio_hist->Write();

        output->Close();
    }