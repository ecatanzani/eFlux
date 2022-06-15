#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>

#include "TF1.h"
#include "TKey.h"
#include "TH1D.h"
#include "TFile.h"
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

inline double wtsydp(const double minene, const double maxene, const double index)
{
    auto dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}

void buildFlux(
    const char* e_counts_input_file,
    const char* e_acc_input_file,
    const char* signal_efficiency,
    const char* background_contamination,
    const char* energy_config_file,
    const char* bdt_best_point_syst_file,
    const char* efficiency_syst_file,
    const double exposure_time = 140143948.911,
    const char* output_file = "flux.root",
    const bool add_systematics = false,
    const bool verbose = true) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        // Extract info from input files
        auto h_e_counts     = extractHistoFromFile(e_counts_input_file, verbose);
        auto h_e_acc        = extractHistoFromFile(e_acc_input_file, verbose);
        auto h_signal_eff   = extractHistoFromFile(signal_efficiency, verbose);
        auto h_back_cont    = extractHistoFromFile(background_contamination, verbose);

        h_e_counts->Sumw2();
        h_e_acc->Sumw2();
        h_signal_eff->Sumw2();
        h_back_cont->Sumw2();

        // Subtract the background contamination
        h_e_counts->Add(h_back_cont.get(), -1);
        
        // Correct for the efficiency
        h_e_counts->Divide(h_signal_eff.get());

        // Divide by the acceptance
        h_e_counts->Divide(h_e_acc.get());

        // Scale by the exposure-time
        h_e_counts->Scale(1/exposure_time);

        // divide by the energy bin width
        for (int bIdx {1}; bIdx<=h_e_counts->GetNbinsX(); ++bIdx) {
            if(h_e_counts->GetBinContent(bIdx)) {
                h_e_counts->SetBinContent(bIdx, static_cast<double>(h_e_counts->GetBinContent(bIdx))/h_e_counts->GetBinWidth(bIdx));
                h_e_counts->SetBinError(bIdx, static_cast<double>(h_e_counts->GetBinError(bIdx))/h_e_counts->GetBinWidth(bIdx));
            }
        }

        // Correct for the efficiency systematics
        auto h_efficiency_syst = extractHistoFromFile(efficiency_syst_file, verbose);
        
        for (int bIdx {1}; bIdx<=h_e_counts->GetNbinsX(); ++bIdx)
            h_e_counts->SetBinError(bIdx, sqrt(pow(h_e_counts->GetBinError(bIdx), 2) + pow(h_e_counts->GetBinContent(bIdx)*h_efficiency_syst->GetBinError(bIdx)/100., 2)));

        // Build flux multiplied by E^3
        auto h_e_counts_E3 = static_cast<TH1D*>(h_e_counts->Clone("h_e_counts_E3"));
        for (int bIdx {1}; bIdx<=h_e_counts_E3->GetNbinsX(); ++bIdx) {
            if(h_e_counts_E3->GetBinContent(bIdx)) {
                h_e_counts_E3->SetBinContent(bIdx, h_e_counts_E3->GetBinContent(bIdx)*pow(h_e_counts_E3->GetXaxis()->GetBinCenter(bIdx), 3));
                h_e_counts_E3->SetBinError(bIdx, h_e_counts_E3->GetBinError(bIdx)*pow(h_e_counts_E3->GetXaxis()->GetBinCenter(bIdx), 3));
            }
        }

        // Correct for the BDT best point efficiency systematic
        auto h_bdt_bestpoint_syst = extractHistoFromFile(bdt_best_point_syst_file, verbose);
        
        for (int bIdx {1}; bIdx<=h_e_counts->GetNbinsX(); ++bIdx)
            h_e_counts_E3->SetBinError(bIdx, sqrt(pow(h_e_counts_E3->GetBinError(bIdx), 2) + pow(h_bdt_bestpoint_syst->GetBinError(bIdx), 2)));

        // Build TGraphs
        std::vector<double> energy                          (h_e_counts->GetNbinsX(), 0);
        std::vector<double> energy_err                      (h_e_counts->GetNbinsX(), 0);

        std::vector<double> flux                            (h_e_counts->GetNbinsX(), 0);
        std::vector<double> flux_err                        (h_e_counts->GetNbinsX(), 0);

        std::vector<double> flux_E3                         (h_e_counts->GetNbinsX(), 0);
        std::vector<double> flux_E3_err                     (h_e_counts->GetNbinsX(), 0);

        for (unsigned int idx=0; idx<energy.size(); ++idx) {
            energy[idx]         = wtsydp(energy_binning[idx], energy_binning[idx+1], -3);
            flux[idx]           = h_e_counts->GetBinContent(idx+1);
            flux_err[idx]       = h_e_counts->GetBinError(idx+1);
            flux_E3[idx]        = h_e_counts_E3->GetBinContent(idx+1);
            flux_E3_err[idx]    = h_e_counts_E3->GetBinError(idx+1);
        }

        TGraphErrors gr_flux(energy.size(), &energy[0], &flux[0], &energy_err[0], &flux_err[0]);
        TGraphErrors gr_flux_E3(energy.size(), &energy[0], &flux_E3[0], &energy_err[0], &flux_E3_err[0]);

        gr_flux.SetName("gr_flux");
        gr_flux.SetTitle("Flux");
        gr_flux.GetXaxis()->SetTitle("Energy [GeV]");
        gr_flux.GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

        gr_flux_E3.SetName("gr_flux_E3");
        gr_flux_E3.SetTitle("Flux");
        gr_flux_E3.GetXaxis()->SetTitle("Energy [GeV]");
        gr_flux_E3.GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

        // Write output file
        TFile *outfile {TFile::Open(output_file, "RECREATE")};
        if (outfile->IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        h_e_counts->SetName("h_all_e_flux");
        h_e_counts->SetTitle("All-Electron DAMPE flux");
        h_e_counts->GetXaxis()->SetTitle("Energy [GeV]");
        h_e_counts->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

        h_e_counts_E3->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

        outfile->cd();

        h_e_counts->Write();
        h_e_counts_E3->Write();

        gr_flux.Write();
        gr_flux_E3.Write();

        outfile->Close();

        if (verbose)
            std::cout << "\n\nOutput file has been written ... [" << output_file << "]\n\n";
    }