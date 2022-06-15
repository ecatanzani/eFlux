#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <tuple>

#include "TF1.h"
#include "TKey.h"
#include "TH1D.h"
#include "TPDF.h"
#include "TPad.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TChain.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TLegendEntry.h"
#include "TFitResultPtr.h"
#include <ROOT/RDataFrame.hxx>

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
        if (!strcmp(key->GetClassName(), "TTree")) {
            if (!strcmp(key->GetName(), "total_tree")) {
                tree_name = static_cast<std::string>(key->GetName());
                break;
            }
        }
    }
    input_file->Close();
    return tree_name;
}

std::shared_ptr<TChain> parse_input_list(const std::string input_list, const bool verbose) {

    auto parse_input_file = [](const std::string input_list) {
        std::ifstream input_file(input_list.c_str());
        if (!input_file.is_open())
        {
            std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
            exit(100);
        }
        std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
        input_file.close();
        return input_string;
    };

    std::istringstream input_stream(parse_input_file(input_list));
    std::shared_ptr<TChain> evtch;
    std::string tmp_str;
    bool first_elm {true};
    while (input_stream >> tmp_str) {
        if (first_elm) {
            evtch = std::make_shared<TChain> (get_tree_name(tmp_str).c_str(), "TMVA data set");
            first_elm = false;
        }
        evtch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }
    return evtch;
}

void mcShift(
    const char* input_file_list,
    const char* output_file,
    const char* energy_config_file,
    const double bdt_cut_value = 0.1,
    const bool mc = false,
    const double energy_th = 20,
    const unsigned int threads=1,
    const bool verbose = true) {

        // Build TChain
        auto chain = parse_input_list(input_file_list, verbose);
        
        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        double skew_th = 0.6;

        // Initialize RDF
        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame fr(*chain);

        std::vector<ROOT::RDF::RResultPtr<TH1D>> h_classifier_bin (energy_nbins);
        std::vector<TH1D*> h_classifier_bin_proton_subtracted (energy_nbins);
        std::vector<TFitResultPtr> h_classifier_bin_fit_result (energy_nbins);

        const double signal_spectral_index = -3;
        auto get_weight = [signal_spectral_index, &energy_binning] (const double energy_gev) -> double {
            return std::pow(energy_gev, signal_spectral_index +1)*std::pow(energy_binning[0], fabs(signal_spectral_index+1));
        };

        // Build the tmva classifier distribution for each bin
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
            auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };
            std::string histo_name = "h_classifier_bin_" + std::to_string(bin_idx);

            if (mc)
                h_classifier_bin[bin_idx-1] = fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("corr_energy_gev", "energy_corr * 0.001")
                                                    .Define("evt_w", get_weight, {"corr_energy_gev"})
                                                    .Histo1D<double, double>({histo_name.c_str(), histo_name.c_str(), 1000, -1, 1}, "tmva_classifier", "evt_w");  
            else {  
                int nbins {1};
                if (bin_idx<=28)
                    nbins = 1000;
                else
                    nbins = 500;

                h_classifier_bin[bin_idx-1] = fr.Filter(bin_filter, {"energy_bin"}).Histo1D<double>({histo_name.c_str(), histo_name.c_str(), nbins, -1, 1}, "tmva_classifier");  
            }                         
        }
        
        std::vector<TF1> gaus_fit_1(energy_nbins), gaus_fit_2(energy_nbins);
        std::vector<TF1> data_proton_linear_fit(energy_nbins);
        
        // Fit proton slope on data
        if (!mc) {
            for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
                std::string tf1_name = "proton_linear_fit_" + std::to_string(bin_idx);
                std::string fit_result_name = "proton_linear_fit_result_" + std::to_string(bin_idx);
                std::string histo_name = std::string(h_classifier_bin[bin_idx-1]->GetName()) + "_proton_subtracted";
                // Define the fit function in the whole BDT range [-1, 1]
                data_proton_linear_fit[bin_idx-1] = TF1(tf1_name.c_str(), "pow(10, [0]+[1]*x)", -1, 1);
                // Clone the data histo - cloned one will be used to subtract proton background TF1
                h_classifier_bin_proton_subtracted[bin_idx-1] = static_cast<TH1D*>(h_classifier_bin[bin_idx-1]->Clone(histo_name.c_str()));
                
                // Set the fit interval accordingly to the energy bin

                /*
                OLD METHOD
                double le {-0.2};
                double he {0};

                if (bin_idx<18) {
                    le = -0.2;
                    he = 0;
                }
                else if (bin_idx>=18 && bin_idx<24) {
                    le = -0.4;
                    he = 0;
                }

                // Hign energy bins
                else if (bin_idx==28) {
                    // Fix this
                    le = -0.3;
                    he = -0.2;
                }
                else if (bin_idx==29) {
                    // Fix this
                    le = -0.3;
                    he = -0.1;
                }
                else if (bin_idx==30) {
                    // Fix this
                    le = -0.3;
                    he = -0.15;
                }
                else if (bin_idx==31) {
                    // Fix this
                    le = -0.3;
                    he = -0.15;
                }
                else if (bin_idx==32) {
                    le = -0.2;
                    he = 0;
                }
                else if (bin_idx==34) {
                    le = -0.25;
                    he = -0.1;
                }
                else if (bin_idx==35) {
                    le = -0.3;
                    he = -0.;
                }
                else if (bin_idx==36) {
                    le = -0.3;
                    he = -0.1;
                }
                else if (bin_idx==37) {
                    le = -0.25;
                    he = -0.15;
                }
                else if (bin_idx==38) {
                    le = -0.3;
                    he = -0.1;
                }
                else if (bin_idx==39) {
                    le = -0.4;
                    he = -0.2;
                }
                else {
                    le = -0.4;
                    he = -0.2;
                }
                */

                double le {-0.2};
                double he {0};

                if (bin_idx<18) {
                    le = -0.2;
                    he = 0;
                }
                else if (bin_idx>=18 && bin_idx<24) {
                    le = -0.4;
                    he = 0;
                }

                // Hign energy bins
                else if (bin_idx==28) {
                    // Fix this
                    le = -0.3;
                    he = -0.2;
                }
                else if (bin_idx==29) {
                    // Fix this
                    le = -0.3;
                    he = -0.1;
                }
                else if (bin_idx==30) {
                    // Fix this
                    le = -0.2;
                    he = 0.2;
                }
                else if (bin_idx==31) {
                    // Fix this
                    le = -0.3;
                    he = -0.15;
                }
                else if (bin_idx==32) {
                    le = -0.2;
                    he = 0.2;
                }
                else if (bin_idx==34) {
                    le = -0.25;
                    he = -0.1;
                }
                else if (bin_idx==35) {
                    le = -0.2;
                    he = 0.2;
                }
                else if (bin_idx==36) {
                    le = -0.4;
                    he = 0;
                }
                else if (bin_idx==37) {
                    le = -0.25;
                    he = -0.15;
                }
                else if (bin_idx==38) {
                    le = -0.4;
                    he = -0.2;
                }
                else if (bin_idx==39) {
                    le = -0.4;
                    he = -0.2;
                }
                else {
                    le = -0.4;
                    he = -0.2;
                }

                TF1 tmp_function("tmp_function", "pow(10, [0]+[1]*x)", -1, 1);
                h_classifier_bin_proton_subtracted[bin_idx-1]->Fit(&tmp_function, "SQN", "", le, he);
                // Move the parameters of the tmp function to the final TF1 as starting parameters
                for (int pidx=0; pidx<tmp_function.GetNpar(); ++pidx)
                    data_proton_linear_fit[bin_idx-1].SetParameter(pidx, tmp_function.GetParameter(pidx));
                // Fit the energy bin with the log-likelihood
                h_classifier_bin_fit_result[bin_idx-1] = static_cast<TFitResultPtr>(h_classifier_bin_proton_subtracted[bin_idx-1]->Fit(&data_proton_linear_fit[bin_idx-1], "SQNL", "", le, he));

                h_classifier_bin_fit_result[bin_idx-1]->SetName(fit_result_name.c_str());
                h_classifier_bin_proton_subtracted[bin_idx-1]->Add(&data_proton_linear_fit[bin_idx-1], -1);
            }
        }

        // First fit
        double mean {0}, rms {0}, le {0}, he {0};
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {

            auto h_pointer = mc ? static_cast<TH1D*>(h_classifier_bin[bin_idx-1]->Clone()) : h_classifier_bin_proton_subtracted[bin_idx-1];

            // Set the correct range in BDT classifier (we want to fit signal peak only)
            h_pointer->GetXaxis()->SetRangeUser(bdt_cut_value, 1);

            mean = h_pointer->GetMean();
            rms = h_pointer->GetRMS();
            double moda {h_pointer->GetBinCenter(h_pointer->GetMaximumBin())};

            if (h_pointer->GetSkewness()>skew_th) {
                le = moda;
                he = moda + 3*rms;
            }
            else {
                le = moda - rms;
                he = moda + rms;
            }

            // Fit
            TF1 gaus_fit("gaus_fit_lv_1", "gaus", le, he);
            !mc ? h_pointer->Fit(&gaus_fit, "RQLN") : h_pointer->Fit(&gaus_fit, "RQWLN");
            gaus_fit_1[bin_idx-1] = gaus_fit;

        }

        // Second fit
        double tf1_mean {0}, tf1_rms {0};
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {

            auto h_pointer = mc ? static_cast<TH1D*>(h_classifier_bin[bin_idx-1]->Clone()) : h_classifier_bin_proton_subtracted[bin_idx-1];
            
            tf1_mean = gaus_fit_1[bin_idx-1].GetParameter(1);
            tf1_rms = gaus_fit_1[bin_idx-1].GetParameter(2);
            if (h_pointer->GetSkewness()>skew_th) {
                le = tf1_mean;
                he = tf1_mean + 3*tf1_rms;
            }
            else {
                le = tf1_mean-tf1_rms;
                he = tf1_mean+tf1_rms;
            }

            // Fit
            TF1 gaus_fit("gaus_fit_lv_2", "gaus", le, he);
            !mc ? h_pointer->Fit(&gaus_fit, "RQLN") : h_pointer->Fit(&gaus_fit, "RQWLN");
            gaus_fit_2[bin_idx-1] = gaus_fit;

        }

        std::vector<double> mean_shift(energy_nbins, 0), mean_shift_error(energy_nbins, 0);
        std::vector<double> rms_shift(energy_nbins, 0), rms_shift_error(energy_nbins, 0);
        std::vector<double> energy_err(energy_nbins, 0), energy(energy_nbins, 0);

        for (int idx=0; idx<energy_nbins; ++idx) {
            energy[idx] = (energy_binning[idx] + energy_binning[idx+1])/2;
        }

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
           
            if (mc || (!mc && bin_idx<28)) {
                mean_shift[bin_idx-1] = gaus_fit_2[bin_idx-1].GetParameter(1);
                rms_shift[bin_idx-1] = gaus_fit_2[bin_idx-1].GetParameter(2);
                mean_shift_error[bin_idx -1] = gaus_fit_2[bin_idx-1].GetParError(1);
                rms_shift_error[bin_idx -1] = gaus_fit_2[bin_idx-1].GetParError(2);
            }
            else {
                // Once the peak is fitted, the mean value and the RMS are extracted directly from the histogram from the fit mean value +/- 3 fit sigmas
                auto mean_value = gaus_fit_2[bin_idx-1].GetParameter(1);
                auto sigma_value = gaus_fit_2[bin_idx-1].GetParameter(2);
                h_classifier_bin_proton_subtracted[bin_idx-1]->GetXaxis()->SetRangeUser(mean_value - 3*sigma_value, mean_value + 3*sigma_value);

                mean_shift[bin_idx-1] = h_classifier_bin_proton_subtracted[bin_idx-1]->GetMean();
                rms_shift[bin_idx-1] = h_classifier_bin_proton_subtracted[bin_idx-1]->GetRMS();
                mean_shift_error[bin_idx -1] = h_classifier_bin_proton_subtracted[bin_idx-1]->GetMeanError();
                rms_shift_error[bin_idx -1] = h_classifier_bin_proton_subtracted[bin_idx-1]->GetRMSError();
            }
        }

        std::vector<double> bins(energy_nbins);
        std::iota(std::begin(bins), std::end(bins), 1);
        
        // Mean graph for bin
        TGraph gr_mean(energy_nbins, &bins[0], &mean_shift[0]);
        gr_mean.SetName("gr_mean");
        gr_mean.GetXaxis()->SetTitle("Energy Bin");
        !mc ? gr_mean.GetYaxis()->SetTitle("data peak position") : gr_mean.GetYaxis()->SetTitle("electron peak position");

        // Mean graph with errors for bin
        TGraphErrors gr_mean_we(energy_nbins, &bins[0], &mean_shift[0], &energy_err[0], &mean_shift_error[0]);
        gr_mean_we.SetName("gr_mean_we");
        gr_mean_we.GetXaxis()->SetTitle("Energy Bin");
        !mc ? gr_mean_we.GetYaxis()->SetTitle("data peak position") : gr_mean_we.GetYaxis()->SetTitle("electron peak position");

        // Sigma graph for bin
        TGraph gr_sigma(energy_nbins, &bins[0], &rms_shift[0]);
        gr_sigma.SetName("gr_sigma");
        gr_sigma.GetXaxis()->SetTitle("Energy Bin");
        !mc ? gr_sigma.GetYaxis()->SetTitle("#sigma_{data}") : gr_sigma.GetYaxis()->SetTitle("#sigma_{electron}");

        // Sigma graph with errors for bin
        TGraphErrors gr_sigma_we(energy_nbins, &bins[0], &rms_shift[0], &energy_err[0], &rms_shift_error[0]);
        gr_sigma_we.SetName("gr_sigma_we");
        gr_sigma_we.GetXaxis()->SetTitle("Energy Bin");
        !mc ? gr_sigma_we.GetYaxis()->SetTitle("#sigma_{data}") : gr_sigma_we.GetYaxis()->SetTitle("#sigma_{electron}");

        // Mean graph full energy range
        TGraphErrors gr_mean_full_interval(energy_nbins, &energy[0], &mean_shift[0], &energy_err[0], &mean_shift_error[0]);
        gr_mean_full_interval.SetName("gr_mean_full_interval");
        gr_mean_full_interval.GetXaxis()->SetTitle("Energy [GeV]");
        !mc ? gr_mean_full_interval.GetYaxis()->SetTitle("data peak position") : gr_mean_full_interval.GetYaxis()->SetTitle("electron peak position");
        TF1 mean_full_interval_fit_func("mean_full_interval_fit_func", "[0] + [1]*log10(x)", energy.front(), energy.back());
        gr_mean_full_interval.Fit(&mean_full_interval_fit_func, "RQ");

        // Sigma graph full energy range
        TGraphErrors gr_sigma_full_interval(energy_nbins, &energy[0], &rms_shift[0], &energy_err[0], &rms_shift_error[0]);
        gr_sigma_full_interval.SetName("gr_sigma_full_interval");
        gr_sigma_full_interval.GetXaxis()->SetTitle("Energy [GeV]");
        !mc ? gr_sigma_full_interval.GetYaxis()->SetTitle("#sigma_{data}") : gr_sigma_full_interval.GetYaxis()->SetTitle("#sigma_{electron}");
        TF1 sigma_full_interval_fit_func("sigma_full_interval_fit_func", "[0] + [1]*log10(x)", energy.front(), energy.back());
        gr_sigma_full_interval.Fit(&sigma_full_interval_fit_func, "RQ");

        std::vector<double> mean_shift_shrink, mean_shift_err_shrink;
        std::vector<double> sigma_shift_shrink, sigma_shift_err_shrink;
        std::vector<double> energy_shrink, energy_err_shrink;
        
        for (int idx=0; idx<energy_nbins; ++idx) {
            if (energy[idx]>energy_th) {
                energy_shrink.push_back(energy[idx]);
                energy_err_shrink.push_back(0);
                mean_shift_shrink.push_back(mean_shift[idx]);
                mean_shift_err_shrink.push_back(mean_shift_error[idx]);
                sigma_shift_shrink.push_back(rms_shift[idx]);
                sigma_shift_err_shrink.push_back(rms_shift_error[idx]);
            }
        }

        // Mean graph small energy range
        TGraphErrors gr_mean_reduced_energy_range((int)energy_shrink.size(), &energy_shrink[0], &mean_shift_shrink[0], &energy_err_shrink[0], &mean_shift_err_shrink[0]);
        gr_mean_reduced_energy_range.SetName("gr_mean_reduced_energy_range");
        gr_mean_reduced_energy_range.GetXaxis()->SetTitle("Energy [GeV]");
        !mc ? gr_mean_reduced_energy_range.GetYaxis()->SetTitle("data peak position") : gr_mean_reduced_energy_range.GetYaxis()->SetTitle("electron peak position");
        TF1 mean_fit_func("mean_fit_func", "[0] + [1]*log10(x)", energy_shrink.front(), energy_shrink.back());
        gr_mean_reduced_energy_range.Fit(&mean_fit_func, "RQ");

        // Sigma graph small energy range
        TGraphErrors gr_sigma_reduced_energy_range((int)energy_shrink.size(), &energy_shrink[0], &sigma_shift_shrink[0], &energy_err_shrink[0], &sigma_shift_err_shrink[0]);
        gr_sigma_reduced_energy_range.SetName("gr_sigma_reduced_energy_range");
        gr_sigma_reduced_energy_range.GetXaxis()->SetTitle("Energy [GeV]");
        !mc ? gr_sigma_reduced_energy_range.GetYaxis()->SetTitle("#sigma_{data}") : gr_sigma_reduced_energy_range.GetYaxis()->SetTitle("#sigma_{electron}");
        TF1 sigma_fit_func("sigma_fit_func", "[0] + [1]*log10(x)", energy_shrink.front(), energy_shrink.back());
        gr_sigma_reduced_energy_range.Fit(&sigma_fit_func, "RQ");

        // Write output file
        TFile outfile(output_file, "RECREATE");
        if (outfile.IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        for (int bidx = 0; bidx < energy_nbins; ++bidx) {
            auto tmp_dir_name = std::string("energybin_") + std::to_string(bidx + 1);
            outfile.mkdir(tmp_dir_name.c_str());
            outfile.cd(tmp_dir_name.c_str());
            h_classifier_bin[bidx]->Write();
            if (!mc) {
                h_classifier_bin_proton_subtracted[bidx]->GetXaxis()->SetRangeUser(-1, 1);
                h_classifier_bin_proton_subtracted[bidx]->Write();
                data_proton_linear_fit[bidx].Write();
                h_classifier_bin_fit_result[bidx]->Write();
            }
            gaus_fit_1[bidx].Write();
            gaus_fit_2[bidx].Write();
        }

        outfile.cd();

        gr_mean.Write();
        gr_mean_we.Write();
        gr_sigma.Write();
        gr_sigma_we.Write();

        gr_mean_full_interval.Write();
        gr_sigma_full_interval.Write();

        gr_mean_reduced_energy_range.Write();  
        gr_sigma_reduced_energy_range.Write();

        mean_full_interval_fit_func.Write();
        mean_fit_func.Write();
        sigma_full_interval_fit_func.Write(); 
        sigma_fit_func.Write();

        outfile.Close();

        // Save linear background functions for each energy bin
        if (!mc) {
            TFile function_outfile("proton_background_function_output.root", "RECREATE");
            if (function_outfile.IsZombie()) {
                std::cerr << "\n\nError writing output ROOT file [" << "proton_background_function_output.root" << "]\n\n";
                exit(100);
            }

            for (int bidx = 0; bidx < energy_nbins; ++bidx) {
                std::string tmp_dir_name = std::string("energybin_") + std::to_string(bidx + 1);
                function_outfile.mkdir(tmp_dir_name.c_str());
                function_outfile.cd(tmp_dir_name.c_str());
                data_proton_linear_fit[bidx].Write();
                h_classifier_bin[bidx]->Write();
            }

            function_outfile.Close();
        }

        // Create PDF file for each energy bin
        
        // Build canvas
        TCanvas print_canvas("print_canvas", "print_canvas");
        print_canvas.SetTicks();

        TPaveLabel label(0.0, 0.95, 0.3, 1, "BDT classifier", "tlNDC");

        for (int bidx = 0; bidx < energy_nbins; ++bidx) {
            
            if (mc) {
                h_classifier_bin[bidx]->SetLineWidth(2);
                h_classifier_bin[bidx]->SetLineColor(kBlue);
                gaus_fit_1[bidx].SetLineWidth(2);
                gaus_fit_1[bidx].SetLineColor(kRed-9);
                gaus_fit_2[bidx].SetLineWidth(2);
                gaus_fit_2[bidx].SetLineColor(kRed+1);

                h_classifier_bin[bidx]->Draw();
                gaus_fit_1[bidx].Draw("same");
                gaus_fit_2[bidx].Draw("same");

            }
            else {
                h_classifier_bin[bidx]->SetLineWidth(2);
                h_classifier_bin[bidx]->SetLineColor(kBlack);
                h_classifier_bin_proton_subtracted[bidx]->SetLineWidth(2);
                h_classifier_bin_proton_subtracted[bidx]->SetLineColor(kBlue);
                data_proton_linear_fit[bidx].SetLineWidth(2);
                data_proton_linear_fit[bidx].SetLineColor(kGreen);
                gaus_fit_1[bidx].SetLineWidth(2);
                gaus_fit_1[bidx].SetLineColor(kRed-9);
                gaus_fit_2[bidx].SetLineWidth(2);
                gaus_fit_2[bidx].SetLineColor(kRed+1);

                h_classifier_bin[bidx]->Draw();
                h_classifier_bin_proton_subtracted[bidx]->Draw("same");
                data_proton_linear_fit[bidx].Draw("same");
                gaus_fit_1[bidx].Draw("same");
                gaus_fit_2[bidx].Draw("same");
            }

            gPad->SetLogy();
            gPad->SetGrid(1,1);
            gStyle->SetOptStat(0);
            
            std::string label_name = "BDT classifier - energy bin " + std::to_string(bidx+1) + " - [" + std::to_string(energy_binning[bidx]) + ", " + std::to_string(energy_binning[bidx+1]) + "] GeV";
            label = TPaveLabel(0.0, 0.95, 0.3, 1, label_name.c_str(), "tlNDC");
            label.Draw();
            gStyle->SetOptTitle(0);
            if (!bidx)
                print_canvas.Print("bdt_classifier_summary.pdf(","Title:BDT classifier");
            else if (bidx==(energy_nbins-1))
                print_canvas.Print("bdt_classifier_summary.pdf)","Title:BDT classifier");
            else
                print_canvas.Print("bdt_classifier_summary.pdf","Title:BDT classifier");
        }   
    }


inline std::shared_ptr<TGraphErrors> extractGR(const char* file, const char* name) {

    TFile* infile = TFile::Open(file, "READ");
    if (infile->IsZombie()) {
        std::cerr << "\n\nError reading input file [" << file << "]" << std::endl;
        exit(100);
    }

    std::shared_ptr<TGraphErrors> gr = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(name)));
    return gr;
}

std::tuple<TGraphErrors, TF1, TGraphErrors, TF1> fitPeakPositionDifference(
    const char* mc_file,
    const char* data_file,
    const double max_energy_gev = 800) {

        auto gr_electron_mc = extractGR(mc_file, "gr_mean_full_interval");
        auto gr_data = extractGR(data_file, "gr_mean_full_interval");

        std::vector<double> energy, energy_err;
        std::vector<double> bins, bins_err;
        std::vector<double> point_difference, point_difference_err;

        for (int idx_p=0; idx_p<gr_electron_mc->GetN(); ++idx_p) {
            energy.push_back(gr_data->GetPointX(idx_p));
            bins.push_back(idx_p+1);
            energy_err.push_back(0);
            bins_err.push_back(0);
            point_difference.push_back(gr_data->GetPointY(idx_p) - gr_electron_mc->GetPointY(idx_p));
            point_difference_err.push_back(sqrt(pow(gr_data->GetErrorY(idx_p), 2) + pow(gr_electron_mc->GetErrorY(idx_p), 2)));
        }

        TGraphErrors gr_difference(energy.size(), &energy[0], &point_difference[0], &energy_err[0], &point_difference_err[0]);
        gr_difference.SetName("gr_difference");
        gr_difference.GetXaxis()->SetTitle("Energy [GeV]");
        gr_difference.GetYaxis()->SetTitle("shift = data peak position - electron peak position");

        TGraphErrors gr_bin_difference(bins.size(), &bins[0], &point_difference[0], &bins_err[0], &point_difference_err[0]);
        gr_difference.SetName("gr_bin_difference");
        gr_difference.GetXaxis()->SetTitle("Energy Bin");
        gr_difference.GetYaxis()->SetTitle("shift = data peak position - electron peak position");

        TF1 shift_fit_function("shift_fit_function", "[0]+[1]*log10(x)", energy.front(), energy.back());
        gr_difference.Fit(&shift_fit_function, "Q", "", energy.front(), max_energy_gev);

        TF1 shift_fit_function_bin("shift_fit_function_bin", "pol1", bins.front(), bins.back());
        gr_bin_difference.Fit(&shift_fit_function_bin, "Q", "", bins.front(), 33);

        return std::tuple<TGraphErrors, TF1, TGraphErrors, TF1>(gr_difference, shift_fit_function, gr_bin_difference, shift_fit_function_bin);

    }

std::tuple<TGraphErrors, TF1, TGraphErrors, TF1> fitShift(
    const char* mc_file,
    const char* data_file,
    const double max_energy_gev,
    const int max_energy_bin) {

        auto gr_mean_electron_mc = extractGR(mc_file, "gr_mean_full_interval");
        auto gr_mean_data = extractGR(data_file, "gr_mean_full_interval");
        auto gr_sigma_electron_mc = extractGR(mc_file, "gr_sigma_full_interval");
        auto gr_sigma_data = extractGR(data_file, "gr_sigma_full_interval");

        std::vector<double> energy, energy_err;
        std::vector<double> bins, bins_err;
        std::vector<double> shift, shift_err;

        for (int idx_p=0; idx_p<gr_mean_electron_mc->GetN(); ++idx_p) {
            energy.push_back(gr_mean_data->GetPointX(idx_p));
            bins.push_back(idx_p+1);
            energy_err.push_back(0);
            bins_err.push_back(0);
            shift.push_back(gr_mean_data->GetPointY(idx_p) - (gr_sigma_data->GetPointY(idx_p)/gr_sigma_electron_mc->GetPointY(idx_p))*gr_mean_electron_mc->GetPointY(idx_p));
            shift_err.push_back(sqrt(
                pow(gr_mean_data->GetErrorY(idx_p), 2) + 
                pow((gr_sigma_data->GetPointY(idx_p)/gr_sigma_electron_mc->GetPointY(idx_p))*gr_mean_electron_mc->GetErrorY(idx_p), 2) +
                pow((gr_mean_electron_mc->GetPointY(idx_p)/gr_sigma_electron_mc->GetPointY(idx_p))*gr_sigma_data->GetErrorY(idx_p) , 2) + 
                pow(((gr_sigma_data->GetPointY(idx_p)*gr_mean_electron_mc->GetPointY(idx_p))/gr_sigma_electron_mc->GetPointY(idx_p))*gr_sigma_electron_mc->GetErrorY(idx_p) , 2)));
        }

        TGraphErrors gr_shift(energy.size(), &energy[0], &shift[0], &energy_err[0], &shift_err[0]);
        gr_shift.SetName("gr_shift");
        gr_shift.GetXaxis()->SetTitle("Energy [GeV]");
        gr_shift.GetYaxis()->SetTitle("shift = data peak position - (#sigma_{data}/#sigma_{mc}) electron peak position");

        TGraphErrors gr_bin_shift(bins.size(), &bins[0], &shift[0], &bins_err[0], &shift_err[0]);
        gr_bin_shift.SetName("gr_bin_shift");
        gr_bin_shift.GetXaxis()->SetTitle("Energy Bin");
        gr_bin_shift.GetYaxis()->SetTitle("shift = data peak position - (#sigma_{data}/#sigma_{mc}) electron peak position");

        TF1 shift_fit_function("shift_fit_function", "[0]+[1]*log10(x)", energy.front(), energy.back());
        gr_shift.Fit(&shift_fit_function, "Q", "", energy.front(), max_energy_gev);
        //gr_shift.Fit(&shift_fit_function, "Q", "", 25, max_energy_gev);

        TF1 shift_fit_function_bin("shift_fit_function_bin", "pol1", bins.front(), bins.back());
        gr_bin_shift.Fit(&shift_fit_function_bin, "Q", "", bins.front(), max_energy_bin);
        //gr_bin_shift.Fit(&shift_fit_function_bin, "Q", "", 8, 33);

        return std::tuple<TGraphErrors, TF1, TGraphErrors, TF1>(gr_shift, shift_fit_function, gr_bin_shift, shift_fit_function_bin);

    }

std::tuple<TGraphErrors, TF1, TGraphErrors, TF1> fitSigmaRatio(
    const char* mc_file,
    const char* data_file,
    const double max_energy_gev,
    const int max_energy_bin) {

        auto gr_electron_mc = extractGR(mc_file, "gr_sigma_full_interval");
        auto gr_data = extractGR(data_file, "gr_sigma_full_interval");

        std::vector<double> energy, energy_err;
        std::vector<double> bins, bins_err;
        std::vector<double> point_ratio, point_ratio_err;


        for (int idx_p=0; idx_p<gr_electron_mc->GetN(); ++idx_p) {
            energy.push_back(gr_data->GetPointX(idx_p));
            bins.push_back(idx_p+1);
            energy_err.push_back(0);
            bins_err.push_back(0);
            point_ratio.push_back(gr_data->GetPointY(idx_p) / gr_electron_mc->GetPointY(idx_p));
            //point_ratio_err.push_back(sqrt(pow(gr_data->GetErrorY(idx_p), 2) + pow(gr_electron_mc->GetErrorY(idx_p), 2)));
            point_ratio_err.push_back((1./gr_electron_mc->GetPointY(idx_p))*sqrt(pow(gr_data->GetErrorY(idx_p), 2) + pow(gr_data->GetPointY(idx_p)*gr_electron_mc->GetErrorY(idx_p), 2)));
        }

        TGraphErrors gr_ratio(energy.size(), &energy[0], &point_ratio[0], &energy_err[0], &point_ratio_err[0]);
        gr_ratio.SetName("gr_ratio");
        gr_ratio.GetXaxis()->SetTitle("Energy [GeV]");
        gr_ratio.GetYaxis()->SetTitle("#sigma_{data}/#sigma_{electron}");

        TGraphErrors gr_bin_ratio(bins.size(), &bins[0], &point_ratio[0], &bins_err[0], &point_ratio_err[0]);
        gr_bin_ratio.SetName("gr_bin_ratio");
        gr_bin_ratio.GetXaxis()->SetTitle("Energy Bin");
        gr_bin_ratio.GetYaxis()->SetTitle("#sigma_{data}/#sigma_{electron}");

        TF1 sigma_ratio_fit_function("sigma_ratio_fit_function", "[0]+[1]*log10(x)", energy.front(), energy.back());
        gr_ratio.Fit(&sigma_ratio_fit_function, "Q", "", energy.front(), max_energy_gev);
        //gr_ratio.Fit(&sigma_ratio_fit_function, "Q", "", 25, max_energy_gev);

        TF1 sigma_ratio_fit_function_bin("sigma_ratio_fit_function_bin", "pol1", bins.front(), bins.back());
        gr_bin_ratio.Fit(&sigma_ratio_fit_function_bin, "Q", "", bins.front(), max_energy_bin);
        //gr_bin_ratio.Fit(&sigma_ratio_fit_function_bin, "Q", "", 8, 33);

        return std::tuple<TGraphErrors, TF1, TGraphErrors, TF1>(gr_ratio, sigma_ratio_fit_function, gr_bin_ratio, sigma_ratio_fit_function_bin);
    }

void calculateCorrection(
    const char* mc_file,
    const char* data_file,
    const char* output_file,
    const double max_energy_gev = 1000,
    const int max_energy_bin = 27) {

        TGraphErrors gr_shift, gr_sigma_ratio;
        TGraphErrors gr_shift_bin, gr_sigma_ratio_bin;
        TF1 shift_fit_func, sigma_ratio_fit_func;
        TF1 shift_fit_func_bin, sigma_ratio_fit_func_bin;

        std::tie(gr_shift, shift_fit_func, gr_shift_bin, shift_fit_func_bin) = fitShift(mc_file, data_file, max_energy_gev, max_energy_bin);
        std::tie(gr_sigma_ratio, sigma_ratio_fit_func, gr_sigma_ratio_bin, sigma_ratio_fit_func_bin) = fitSigmaRatio(mc_file, data_file, max_energy_gev, max_energy_bin);

        TFile outfile(output_file, "RECREATE");
        if (outfile.IsZombie()) {
            std::cerr << "\n\nError writing output file [" << output_file << "]" << std::endl;
            exit(100);
        }

        gr_shift.Write();
        gr_shift_bin.Write();
        gr_sigma_ratio.Write();
        gr_sigma_ratio_bin.Write();

        shift_fit_func.Write();
        shift_fit_func_bin.Write();
        sigma_ratio_fit_func.Write();
        sigma_ratio_fit_func_bin.Write();

        outfile.Close();
    }