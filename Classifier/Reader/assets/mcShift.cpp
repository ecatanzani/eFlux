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
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <ROOT/RDataFrame.hxx>

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

void mcShift(
    const char* input_file,
    const char* output_file,
    const char* energy_config_file,
    const bool verbose = true,
    const bool mc = false,
    const double energy_th = 20,
    const unsigned int threads=1) {

        // Extract energy binning from config file
        auto energy_binning = parse_energy_config(energy_config_file);
        auto energy_nbins = (int)energy_binning.size() - 1;

        double skew_th = 0.6;

        // Read input file
        TFile* infile {TFile::Open(input_file, "READ")};
        if (infile->IsZombie()) {
            std::cerr << "\n\nError opening input file [" << input_file << "]" << std::endl;
            exit(100);
        }

        TIter nextkey(infile->GetListOfKeys());
        TKey *key {nullptr};
        std::unique_ptr<TTree> mytree;
        while ((key=static_cast<TKey*>(nextkey()))) 
        {
            TObject *obj {key->ReadObj()};
            if (obj->IsA()->InheritsFrom(TTree::Class())) {
                mytree = std::unique_ptr<TTree>(static_cast<TTree*>(obj));
                break;
            }
        }

        if (verbose)
            std::cout << "\nFound TTree in input file [" << mytree->GetName() << "]";

        // Initialize RDF
        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame fr(*mytree);

        std::vector<ROOT::RDF::RResultPtr<TH1D>> h_classifier_bin(energy_nbins);

        const double signal_spectral_index = -3;
        auto get_weight = [signal_spectral_index, &energy_binning] (const double energy_gev) -> double {
            return std::pow(energy_gev, signal_spectral_index +1)*std::pow(energy_binning[0], fabs(signal_spectral_index+1));
        };

        // Build the tmva classifier distribution for each bin
        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
            auto bin_filter = [=](int energy_bin) -> bool { return energy_bin == bin_idx; };
            h_classifier_bin[bin_idx-1] = !mc ? fr.Filter(bin_filter, {"energy_bin"}).Histo1D<double>("tmva_classifier") : 
                                                fr.Filter(bin_filter, {"energy_bin"})
                                                    .Define("simu_energy_gev", [](const double energy) -> double {return energy*0.001;}, {"simu_energy"})
                                                    .Define("evt_w", get_weight, {"simu_energy_gev"})
                                                    .Histo1D<double, double>("tmva_classifier", "evt_w");
        }

        std::vector<TF1> gaus_fit_1(energy_nbins), gaus_fit_2(energy_nbins);
        
        // First fit
        double mean {0}, rms {0}, le {0}, he {0};
        for (auto it=std::begin(h_classifier_bin); it != std::end(h_classifier_bin); ++it) {
            mean = it->GetPtr()->GetMean();
            rms = it->GetPtr()->GetRMS();
            double moda {it->GetPtr()->GetBinCenter(it->GetPtr()->GetMaximumBin())};
            if (it->GetPtr()->GetSkewness()>skew_th) {
                le = moda;
                he = moda + 3*rms;
            }
            else {
                le = moda - rms;
                he = moda + rms;
            }
            TF1 gaus_fit("gaus_fit_lv_1", "gaus", le, he);
            !mc ? it->GetPtr()->Fit(&gaus_fit, "RQLN") : it->GetPtr()->Fit(&gaus_fit, "RQWLN");
            gaus_fit_1[std::distance(std::begin(h_classifier_bin), it)] = gaus_fit;
        }

        // Second fit
        double tf1_mean {0}, tf1_rms {0};
        for (auto it=std::begin(h_classifier_bin); it != std::end(h_classifier_bin); ++it) {
            tf1_mean = gaus_fit_1[std::distance(std::begin(h_classifier_bin), it)].GetParameter(1);
            tf1_rms = gaus_fit_1[std::distance(std::begin(h_classifier_bin), it)].GetParameter(2);
            if (it->GetPtr()->GetSkewness()>skew_th) {
                le = tf1_mean;
                he = tf1_mean + 3*tf1_rms;
            }
            else {
                le = tf1_mean-tf1_rms;
                he = tf1_mean+tf1_rms;
            }
            TF1 gaus_fit("gaus_fit_lv_2", "gaus", le, he);
            !mc ? it->GetPtr()->Fit(&gaus_fit, "RQLN") : it->GetPtr()->Fit(&gaus_fit, "RQWLN");
            gaus_fit_2[std::distance(std::begin(h_classifier_bin), it)] = gaus_fit;
        }

        std::vector<double> mean_shift(energy_nbins, 0), mean_shift_error(energy_nbins, 0);
        std::vector<double> rms_shift(energy_nbins, 0), rms_shift_error(energy_nbins, 0);
        std::vector<double> energy_err(energy_nbins, 0), energy(energy_nbins, 0);

        for (int idx=0; idx<energy_nbins; ++idx) {
            energy[idx] = (energy_binning[idx] + energy_binning[idx+1])/2;
        }

        for (int bin_idx = 1; bin_idx <= energy_nbins; ++bin_idx) {
            mean_shift[bin_idx-1] = gaus_fit_2[bin_idx-1].GetParameter(1);
            rms_shift[bin_idx-1] = gaus_fit_2[bin_idx-1].GetParameter(2);
            mean_shift_error[bin_idx -1] = gaus_fit_2[bin_idx-1].GetParError(1);
            rms_shift_error[bin_idx -1] = gaus_fit_2[bin_idx-1].GetParError(2);
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
    const double max_energy_gev = 800) {

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
        gr_shift.GetYaxis()->SetTitle("shift = data peak position - electron peak position");

        TGraphErrors gr_bin_shift(bins.size(), &bins[0], &shift[0], &bins_err[0], &shift_err[0]);
        gr_bin_shift.SetName("gr_bin_shift");
        gr_bin_shift.GetXaxis()->SetTitle("Energy Bin");
        gr_bin_shift.GetYaxis()->SetTitle("shift = data peak position - electron peak position");

        TF1 shift_fit_function("shift_fit_function", "[0]+[1]*log10(x)", energy.front(), energy.back());
        gr_shift.Fit(&shift_fit_function, "Q", "", energy.front(), max_energy_gev);

        TF1 shift_fit_function_bin("shift_fit_function_bin", "pol1", bins.front(), bins.back());
        gr_bin_shift.Fit(&shift_fit_function_bin, "Q", "", bins.front(), 33);

        return std::tuple<TGraphErrors, TF1, TGraphErrors, TF1>(gr_shift, shift_fit_function, gr_bin_shift, shift_fit_function_bin);

    }

std::tuple<TGraphErrors, TF1, TGraphErrors, TF1> fitSigmaRatio(
    const char* mc_file,
    const char* data_file,
    const double max_energy_gev = 800) {

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
            point_ratio_err.push_back(sqrt(pow(gr_data->GetErrorY(idx_p), 2) + pow(gr_electron_mc->GetErrorY(idx_p), 2)));
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

        TF1 sigma_ratio_fit_function_bin("sigma_ratio_fit_function_bin", "pol1", bins.front(), bins.back());
        gr_bin_ratio.Fit(&sigma_ratio_fit_function_bin, "Q", "", bins.front(), 33);

        return std::tuple<TGraphErrors, TF1, TGraphErrors, TF1>(gr_ratio, sigma_ratio_fit_function, gr_bin_ratio, sigma_ratio_fit_function_bin);
    }

void calculateCorrection(
    const char* mc_file,
    const char* data_file,
    const char* output_file,
    const double max_energy_gev = 800) {

        TGraphErrors gr_shift, gr_sigma_ratio;
        TGraphErrors gr_shift_bin, gr_sigma_ratio_bin;
        TF1 shift_fit_func, sigma_ratio_fit_func;
        TF1 shift_fit_func_bin, sigma_ratio_fit_func_bin;

        std::tie(gr_shift, shift_fit_func, gr_shift_bin, shift_fit_func_bin) = fitShift(mc_file, data_file, max_energy_gev);
        std::tie(gr_sigma_ratio, sigma_ratio_fit_func, gr_sigma_ratio_bin, sigma_ratio_fit_func_bin) = fitSigmaRatio(mc_file, data_file, max_energy_gev);

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