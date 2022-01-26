#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TH1D.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

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

inline std::vector<float> shrink_binning(const std::vector<float>& complete_binning, const double energy_th) {
    std::vector<float> new_binning;
    for (auto&& elm : complete_binning)
        if (elm>energy_th)
            new_binning.push_back(elm);
    return new_binning;
}

inline std::shared_ptr<TH1D> shrink_histo(TH1D* histo, const std::vector<float>& binning) {
    std::shared_ptr<TH1D> new_histo = std::make_shared<TH1D>(histo->GetName(), histo->GetTitle(), binning.size()-1, &binning[0]);

    std::vector<double> bin_content, bin_error;

    for (int bidx = 1; bidx <= histo->GetNbinsX(); ++bidx) {
        if (histo->GetBinCenter(bidx) > binning[0]) {
            bin_content.push_back(histo->GetBinContent(bidx));
            bin_error.push_back(histo->GetBinError(bidx));
        }
    }

    for (int bidx = 1; bidx <= new_histo->GetNbinsX(); ++bidx) {
        new_histo->SetBinContent(bidx, bin_content[bidx-1]);
        new_histo->SetBinError(bidx, bin_error[bidx-1]);
    }

    return new_histo;
}

inline std::shared_ptr<TGraphAsymmErrors> fromTH1toGR(std::shared_ptr<TH1> histo) {
    int points = histo->GetNbinsX();
    std::vector<double> x_values, y_values, x_err, y_err;

    for (int bidx = 1; bidx <= points; ++bidx) {
        x_values.push_back(histo->GetBinCenter(bidx));
        y_values.push_back(histo->GetBinContent(bidx));
        x_err.push_back(0);
        y_err.push_back(histo->GetBinError(bidx));
    }

    std::shared_ptr<TGraphAsymmErrors> gr = std::make_shared<TGraphAsymmErrors>(points, &x_values[0], &y_values[0], &x_err[0], &y_err[0]);
    gr->SetName(histo->GetName());
    gr->SetTitle(histo->GetTitle());
    gr->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    gr->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());

    return gr;
}

void buildEfficiency(
    const char* input_file_path, 
    const char* output_file_path,
    const char* energy_config_file,
    const double energy_th = 20) {

        auto energy_binning = shrink_binning(parse_energy_config(energy_config_file), energy_th);
        
        TFile *input_file = TFile::Open(input_file_path, "READ");
        if (!input_file->IsOpen()) {
            std::cerr << "\n\nError opening input ROOT file [" << input_file_path << "]\n\n";
            exit(100);
        }

        /*
        auto h_maxelayer_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_maxelayer_lastcut_pass")), energy_binning);
        auto h_maxelayer_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_maxelayer_lastcut")), energy_binning);
        auto h_maxbarlayer_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_maxbarlayer_lastcut_pass")), energy_binning);
        auto h_maxbarlayer_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_maxbarlayer_lastcut")), energy_binning);
        auto h_bgotrack_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_bgotrack_lastcut_pass")), energy_binning);
        auto h_bgotrack_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_bgotrack_lastcut")), energy_binning);
        auto h_bgofiducial_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_bgofiducial_lastcut_pass")), energy_binning);
        auto h_bgofiducial_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_bgofiducial_lastcut")), energy_binning);
        */
       
        auto h_nbarlayer13_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_nbarlayer13_lastcut_pass")), energy_binning);
        auto h_nbarlayer13_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_nbarlayer13_lastcut")), energy_binning);
        auto h_maxrms_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_maxrms_lastcut_pass")), energy_binning);
        auto h_maxrms_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_maxrms_lastcut")), energy_binning);
        auto h_trackselection_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_trackselection_lastcut_pass")), energy_binning);
        auto h_trackselection_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_trackselection_lastcut")), energy_binning);
        auto h_psdstkmatch_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_psdstkmatch_lastcut_pass")), energy_binning);
        auto h_psdstkmatch_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_psdstkmatch_lastcut")), energy_binning);
        auto h_psdcharge_lastcut_pass = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_psdcharge_lastcut_pass")), energy_binning);
        auto h_psdcharge_lastcut = shrink_histo(static_cast<TH1D*>(input_file->Get("Stats/h_psdcharge_lastcut")), energy_binning);
        
        /*
        h_maxelayer_lastcut_pass->SetDirectory(0);
        h_maxelayer_lastcut->SetDirectory(0);
        h_maxbarlayer_lastcut_pass->SetDirectory(0);
        h_maxbarlayer_lastcut->SetDirectory(0);
        h_bgotrack_lastcut_pass->SetDirectory(0);
        h_bgotrack_lastcut->SetDirectory(0);
        h_bgofiducial_lastcut_pass->SetDirectory(0);
        h_bgofiducial_lastcut->SetDirectory(0);
        */

        h_nbarlayer13_lastcut_pass->SetDirectory(0);
        h_nbarlayer13_lastcut->SetDirectory(0);
        h_maxrms_lastcut_pass->SetDirectory(0);
        h_maxrms_lastcut->SetDirectory(0);
        h_trackselection_lastcut_pass->SetDirectory(0);
        h_trackselection_lastcut->SetDirectory(0);
        h_psdstkmatch_lastcut_pass->SetDirectory(0);
        h_psdstkmatch_lastcut->SetDirectory(0);
        h_psdcharge_lastcut_pass->SetDirectory(0);
        h_psdcharge_lastcut->SetDirectory(0);

        input_file->Close();

        /*
        std::unique_ptr<TEfficiency> maxelayer_eff;
        std::unique_ptr<TEfficiency> maxbarlayer_eff;
        std::unique_ptr<TEfficiency> bgotrack_eff;
        std::unique_ptr<TEfficiency> bgofiducial_eff;
        */
        
        std::unique_ptr<TEfficiency> nbarlayer13_eff;
        std::unique_ptr<TEfficiency> maxrms_eff;
        std::unique_ptr<TEfficiency> trackselection_eff;
        std::unique_ptr<TEfficiency> psdstkmatch_eff;
        std::unique_ptr<TEfficiency> psdcharge_eff;
        std::unique_ptr<TEfficiency> stkcharge_eff;

        /*
        if (TEfficiency::CheckConsistency(*h_maxelayer_lastcut_pass, *h_maxelayer_lastcut))
            maxelayer_eff = std::make_unique<TEfficiency>(*h_maxelayer_lastcut_pass, *h_maxelayer_lastcut);
        if (TEfficiency::CheckConsistency(*h_maxbarlayer_lastcut_pass, *h_maxbarlayer_lastcut))
            maxbarlayer_eff = std::make_unique<TEfficiency>(*h_maxbarlayer_lastcut_pass, *h_maxbarlayer_lastcut);
        if (TEfficiency::CheckConsistency(*h_bgotrack_lastcut_pass, *h_bgotrack_lastcut))
            bgotrack_eff = std::make_unique<TEfficiency>(*h_bgotrack_lastcut_pass, *h_bgotrack_lastcut);
        if (TEfficiency::CheckConsistency(*h_bgofiducial_lastcut_pass, *h_bgofiducial_lastcut))
            bgofiducial_eff = std::make_unique<TEfficiency>(*h_bgofiducial_lastcut_pass, *h_bgofiducial_lastcut);
        */

        if (TEfficiency::CheckConsistency(*h_nbarlayer13_lastcut_pass, *h_nbarlayer13_lastcut))
            nbarlayer13_eff = std::make_unique<TEfficiency>(*h_nbarlayer13_lastcut_pass, *h_nbarlayer13_lastcut);
        if (TEfficiency::CheckConsistency(*h_maxrms_lastcut_pass, *h_maxrms_lastcut))
            maxrms_eff = std::make_unique<TEfficiency>(*h_maxrms_lastcut_pass, *h_maxrms_lastcut);
        if (TEfficiency::CheckConsistency(*h_trackselection_lastcut_pass, *h_trackselection_lastcut))
            trackselection_eff = std::make_unique<TEfficiency>(*h_trackselection_lastcut_pass, *h_trackselection_lastcut);
        if (TEfficiency::CheckConsistency(*h_psdstkmatch_lastcut_pass, *h_psdstkmatch_lastcut))
            psdstkmatch_eff = std::make_unique<TEfficiency>(*h_psdstkmatch_lastcut_pass, *h_psdstkmatch_lastcut);
        if (TEfficiency::CheckConsistency(*h_psdcharge_lastcut_pass, *h_psdcharge_lastcut))
            psdcharge_eff = std::make_unique<TEfficiency>(*h_psdcharge_lastcut_pass, *h_psdcharge_lastcut);

        /*
        maxelayer_eff->SetStatisticOption(TEfficiency::kBUniform);
        maxbarlayer_eff->SetStatisticOption(TEfficiency::kBUniform);
        bgotrack_eff->SetStatisticOption(TEfficiency::kBUniform);
        bgofiducial_eff->SetStatisticOption(TEfficiency::kBUniform);
        */

        nbarlayer13_eff->SetStatisticOption(TEfficiency::kBUniform);
        maxrms_eff->SetStatisticOption(TEfficiency::kBUniform);
        trackselection_eff->SetStatisticOption(TEfficiency::kBUniform);
        psdstkmatch_eff->SetStatisticOption(TEfficiency::kBUniform);
        psdcharge_eff->SetStatisticOption(TEfficiency::kBUniform);

        /*
        maxelayer_eff->SetName("maxelayer_eff");
        maxbarlayer_eff->SetName("maxbarlayer_eff");
        bgotrack_eff->SetName("bgotrack_eff");
        bgofiducial_eff->SetName("bgofiducial_eff");
        */
       
        nbarlayer13_eff->SetName("nbarlayer13_eff");
        maxrms_eff->SetName("maxrms_eff");
        trackselection_eff->SetName("trackselection_eff");
        psdstkmatch_eff->SetName("psdstkmatch_eff");
        psdcharge_eff->SetName("psdcharge_eff");

        /*
        maxelayer_eff->SetTitle("MaxELayer efficiency");
        maxbarlayer_eff->SetTitle("MaxBarLayer efficiency");
        bgotrack_eff->SetTitle("BGO Track efficiency");
        bgofiducial_eff->SetTitle("BGO Fiducial efficiency");
        */
       
        nbarlayer13_eff->SetTitle("nBarLayer13 efficiency");
        maxrms_eff->SetTitle("max RMS efficiency");
        trackselection_eff->SetTitle("Track Selection efficiency");
        psdstkmatch_eff->SetTitle("PSD/STK Match efficiency");
        psdcharge_eff->SetTitle("PSD charge efficiency");

        TFile *output_file = TFile::Open(output_file_path, "RECREATE");
        if (!output_file->IsOpen()) {
            std::cerr << "\n\nError opening output ROOT file [" << output_file_path << "]\n\n";
            exit(100);
        }

        /*
        maxelayer_eff->Write();
        maxbarlayer_eff->Write();
        bgotrack_eff->Write();
        bgofiducial_eff->Write();
        */
       
        nbarlayer13_eff->Write();
        maxrms_eff->Write();
        trackselection_eff->Write();
        psdstkmatch_eff->Write();
        psdcharge_eff->Write();

        output_file->Close();
    }