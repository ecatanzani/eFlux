#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>
#include <string>
#include <tuple>

#include "TPad.h"
#include "TF1.h"
#include "TKey.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

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

inline std::shared_ptr<TGraphErrors> extractGraphFromFile(const char* file, const bool verbose) {
    TFile* infile {TFile::Open(file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key {nullptr};
    std::shared_ptr<TGraphErrors> mygraph;
    while ((key=static_cast<TKey*>(nextkey())))  {
        TObject *obj {key->ReadObj()};
        if (obj->IsA()->InheritsFrom(TGraphErrors::Class())) {
            mygraph = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound TGraph in input file [" << mygraph->GetName() << "] --> [" << file << "]";
    
    return mygraph;
}

inline std::shared_ptr<TGraph> extractGraphNEFromFile(const char* file, const bool verbose) {
    TFile* infile {TFile::Open(file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key {nullptr};
    std::shared_ptr<TGraph> mygraph;
    while ((key=static_cast<TKey*>(nextkey())))  {
        TObject *obj {key->ReadObj()};
        if (obj->IsA()->InheritsFrom(TGraph::Class())) {
            mygraph = std::shared_ptr<TGraph>(static_cast<TGraph*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound TGraph in input file [" << mygraph->GetName() << "] --> [" << file << "]";
    
    return mygraph;
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
        std::cout << "\nFound histo in input file [" << myhisto->GetName() << "] --> [" << file << "]";
    
    return myhisto;
}

inline std::tuple<std::shared_ptr<TGraph>, std::shared_ptr<TGraph>> evaluate_err(std::shared_ptr<TGraphErrors> flux, std::shared_ptr<TH1D> histo) {
    std::vector<double> energy (flux->GetN(), 0);
    std::vector<double> bin (flux->GetN(), 0);
    std::vector<double> stat (flux->GetN(), 0);

    double bin_error {0};
    double bin_content {0};
    for (int pidx=0; pidx<flux->GetN(); ++pidx) {
        bin[pidx] = pidx+1;
        energy[pidx] = flux->GetPointX(pidx);
        bin_error = histo->GetBinError(histo->FindBin(flux->GetPointX(pidx)));
        bin_content = histo->GetBinContent(histo->FindBin(flux->GetPointX(pidx)));
        stat[pidx] = (bin_error/bin_content)*100;
    }

    std::shared_ptr<TGraph> gr_stat_err = std::make_shared<TGraph>(static_cast<int>(energy.size()), energy.data(), stat.data());
    gr_stat_err->SetName("gr_flux_stat_err");
    gr_stat_err->SetTitle("Statistical error on flux");
    gr_stat_err->GetXaxis()->SetTitle("Energy [GeV]");
    gr_stat_err->GetYaxis()->SetTitle("Stat [%]");

    std::shared_ptr<TGraph> gr_bin_stat_err = std::make_shared<TGraph>(static_cast<int>(bin.size()), bin.data(), stat.data());
    return std::make_tuple(gr_stat_err, gr_bin_stat_err);
}

void buildFluxStat(
    const char* flux_file,
    const char* data_selection_file,
    const char* background_estimation_file,
    const char* output_file,
    const bool verbose = true) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;
        
        auto flux = extractGraphFromFile(flux_file, verbose);
        auto data_selection = extractHistoFromFile(data_selection_file, verbose);
        auto background_estimation = extractHistoFromFile(background_estimation_file, verbose);

        data_selection->Sumw2();
        background_estimation->Sumw2();

        data_selection->Add(background_estimation.get(), -1);

        auto flux_gr_stat = evaluate_err(flux, data_selection);
        
        std::unique_ptr<TFile> file = std::make_unique<TFile>(output_file, "RECREATE");
        if (file->IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]" << std::endl;
            exit(100);
        }

        std::get<0>(flux_gr_stat)->SetName("gr_flux_stat_err");
        std::get<0>(flux_gr_stat)->SetTitle("Statistical error on flux");
        std::get<0>(flux_gr_stat)->GetXaxis()->SetTitle("Energy [GeV]");
        std::get<0>(flux_gr_stat)->GetYaxis()->SetTitle("Stat [%]");

        std::get<1>(flux_gr_stat)->SetName("gr_bin_flux_stat_err");
        std::get<1>(flux_gr_stat)->SetTitle("Statistical error on flux");
        std::get<1>(flux_gr_stat)->GetXaxis()->SetTitle("Energy bin");
        std::get<1>(flux_gr_stat)->GetYaxis()->SetTitle("Stat [%]");

        std::get<0>(flux_gr_stat)->Write();
        std::get<1>(flux_gr_stat)->Write();

        file->Close();
    
    }

void buildSignalEfficiencyStat(
    const char* flux_file,
    const char* signal_efficiency_file,
    const char* output_file,
    const bool verbose = true) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;
        
        auto flux = extractGraphFromFile(flux_file, verbose);
        auto signal_efficiency = extractHistoFromFile(signal_efficiency_file, verbose);

        signal_efficiency->Sumw2();
        
        auto signal_efficiency_gr_stat = evaluate_err(flux, signal_efficiency);
        
        std::unique_ptr<TFile> file = std::make_unique<TFile>(output_file, "RECREATE");
        if (file->IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]" << std::endl;
            exit(100);
        }

        std::get<0>(signal_efficiency_gr_stat)->SetName("gr_signal_efficiency_stat_err");
        std::get<0>(signal_efficiency_gr_stat)->SetTitle("Statistical error on signal efficiency");
        std::get<0>(signal_efficiency_gr_stat)->GetXaxis()->SetTitle("Energy [GeV]");
        std::get<0>(signal_efficiency_gr_stat)->GetYaxis()->SetTitle("Stat [%]");

        std::get<1>(signal_efficiency_gr_stat)->SetName("gr_bin_signal_efficiency_stat_err");
        std::get<1>(signal_efficiency_gr_stat)->SetTitle("Statistical error on signal efficiency");
        std::get<1>(signal_efficiency_gr_stat)->GetXaxis()->SetTitle("Energy bin");
        std::get<1>(signal_efficiency_gr_stat)->GetYaxis()->SetTitle("Stat [%]");

        std::get<0>(signal_efficiency_gr_stat)->Write();
        std::get<1>(signal_efficiency_gr_stat)->Write();

        file->Close();
    
    }

void buildAcceptanceStat(
    const char* flux_file,
    const char* acceptance_file,
    const char* output_file,
    const bool verbose = true) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;
        
        auto flux = extractGraphFromFile(flux_file, verbose);
        auto acceptance = extractHistoFromFile(acceptance_file, verbose);

        acceptance->Sumw2();
        
        auto acceptance_gr_stat = evaluate_err(flux, acceptance);
        
        std::unique_ptr<TFile> file = std::make_unique<TFile>(output_file, "RECREATE");
        if (file->IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]" << std::endl;
            exit(100);
        }

        std::get<0>(acceptance_gr_stat)->SetName("gr_acceptance_stat_err");
        std::get<0>(acceptance_gr_stat)->SetTitle("Statistical error on acceptance");
        std::get<0>(acceptance_gr_stat)->GetXaxis()->SetTitle("Energy [GeV]");
        std::get<0>(acceptance_gr_stat)->GetYaxis()->SetTitle("Stat [%]");

        std::get<1>(acceptance_gr_stat)->SetName("gr_bin_acceptance_stat_err");
        std::get<1>(acceptance_gr_stat)->SetTitle("Statistical error on acceptance");
        std::get<1>(acceptance_gr_stat)->GetXaxis()->SetTitle("Energy bin");
        std::get<1>(acceptance_gr_stat)->GetYaxis()->SetTitle("Stat [%]");

        std::get<0>(acceptance_gr_stat)->Write();
        std::get<1>(acceptance_gr_stat)->Write();

        file->Close();
    
    }

std::tuple<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> build_histo_from_tree(const char* input_tree_file, const char* energy_config_file) {
    
    // Extract energy binning from config file
    auto energy_binning = parse_energy_config(energy_config_file);
    auto energy_nbins = (int)energy_binning.size() - 1;

    // Build histo
    std::shared_ptr<TH1D> histo = std::make_shared<TH1D>("bdt_best_point_syst", "bdt_best_point_syst", energy_nbins, energy_binning.data());
    std::shared_ptr<TH1D> histo_err = std::make_shared<TH1D>("bdt_best_point_syst_err", "bdt_best_point_syst_err", energy_nbins, energy_binning.data());

    std::unique_ptr<TFile> bdt_best_point_file = std::make_unique<TFile>(input_tree_file, "READ");
    if (bdt_best_point_file->IsZombie()) {
        std::cerr << "\n\nError reading input ROOT file [" << input_tree_file << "]" << std::endl;
        exit(100);
    }

    TTreeReader myReader("best_cuts_tree", bdt_best_point_file.get());
    
    TTreeReaderValue<double> f_ec_b_sub_Y(myReader, "f_ec_b_sub_Y");
    TTreeReaderValue<double> f_ec_b_sub_err(myReader, "f_ec_b_sub_err");
    TTreeReaderValue<int> energy_bin(myReader, "energy_bin");
    
    
    while (myReader.Next()) {
        histo->SetBinContent(*energy_bin, *f_ec_b_sub_Y);
        histo->SetBinError(*energy_bin, *f_ec_b_sub_err);

        histo_err->SetBinContent(*energy_bin, 0);
        histo_err->SetBinError(*energy_bin, *f_ec_b_sub_err);
    }

    bdt_best_point_file->Close();

    return std::make_tuple(histo, histo_err);
}

void buildBDTBestPointSyst(
    const char* flux_file,
    const char* bdt_best_point_tree_file,
    const char* energy_config_file,
    const char* output_file,
    const bool verbose = true) {

        if (verbose)
            std::cout << "\n\nReading input files... " << std::endl;
        
        auto flux = extractGraphFromFile(flux_file, verbose);

        std::shared_ptr<TH1D> bdt_syst, bdt_syst_err;
        std::tie(bdt_syst, bdt_syst_err) = build_histo_from_tree(bdt_best_point_tree_file, energy_config_file);

        auto bdt_syst_gr_syst = evaluate_err(flux, bdt_syst);
        
        std::unique_ptr<TFile> file = std::make_unique<TFile>(output_file, "RECREATE");
        if (file->IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]" << std::endl;
            exit(100);
        }

        std::get<0>(bdt_syst_gr_syst)->SetName("gr_bdt_bestopint_syst_err");
        std::get<0>(bdt_syst_gr_syst)->SetTitle("Systematic error on BDT best point");
        std::get<0>(bdt_syst_gr_syst)->GetXaxis()->SetTitle("Energy [GeV]");
        std::get<0>(bdt_syst_gr_syst)->GetYaxis()->SetTitle("Stat [%]");

        std::get<1>(bdt_syst_gr_syst)->SetName("gr_bin_bdt_bestopint_syst_err");
        std::get<1>(bdt_syst_gr_syst)->SetTitle("Systematic error on BDT best point");
        std::get<1>(bdt_syst_gr_syst)->GetXaxis()->SetTitle("Energy bin");
        std::get<1>(bdt_syst_gr_syst)->GetYaxis()->SetTitle("Stat [%]");

        std::get<0>(bdt_syst_gr_syst)->Write();
        std::get<1>(bdt_syst_gr_syst)->Write();

        bdt_syst->Write();
        bdt_syst_err->Write();

        file->Close();
    
    }

void buildSystContributionsPlot(
    const char* syst_signal_efficiency_file,
    const char* syst_acceptance_file,
    const char* syst_bdt_best_point_file,
    const char* syst_preselection_eff_file,
    const char* stat_flux_file,
    const char* output_file = "err_comparison.root",
    const bool verbose = true) {

        auto root_sum_square = [](
            std::shared_ptr<TGraph> g1, 
            std::shared_ptr<TGraph> g2, 
            std::shared_ptr<TGraph> g3, 
            std::shared_ptr<TGraph> g4) -> std::shared_ptr<TGraph> {
                    
                std::vector<double> energy (g1->GetN(), 0);
                std::vector<double> syst (g1->GetN(), 0);

                double syst_1 {0};
                double syst_2 {0};
                double syst_3 {0};
                double syst_4 {0};
                double syst_tot {0};

                for (int pidx=0; pidx<g1->GetN(); ++pidx) {
                    energy[pidx] = g1->GetPointX(pidx);
                    
                    syst_1 = g1->GetPointY(pidx);
                    syst_2 = g2->GetPointY(pidx);
                    syst_3 = g3->GetPointY(pidx);
                    syst_4 = g4->GetPointY(pidx);
                    syst_tot = sqrt(pow(syst_1, 2) + pow(syst_2, 2) + pow(syst_3, 2) + pow(syst_4, 2));
                    syst[pidx] = syst_tot;
                }

                std::shared_ptr<TGraph> gr_syst = std::make_shared<TGraph>(energy.size(), energy.data(), syst.data());
                gr_syst->SetName("gr_syst_total");
                gr_syst->SetTitle("Total systematic error");
                gr_syst->GetXaxis()->SetTitle("Energy [GeV]");
                gr_syst->GetYaxis()->SetTitle("syst [%]");

                return gr_syst;
            };

        auto gr_signal_efficiency = extractGraphNEFromFile(syst_signal_efficiency_file, verbose);
        auto gr_acceptance = extractGraphNEFromFile(syst_acceptance_file, verbose);
        auto gr_bdt_best_point = extractGraphNEFromFile(syst_bdt_best_point_file, verbose);
        auto gr_preselection_eff = extractGraphNEFromFile(syst_preselection_eff_file, verbose);
        auto gr_flux_stat = extractGraphNEFromFile(stat_flux_file, verbose);

        gr_signal_efficiency->SetTitle("signal efficiency");
        gr_acceptance->SetTitle("acceptance");
        gr_bdt_best_point->SetTitle("BDT best point");
        gr_preselection_eff->SetTitle("preselection efficiency");
        
        gr_flux_stat->SetTitle("Statistic error");
        gr_flux_stat->GetXaxis()->SetTitle("Energy [GeV]");
        gr_flux_stat->GetYaxis()->SetTitle("uncertainty contribution [%]");

        auto gr_syst_total = root_sum_square(gr_signal_efficiency, gr_acceptance, gr_bdt_best_point, gr_preselection_eff);

        TCanvas c_syst("c_syst", "c_syst", 500, 500);
        c_syst.SetTicks();

        gr_signal_efficiency->SetLineWidth(2);
        gr_acceptance->SetLineWidth(2);
        gr_bdt_best_point->SetLineWidth(2);
        gr_preselection_eff->SetLineWidth(2);

        gr_signal_efficiency->SetMarkerStyle(20);
        gr_acceptance->SetMarkerStyle(20);
        gr_bdt_best_point->SetMarkerStyle(20);
        gr_preselection_eff->SetMarkerStyle(20);

        gr_signal_efficiency->SetLineColor(kCyan+1);
        gr_acceptance->SetLineColor(kGreen+1);
        gr_bdt_best_point->SetLineColor(kMagenta+1);
        gr_preselection_eff->SetLineColor(kBlue+1);

        gr_signal_efficiency->SetMarkerColor(kCyan+1);
        gr_acceptance->SetMarkerColor(kGreen+1);
        gr_bdt_best_point->SetMarkerColor(kMagenta+1);
        gr_preselection_eff->SetMarkerColor(kBlue+1);

        gr_bdt_best_point->GetXaxis()->SetTitle("Energy [GeV]");
        gr_bdt_best_point->GetYaxis()->SetTitle("syst [%]");

        gr_bdt_best_point->Draw("ALP");
        gr_preselection_eff->Draw("LP");
        gr_signal_efficiency->Draw("LP");
        gr_acceptance->Draw("LP");

        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGrid(1,1);

        gPad->Update(); 
        gr_bdt_best_point->SetMinimum(1e-3);
        gr_bdt_best_point->SetMaximum(1e+2); 
        gPad->Update();

        auto legend = c_syst.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("lp");
        }

        TCanvas c_syst_comparison("c_syst_comparison", "c_syst_comparison", 500, 500);
        c_syst_comparison.SetTicks();

        gr_syst_total->SetLineWidth(2);
        gr_syst_total->SetMarkerStyle(20);

        gr_flux_stat->SetLineWidth(2);
        gr_flux_stat->SetMarkerStyle(20);

        gr_syst_total->SetLineColor(kBlack);
        gr_syst_total->SetMarkerColor(kBlack);

        gr_flux_stat->SetLineColor(kRed);
        gr_flux_stat->SetMarkerColor(kRed);

        gr_flux_stat->Draw("ALP");
        gr_syst_total->Draw("LP");

        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGrid(1,1);

        gPad->Update(); 
        gr_flux_stat->SetMinimum(1e-2);
        gr_flux_stat->SetMaximum(1e+2); 
        gPad->Update();

        legend = c_syst_comparison.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("lp");
        }


        TFile file(output_file, "RECREATE");
        if (!file.IsOpen())
        {
            std::cout << "Error: cannot write output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        c_syst.Write();
        c_syst_comparison.Write();

        file.Close();

    }