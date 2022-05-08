#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include "TF1.h"
#include "TH1D.h"
#include "TPDF.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TPaveLabel.h"
#include "TEfficiency.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"

void buildEfficiencyRatios(
    const char* data_eff_file,
    const char* electron_mc_eff_file,
    const char* output_file = "efficiency_ratios.root",
    const double fit_low_energy_gev = 15,
    const double fit_high_energy_gev = 800,
    const bool verbose = true)
    {
        TFile* datafile = TFile::Open(data_eff_file, "READ");
        if (datafile->IsZombie()) {
            std::cerr << "\n\nError opening input DATA ROOT file [" << data_eff_file << "]\n\n";
            exit(100);
        }

        auto data_trigger_eff_het_over_let_bdt                      = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_let_bdt"));
        auto data_trigger_eff_het_over_unb_bdt                      = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_unb_bdt"));
        auto data_maxrms_eff_bdt                                    = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_eff_bdt"));
        auto data_nbarlayer13_eff_bdt                               = static_cast<TEfficiency*>(datafile->Get("efficiencies/nbarlayer13_eff_bdt"));
        auto data_maxrms_and_nbarlayer13_eff_bdt                    = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_and_nbarlayer13_eff_bdt"));
        auto data_sumrms_low_energy_eff_bdt                         = static_cast<TEfficiency*>(datafile->Get("efficiencies/sumrms_low_energy_eff_bdt"));
        auto data_track_selection_eff_bdt                           = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_bdt"));
        auto data_track_selection_eff_within_stk_fvolume_bdt        = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_within_stk_fvolume_bdt"));
        auto data_stk_1_rm_eff_bdt                                  = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_1_rm_eff_bdt"));
        auto data_psd_stk_match_eff_bdt                             = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_stk_match_eff_bdt"));
        auto data_psd_charge_eff_bdt                                = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_charge_eff_bdt"));
        auto data_stk_charge_eff_bdt                                = static_cast<TEfficiency*>(datafile->Get("efficiencies/stk_charge_eff_bdt"));

        data_trigger_eff_het_over_let_bdt                       ->SetDirectory(0);
        data_trigger_eff_het_over_unb_bdt                       ->SetDirectory(0);
        data_maxrms_eff_bdt                                     ->SetDirectory(0);
        data_nbarlayer13_eff_bdt                                ->SetDirectory(0);
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetDirectory(0);
        data_sumrms_low_energy_eff_bdt                          ->SetDirectory(0);
        data_track_selection_eff_bdt                            ->SetDirectory(0);
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetDirectory(0);
        data_stk_1_rm_eff_bdt                                   ->SetDirectory(0);
        data_psd_stk_match_eff_bdt                              ->SetDirectory(0);
        data_psd_charge_eff_bdt                                 ->SetDirectory(0);
        data_stk_charge_eff_bdt                                 ->SetDirectory(0);

        datafile->Close();

        TFile* mcfile = TFile::Open(electron_mc_eff_file, "READ");
        if (mcfile->IsZombie()) {
            std::cerr << "\n\nError opening input MC ROOT file [" << electron_mc_eff_file << "]\n\n";
            exit(100);
        }

        auto mc_trigger_eff_het_over_let_bdt                        = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_let_bdt"));
        auto mc_trigger_eff_het_over_unb_bdt                        = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_unb_bdt"));
        auto mc_maxrms_eff_bdt                                      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_eff_bdt"));
        auto mc_nbarlayer13_eff_bdt                                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/nbarlayer13_eff_bdt"));
        auto mc_maxrms_and_nbarlayer13_eff_bdt                      = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_and_nbarlayer13_eff_bdt"));
        auto mc_sumrms_low_energy_eff_bdt                           = static_cast<TEfficiency*>(mcfile->Get("efficiencies/sumrms_low_energy_eff_bdt"));
        auto mc_track_selection_eff_bdt                             = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_bdt"));
        auto mc_track_selection_eff_within_stk_fvolume_bdt          = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_within_stk_fvolume_bdt"));
        auto mc_stk_1_rm_eff_bdt                                    = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_1_rm_eff_bdt"));
        auto mc_psd_stk_match_eff_bdt                               = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_stk_match_eff_bdt"));
        auto mc_psd_charge_eff_bdt                                  = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_charge_eff_bdt"));
        auto mc_stk_charge_eff_bdt                                  = static_cast<TEfficiency*>(mcfile->Get("efficiencies/stk_charge_eff_bdt"));

        mc_trigger_eff_het_over_let_bdt                         ->SetDirectory(0);
        mc_trigger_eff_het_over_unb_bdt                         ->SetDirectory(0);
        mc_maxrms_eff_bdt                                       ->SetDirectory(0);
        mc_nbarlayer13_eff_bdt                                  ->SetDirectory(0);
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetDirectory(0);
        mc_sumrms_low_energy_eff_bdt                            ->SetDirectory(0);
        mc_track_selection_eff_bdt                              ->SetDirectory(0);
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetDirectory(0);
        mc_stk_1_rm_eff_bdt                                     ->SetDirectory(0);
        mc_psd_stk_match_eff_bdt                                ->SetDirectory(0);
        mc_psd_charge_eff_bdt                                   ->SetDirectory(0);
        mc_stk_charge_eff_bdt                                   ->SetDirectory(0);

        mcfile->Close();

        // Convert efficiencies from TEfficiency to TH1D
        auto get_efficiency_histo = [](TEfficiency* eff) -> std::shared_ptr<TH1D> {
            auto h_eff = std::shared_ptr<TH1D>(static_cast<TH1D*>(eff->GetPassedHistogram()->Clone()));
            auto h_total = std::shared_ptr<TH1D>(static_cast<TH1D*>(eff->GetTotalHistogram()->Clone()));
            h_eff->Divide(h_total.get());
            return h_eff;
        };

        auto get_efficiency_histo_he = [](TEfficiency* eff) -> std::shared_ptr<TH1D> {
            auto gr = std::make_shared<TGraphAsymmErrors>(*static_cast<TGraphAsymmErrors*>(eff->CreateGraph()));
            auto h_eff = std::shared_ptr<TH1D>(static_cast<TH1D*>(eff->GetPassedHistogram()->Clone()));
            h_eff->Reset();
            for (int i = 0; i < gr->GetN(); ++i) {
                h_eff->SetBinContent(i+1, gr->GetY()[i]);
                h_eff->SetBinError(i+1, std::max(gr->GetEYhigh()[i], gr->GetEYlow()[i]));
            }
            return h_eff;
        };

        auto h_data_trigger_eff_het_over_let_bdt = get_efficiency_histo_he(data_trigger_eff_het_over_let_bdt);
        auto h_data_trigger_eff_het_over_unb_bdt = get_efficiency_histo_he(data_trigger_eff_het_over_unb_bdt);
        auto h_data_maxrms_eff_bdt = get_efficiency_histo_he(data_maxrms_eff_bdt);
        auto h_data_nbarlayer13_eff_bdt = get_efficiency_histo_he(data_nbarlayer13_eff_bdt);
        auto h_data_maxrms_and_nbarlayer13_eff_bdt = get_efficiency_histo_he(data_maxrms_and_nbarlayer13_eff_bdt);
        auto h_data_sumrms_low_energy_eff_bdt = get_efficiency_histo_he(data_sumrms_low_energy_eff_bdt);
        auto h_data_track_selection_eff_bdt = get_efficiency_histo_he(data_track_selection_eff_bdt);
        auto h_data_track_selection_eff_within_stk_fvolume_bdt = get_efficiency_histo_he(data_track_selection_eff_within_stk_fvolume_bdt);
        auto h_data_stk_1_rm_eff_bdt = get_efficiency_histo_he(data_stk_1_rm_eff_bdt);
        auto h_data_psd_stk_match_eff_bdt = get_efficiency_histo_he(data_psd_stk_match_eff_bdt);
        auto h_data_psd_charge_eff_bdt = get_efficiency_histo_he(data_psd_charge_eff_bdt);
        auto h_data_stk_charge_eff_bdt = get_efficiency_histo_he(data_stk_charge_eff_bdt);

        auto h_mc_trigger_eff_het_over_let_bdt = get_efficiency_histo_he(mc_trigger_eff_het_over_let_bdt);
        auto h_mc_trigger_eff_het_over_unb_bdt = get_efficiency_histo_he(mc_trigger_eff_het_over_unb_bdt);
        auto h_mc_maxrms_eff_bdt = get_efficiency_histo_he(mc_maxrms_eff_bdt);
        auto h_mc_nbarlayer13_eff_bdt = get_efficiency_histo_he(mc_nbarlayer13_eff_bdt);
        auto h_mc_maxrms_and_nbarlayer13_eff_bdt = get_efficiency_histo_he(mc_maxrms_and_nbarlayer13_eff_bdt);
        auto h_mc_sumrms_low_energy_eff_bdt = get_efficiency_histo_he(mc_sumrms_low_energy_eff_bdt);
        auto h_mc_track_selection_eff_bdt = get_efficiency_histo_he(mc_track_selection_eff_bdt);
        auto h_mc_track_selection_eff_within_stk_fvolume_bdt = get_efficiency_histo_he(mc_track_selection_eff_within_stk_fvolume_bdt);
        auto h_mc_stk_1_rm_eff_bdt = get_efficiency_histo_he(mc_stk_1_rm_eff_bdt);
        auto h_mc_psd_stk_match_eff_bdt = get_efficiency_histo_he(mc_psd_stk_match_eff_bdt);
        auto h_mc_psd_charge_eff_bdt = get_efficiency_histo_he(mc_psd_charge_eff_bdt);
        auto h_mc_stk_charge_eff_bdt = get_efficiency_histo_he(mc_stk_charge_eff_bdt);

        // Clone data histos for the future ratio
        auto h_ratio_trigger_eff_het_over_let_bdt = static_cast<TH1D*>(h_data_trigger_eff_het_over_let_bdt->Clone("h_ratio_trigger_eff_het_over_let_bdt"));
        auto h_ratio_trigger_eff_het_over_unb_bdt = static_cast<TH1D*>(h_data_trigger_eff_het_over_unb_bdt->Clone("h_ratio_trigger_eff_het_over_unb_bdt"));
        auto h_ratio_maxrms_eff_bdt = static_cast<TH1D*>(h_data_maxrms_eff_bdt->Clone("h_ratio_maxrms_eff_bdt"));
        auto h_ratio_nbarlayer13_eff_bdt = static_cast<TH1D*>(h_data_nbarlayer13_eff_bdt->Clone("h_ratio_nbarlayer13_eff_bdt"));
        auto h_ratio_maxrms_and_nbarlayer13_eff_bdt = static_cast<TH1D*>(h_data_maxrms_and_nbarlayer13_eff_bdt->Clone("h_ratio_maxrms_and_nbarlayer13_eff_bdt"));
        auto h_ratio_sumrms_low_energy_eff_bdt = static_cast<TH1D*>(h_data_sumrms_low_energy_eff_bdt->Clone("h_ratio_sumrms_low_energy_eff_bdt"));
        auto h_ratio_track_selection_eff_bdt = static_cast<TH1D*>(h_data_track_selection_eff_bdt->Clone("h_ratio_track_selection_eff_bdt"));
        auto h_ratio_track_selection_eff_within_stk_fvolume_bdt = static_cast<TH1D*>(h_data_track_selection_eff_within_stk_fvolume_bdt->Clone("h_ratio_track_selection_eff_within_stk_fvolume_bdt"));
        auto h_ratio_stk_1_rm_eff_bdt = static_cast<TH1D*>(h_data_stk_1_rm_eff_bdt->Clone("h_ratio_stk_1_rm_eff_bdt"));
        auto h_ratio_psd_stk_match_eff_bdt = static_cast<TH1D*>(h_data_psd_stk_match_eff_bdt->Clone("h_ratio_psd_stk_match_eff_bdt"));
        auto h_ratio_psd_charge_eff_bdt = static_cast<TH1D*>(h_data_psd_charge_eff_bdt->Clone("h_ratio_psd_charge_eff_bdt"));
        auto h_ratio_stk_charge_eff_bdt = static_cast<TH1D*>(h_data_stk_charge_eff_bdt->Clone("h_ratio_stk_charge_eff_bdt"));

        // Divide data over mc efficiencies
        h_ratio_trigger_eff_het_over_let_bdt->Divide(h_mc_trigger_eff_het_over_let_bdt.get());
        h_ratio_trigger_eff_het_over_unb_bdt->Divide(h_mc_trigger_eff_het_over_unb_bdt.get());
        h_ratio_maxrms_eff_bdt->Divide(h_mc_maxrms_eff_bdt.get());
        h_ratio_nbarlayer13_eff_bdt->Divide(h_mc_nbarlayer13_eff_bdt.get());
        h_ratio_maxrms_and_nbarlayer13_eff_bdt->Divide(h_mc_maxrms_and_nbarlayer13_eff_bdt.get());
        h_ratio_sumrms_low_energy_eff_bdt->Divide(h_mc_sumrms_low_energy_eff_bdt.get());
        h_ratio_track_selection_eff_bdt->Divide(h_mc_track_selection_eff_bdt.get());
        h_ratio_track_selection_eff_within_stk_fvolume_bdt->Divide(h_mc_track_selection_eff_within_stk_fvolume_bdt.get());
        h_ratio_stk_1_rm_eff_bdt->Divide(h_mc_stk_1_rm_eff_bdt.get());
        h_ratio_psd_stk_match_eff_bdt->Divide(h_mc_psd_stk_match_eff_bdt.get());
        h_ratio_psd_charge_eff_bdt->Divide(h_mc_psd_charge_eff_bdt.get());
        h_ratio_stk_charge_eff_bdt->Divide(h_mc_stk_charge_eff_bdt.get());

        // Fit histos with a pol1
        TF1 f_ratio_trigger_eff_het_over_let_bdt("f_ratio_trigger_eff_het_over_let_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_trigger_eff_het_over_unb_bdt("f_ratio_trigger_eff_het_over_unb_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_maxrms_eff_bdt("f_ratio_maxrms_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_nbarlayer13_eff_bdt("f_ratio_nbarlayer13_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_maxrms_and_nbarlayer13_eff_bdt("f_ratio_maxrms_and_nbarlayer13_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_sumrms_low_energy_eff_bdt("f_ratio_sumrms_low_energy_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        //TF1 f_ratio_track_selection_eff_bdt("f_ratio_track_selection_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_track_selection_eff_bdt("f_ratio_track_selection_eff_bdt", "[0]", 10, 1e+4);
        TF1 f_ratio_track_selection_eff_within_stk_fvolume_bdt("f_ratio_track_selection_eff_within_stk_fvolume_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_stk_1_rm_eff_bdt("f_ratio_stk_1_rm_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_psd_stk_match_eff_bdt("f_ratio_psd_stk_match_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_psd_charge_eff_bdt("f_ratio_psd_charge_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_stk_charge_eff_bdt("f_ratio_stk_charge_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);

        h_ratio_trigger_eff_het_over_let_bdt->Fit(&f_ratio_trigger_eff_het_over_let_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_trigger_eff_het_over_unb_bdt->Fit(&f_ratio_trigger_eff_het_over_unb_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_maxrms_eff_bdt->Fit(&f_ratio_maxrms_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_nbarlayer13_eff_bdt->Fit(&f_ratio_nbarlayer13_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_maxrms_and_nbarlayer13_eff_bdt->Fit(&f_ratio_maxrms_and_nbarlayer13_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_sumrms_low_energy_eff_bdt->Fit(&f_ratio_sumrms_low_energy_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_track_selection_eff_bdt->Fit(&f_ratio_track_selection_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_track_selection_eff_within_stk_fvolume_bdt->Fit(&f_ratio_track_selection_eff_within_stk_fvolume_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_stk_1_rm_eff_bdt->Fit(&f_ratio_stk_1_rm_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_psd_stk_match_eff_bdt->Fit(&f_ratio_psd_stk_match_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_psd_charge_eff_bdt->Fit(&f_ratio_psd_charge_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        h_ratio_stk_charge_eff_bdt->Fit(&f_ratio_stk_charge_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);

        // Save results to the output file
        TFile output_eff_ratio_file(output_file, "RECREATE");
        if (output_eff_ratio_file.IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        h_ratio_trigger_eff_het_over_let_bdt->Write();
        h_ratio_trigger_eff_het_over_unb_bdt->Write();
        h_ratio_maxrms_eff_bdt->Write();
        h_ratio_nbarlayer13_eff_bdt->Write();
        h_ratio_maxrms_and_nbarlayer13_eff_bdt->Write();
        h_ratio_sumrms_low_energy_eff_bdt->Write();
        h_ratio_track_selection_eff_bdt->Write();
        h_ratio_track_selection_eff_within_stk_fvolume_bdt->Write();
        h_ratio_stk_1_rm_eff_bdt->Write();
        h_ratio_psd_stk_match_eff_bdt->Write();
        h_ratio_psd_charge_eff_bdt->Write();
        h_ratio_stk_charge_eff_bdt->Write();

        f_ratio_trigger_eff_het_over_let_bdt.Write();
        f_ratio_trigger_eff_het_over_unb_bdt.Write();
        f_ratio_maxrms_eff_bdt.Write();
        f_ratio_nbarlayer13_eff_bdt.Write();
        f_ratio_maxrms_and_nbarlayer13_eff_bdt.Write();
        f_ratio_sumrms_low_energy_eff_bdt.Write();
        f_ratio_track_selection_eff_bdt.Write();
        f_ratio_track_selection_eff_within_stk_fvolume_bdt.Write();
        f_ratio_stk_1_rm_eff_bdt.Write();
        f_ratio_psd_stk_match_eff_bdt.Write();
        f_ratio_psd_charge_eff_bdt.Write();
        f_ratio_stk_charge_eff_bdt.Write();

        output_eff_ratio_file.Close();

        // Transform TH1Ds in TGraphs
        auto build_graph = [](TH1D* h) -> std::shared_ptr<TGraphErrors> {

            std::vector<double> x_values (h->GetNbinsX());
            std::vector<double> x_errors (h->GetNbinsX(), 0);
            std::vector<double> y_values (h->GetNbinsX());
            std::vector<double> y_errors (h->GetNbinsX());

            for (int bidx=1; bidx<=h->GetNbinsX(); ++bidx) {
                x_values[bidx-1] = h->GetBinCenter(bidx);
                y_values[bidx-1] = h->GetBinContent(bidx);
                y_errors[bidx-1] = h->GetBinError(bidx);
            }

            std::shared_ptr<TGraphErrors> gr = std::make_shared<TGraphErrors>(h->GetNbinsX(), x_values.data(), y_values.data(), x_errors.data(), y_errors.data());
            
            gr->SetLineColor(kBlue+2);
            gr->SetMarkerColor(kBlue+2);
            gr->SetMarkerStyle(kFullDotMedium);

            gr->GetXaxis()->SetTitle("Energy [GeV]");
            gr->GetYaxis()->SetTitle("Efficiency ratio data/mc");

            return gr;
        };

        // Print summary PDF
        TCanvas print_canvas("print_canvas", "print_canvas");
        print_canvas.SetTicks();

        // Trigger Efficiency - HET over LET
        auto gr_ratio_trigger_eff_het_over_let_bdt = build_graph(h_ratio_trigger_eff_het_over_let_bdt);

        gr_ratio_trigger_eff_het_over_let_bdt->Draw("AP");
        
        f_ratio_trigger_eff_het_over_let_bdt.SetLineWidth(2);
        f_ratio_trigger_eff_het_over_let_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_trigger_eff_het_over_let_bdt->SetMinimum(0);
        gr_ratio_trigger_eff_het_over_let_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label(0.0, 0.95, 0.3, 1, "Trigger Efficiency - HET over LET", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf(","Title:Trigger Efficiency - HET over LET");

        // Trigger Efficiency - HET over UNB
        auto gr_ratio_trigger_eff_het_over_unb_bdt = build_graph(h_ratio_trigger_eff_het_over_unb_bdt);

        gr_ratio_trigger_eff_het_over_unb_bdt->Draw("AP");
        
        f_ratio_trigger_eff_het_over_unb_bdt.SetLineWidth(2);
        f_ratio_trigger_eff_het_over_unb_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_trigger_eff_het_over_unb_bdt->SetMinimum(0);
        gr_ratio_trigger_eff_het_over_unb_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "Trigger Efficiency - HET over UNB", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:Trigger Efficiency - HET over UNB");

        // maxRMS Efficiency
        auto gr_ratio_maxrms_eff_bdt = build_graph(h_ratio_maxrms_eff_bdt);

        gr_ratio_maxrms_eff_bdt->Draw("AP");
        
        f_ratio_maxrms_eff_bdt.SetLineWidth(2);
        f_ratio_maxrms_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_maxrms_eff_bdt->SetMinimum(0);
        gr_ratio_maxrms_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "maxRMS Efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:maxRMS Efficiency");

        // nbarlayer13 Efficiency
        auto gr_ratio_nbarlayer13_eff_bdt = build_graph(h_ratio_nbarlayer13_eff_bdt);

        gr_ratio_nbarlayer13_eff_bdt->Draw("AP");

        f_ratio_nbarlayer13_eff_bdt.SetLineWidth(2);
        f_ratio_nbarlayer13_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_nbarlayer13_eff_bdt->SetMinimum(0);
        gr_ratio_nbarlayer13_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "nbarlayer13 Efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:nbarlayer13 Efficiency");

        // maxrms and nbarlayer13 efficiency
        auto gr_ratio_maxrms_and_nbarlayer13_eff_bdt = build_graph(h_ratio_maxrms_and_nbarlayer13_eff_bdt);

        gr_ratio_maxrms_and_nbarlayer13_eff_bdt->Draw("AP");

        f_ratio_maxrms_and_nbarlayer13_eff_bdt.SetLineWidth(2);
        f_ratio_maxrms_and_nbarlayer13_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_maxrms_and_nbarlayer13_eff_bdt->SetMinimum(0);
        gr_ratio_maxrms_and_nbarlayer13_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "maxrms and nbarlayer13 efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:maxrms and nbarlayer13 efficiency");

        // sumrms low energy efficiency
        auto gr_ratio_sumrms_low_energy_eff_bdt = build_graph(h_ratio_sumrms_low_energy_eff_bdt);

        gr_ratio_sumrms_low_energy_eff_bdt->Draw("AP");

        f_ratio_sumrms_low_energy_eff_bdt.SetLineWidth(2);
        f_ratio_sumrms_low_energy_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_sumrms_low_energy_eff_bdt->SetMinimum(0);
        gr_ratio_sumrms_low_energy_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "sumrms low energy efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:sumrms low energy efficiency");      

        // trackselection efficiency
        auto gr_ratio_track_selection_eff_bdt = build_graph(h_ratio_track_selection_eff_bdt);

        gr_ratio_track_selection_eff_bdt->Draw("AP");

        f_ratio_track_selection_eff_bdt.SetLineWidth(2);
        f_ratio_track_selection_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_track_selection_eff_bdt->SetMinimum(0);
        gr_ratio_track_selection_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "trackselection efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:trackselection efficiency");   

        // trackselection efficiency within STK fiducial volume
        auto gr_ratio_track_selection_eff_within_stk_fvolume_bdt = build_graph(h_ratio_track_selection_eff_within_stk_fvolume_bdt);

        gr_ratio_track_selection_eff_within_stk_fvolume_bdt->Draw("AP");

        f_ratio_track_selection_eff_within_stk_fvolume_bdt.SetLineWidth(2);
        f_ratio_track_selection_eff_within_stk_fvolume_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_track_selection_eff_within_stk_fvolume_bdt->SetMinimum(0);
        gr_ratio_track_selection_eff_within_stk_fvolume_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "trackselection efficiency within STK fiducial volume", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:trackselection efficiency within STK fiducial volume");

        // STK 1 RM efficiency
        auto gr_ratio_stk_1_rm_eff_bdt = build_graph(h_ratio_stk_1_rm_eff_bdt);

        gr_ratio_stk_1_rm_eff_bdt->Draw("AP");

        f_ratio_stk_1_rm_eff_bdt.SetLineWidth(2);
        f_ratio_stk_1_rm_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_stk_1_rm_eff_bdt->SetMinimum(0);
        gr_ratio_stk_1_rm_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK 1 RM efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:STK 1 RM efficiency");

        // psd-stk match efficiency
        auto gr_ratio_psd_stk_match_eff_bdt = build_graph(h_ratio_psd_stk_match_eff_bdt);

        gr_ratio_psd_stk_match_eff_bdt->Draw("AP");

        f_ratio_psd_stk_match_eff_bdt.SetLineWidth(2);
        f_ratio_psd_stk_match_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_psd_stk_match_eff_bdt->SetMinimum(0);
        gr_ratio_psd_stk_match_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "psd-stk match efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:psd-stk match efficiency");

        // psd charge efficiency
        auto gr_ratio_psd_charge_eff_bdt = build_graph(h_ratio_psd_charge_eff_bdt);

        gr_ratio_psd_charge_eff_bdt->Draw("AP");

        f_ratio_psd_charge_eff_bdt.SetLineWidth(2);
        f_ratio_psd_charge_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_psd_charge_eff_bdt->SetMinimum(0);
        gr_ratio_psd_charge_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "psd charge efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf","Title:psd charge efficiency");

        // stk charge efficiency
        auto gr_ratio_stk_charge_eff_bdt = build_graph(h_ratio_stk_charge_eff_bdt);

        gr_ratio_stk_charge_eff_bdt->Draw("AP");

        f_ratio_stk_charge_eff_bdt.SetLineWidth(2);
        f_ratio_stk_charge_eff_bdt.Draw("same");

        gPad->Update(); 
        gr_ratio_stk_charge_eff_bdt->SetMinimum(0);
        gr_ratio_stk_charge_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        label = TPaveLabel(0.0, 0.95, 0.3, 1, "stk charge efficiency", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("efficiency_ratios.pdf)","Title:stk charge efficiency");

    }