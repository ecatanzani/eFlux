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

        auto get_efficiency_ratio = [](TEfficiency* eff_data, TEfficiency* eff_mc) -> std::shared_ptr<TGraphAsymmErrors> {
            std::vector<double> energy;
            std::vector<double> mode_ratio;
            std::vector<double> up_err_effratio;
            std::vector<double> down_err_effratio;
            std::vector<double> up_err_energy;
            std::vector<double> down_err_energy;

            std::unique_ptr<TGraphAsymmErrors> gr_eff_data = std::unique_ptr<TGraphAsymmErrors>(static_cast<TGraphAsymmErrors*>(eff_data->CreateGraph()));
            std::unique_ptr<TGraphAsymmErrors> gr_eff_mc = std::unique_ptr<TGraphAsymmErrors>(static_cast<TGraphAsymmErrors*>(eff_mc->CreateGraph()));

            double up_err {0};
            double err_down {0};

            for (int pidx=0; pidx<gr_eff_data->GetN(); ++pidx) {
                energy.push_back(gr_eff_data->GetPointX(pidx));
                mode_ratio.push_back(gr_eff_data->GetPointY(pidx)/gr_eff_mc->GetPointY(pidx));

                up_err = ((gr_eff_data->GetPointY(pidx) + gr_eff_data->GetErrorYhigh(pidx))/(gr_eff_mc->GetPointY(pidx) - gr_eff_mc->GetErrorYlow(pidx))) - (gr_eff_data->GetPointY(pidx)/gr_eff_mc->GetPointY(pidx));
                    
                err_down = (gr_eff_data->GetPointY(pidx)/gr_eff_mc->GetPointY(pidx)) - ((gr_eff_data->GetPointY(pidx) - gr_eff_data->GetErrorYlow(pidx))/(gr_eff_mc->GetPointY(pidx) + gr_eff_mc->GetErrorYhigh(pidx)));

                up_err_effratio.push_back(up_err);
                down_err_effratio.push_back(err_down);

                up_err_energy.push_back(0);
                down_err_energy.push_back(0);
            }

            std::shared_ptr<TGraphAsymmErrors> ratio = std::make_shared<TGraphAsymmErrors>(
                    static_cast<int>(energy.size()),
                    energy.data(),
                    mode_ratio.data(),
                    down_err_energy.data(),
                    up_err_energy.data(),
                    down_err_effratio.data(),
                    up_err_effratio.data()
                );
            
            ratio->SetName((std::string("ratio_") + std::string(eff_data->GetName())).c_str());
            ratio->GetXaxis()->SetTitle("Energy [GeV]");
            ratio->GetYaxis()->SetTitle("eff_{DATA}/eff_{MC}");

            return ratio;
        };

        // Divide data over mc efficiencies
        auto ratio_trigger_eff_het_over_let_bdt = get_efficiency_ratio(data_trigger_eff_het_over_let_bdt, mc_trigger_eff_het_over_let_bdt);
        auto ratio_trigger_eff_het_over_unb_bdt = get_efficiency_ratio(data_trigger_eff_het_over_unb_bdt, mc_trigger_eff_het_over_unb_bdt);
        auto ratio_maxrms_eff_bdt = get_efficiency_ratio(data_maxrms_eff_bdt, mc_maxrms_eff_bdt);
        auto ratio_nbarlayer13_eff_bdt = get_efficiency_ratio(data_nbarlayer13_eff_bdt, mc_nbarlayer13_eff_bdt);
        auto ratio_maxrms_and_nbarlayer13_eff_bdt = get_efficiency_ratio(data_maxrms_and_nbarlayer13_eff_bdt, mc_maxrms_and_nbarlayer13_eff_bdt);
        auto ratio_sumrms_low_energy_eff_bdt = get_efficiency_ratio(data_sumrms_low_energy_eff_bdt, mc_sumrms_low_energy_eff_bdt);
        auto ratio_track_selection_eff_bdt = get_efficiency_ratio(data_track_selection_eff_bdt, mc_track_selection_eff_bdt);
        auto ratio_track_selection_eff_within_stk_fvolume_bdt = get_efficiency_ratio(data_track_selection_eff_within_stk_fvolume_bdt, mc_track_selection_eff_within_stk_fvolume_bdt);
        auto ratio_stk_1_rm_eff_bdt = get_efficiency_ratio(data_stk_1_rm_eff_bdt, mc_stk_1_rm_eff_bdt); 
        auto ratio_psd_stk_match_eff_bdt = get_efficiency_ratio(data_psd_stk_match_eff_bdt, mc_psd_stk_match_eff_bdt);
        auto ratio_psd_charge_eff_bdt = get_efficiency_ratio(data_psd_charge_eff_bdt, mc_psd_charge_eff_bdt);
        auto ratio_stk_charge_eff_bdt = get_efficiency_ratio(data_stk_charge_eff_bdt, mc_stk_charge_eff_bdt);

        // Fit histos with a pol1
        TF1 f_ratio_trigger_eff_het_over_let_bdt("f_ratio_trigger_eff_het_over_let_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_trigger_eff_het_over_unb_bdt("f_ratio_trigger_eff_het_over_unb_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_maxrms_eff_bdt("f_ratio_maxrms_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_nbarlayer13_eff_bdt("f_ratio_nbarlayer13_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_maxrms_and_nbarlayer13_eff_bdt("f_ratio_maxrms_and_nbarlayer13_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_sumrms_low_energy_eff_bdt("f_ratio_sumrms_low_energy_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_track_selection_eff_bdt("f_ratio_track_selection_eff_bdt", "[0]", 10, 1e+4);
        TF1 f_ratio_track_selection_eff_within_stk_fvolume_bdt("f_ratio_track_selection_eff_within_stk_fvolume_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_stk_1_rm_eff_bdt("f_ratio_stk_1_rm_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_psd_stk_match_eff_bdt("f_ratio_psd_stk_match_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_psd_charge_eff_bdt("f_ratio_psd_charge_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);
        TF1 f_ratio_stk_charge_eff_bdt("f_ratio_stk_charge_eff_bdt", "[0]+[1]*log10(x)", 10, 1e+4);

        ratio_trigger_eff_het_over_let_bdt->Fit(&f_ratio_trigger_eff_het_over_let_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_trigger_eff_het_over_unb_bdt->Fit(&f_ratio_trigger_eff_het_over_unb_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_maxrms_eff_bdt->Fit(&f_ratio_maxrms_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_nbarlayer13_eff_bdt->Fit(&f_ratio_nbarlayer13_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_maxrms_and_nbarlayer13_eff_bdt->Fit(&f_ratio_maxrms_and_nbarlayer13_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_sumrms_low_energy_eff_bdt->Fit(&f_ratio_sumrms_low_energy_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_track_selection_eff_bdt->Fit(&f_ratio_track_selection_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_track_selection_eff_within_stk_fvolume_bdt->Fit(&f_ratio_track_selection_eff_within_stk_fvolume_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_stk_1_rm_eff_bdt->Fit(&f_ratio_stk_1_rm_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_psd_stk_match_eff_bdt->Fit(&f_ratio_psd_stk_match_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_psd_charge_eff_bdt->Fit(&f_ratio_psd_charge_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);
        ratio_stk_charge_eff_bdt->Fit(&f_ratio_stk_charge_eff_bdt, "QN", "", fit_low_energy_gev, fit_high_energy_gev);



        // Save results to the output file
        TFile output_eff_ratio_file(output_file, "RECREATE");
        if (output_eff_ratio_file.IsZombie()) {
            std::cerr << "\n\nError writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        ratio_trigger_eff_het_over_let_bdt->Write();
        ratio_trigger_eff_het_over_unb_bdt->Write();
        ratio_maxrms_eff_bdt->Write();
        ratio_nbarlayer13_eff_bdt->Write();
        ratio_maxrms_and_nbarlayer13_eff_bdt->Write();
        ratio_sumrms_low_energy_eff_bdt->Write();
        ratio_track_selection_eff_bdt->Write();
        ratio_track_selection_eff_within_stk_fvolume_bdt->Write();
        ratio_stk_1_rm_eff_bdt->Write();
        ratio_psd_stk_match_eff_bdt->Write();
        ratio_psd_charge_eff_bdt->Write();
        ratio_stk_charge_eff_bdt->Write();

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

        // Build canvases
        ratio_trigger_eff_het_over_let_bdt->SetMarkerStyle(8);
        ratio_trigger_eff_het_over_unb_bdt->SetMarkerStyle(8);
        ratio_maxrms_eff_bdt->SetMarkerStyle(8);
        ratio_nbarlayer13_eff_bdt->SetMarkerStyle(8);
        ratio_maxrms_and_nbarlayer13_eff_bdt->SetMarkerStyle(8);
        ratio_sumrms_low_energy_eff_bdt->SetMarkerStyle(8);
        ratio_track_selection_eff_bdt->SetMarkerStyle(8);
        ratio_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(8);
        ratio_stk_1_rm_eff_bdt->SetMarkerStyle(8);
        ratio_psd_stk_match_eff_bdt->SetMarkerStyle(8);
        ratio_psd_charge_eff_bdt->SetMarkerStyle(8);
        ratio_stk_charge_eff_bdt->SetMarkerStyle(8);

        ratio_trigger_eff_het_over_let_bdt->SetLineWidth(2);
        ratio_trigger_eff_het_over_unb_bdt->SetLineWidth(2);
        ratio_maxrms_eff_bdt->SetLineWidth(2);
        ratio_nbarlayer13_eff_bdt->SetLineWidth(2);
        ratio_maxrms_and_nbarlayer13_eff_bdt->SetLineWidth(2);
        ratio_sumrms_low_energy_eff_bdt->SetLineWidth(2);
        ratio_track_selection_eff_bdt->SetLineWidth(2);
        ratio_track_selection_eff_within_stk_fvolume_bdt->SetLineWidth(2);
        ratio_stk_1_rm_eff_bdt->SetLineWidth(2);
        ratio_psd_stk_match_eff_bdt->SetLineWidth(2);
        ratio_psd_charge_eff_bdt->SetLineWidth(2);
        ratio_stk_charge_eff_bdt->SetLineWidth(2);

        f_ratio_trigger_eff_het_over_let_bdt.SetLineWidth(3);
        f_ratio_trigger_eff_het_over_unb_bdt.SetLineWidth(3);
        f_ratio_maxrms_eff_bdt.SetLineWidth(3);
        f_ratio_nbarlayer13_eff_bdt.SetLineWidth(3);
        f_ratio_maxrms_and_nbarlayer13_eff_bdt.SetLineWidth(3);
        f_ratio_sumrms_low_energy_eff_bdt.SetLineWidth(3);
        f_ratio_track_selection_eff_bdt.SetLineWidth(3);
        f_ratio_track_selection_eff_within_stk_fvolume_bdt.SetLineWidth(3);
        f_ratio_stk_1_rm_eff_bdt.SetLineWidth(3);
        f_ratio_psd_stk_match_eff_bdt.SetLineWidth(3);
        f_ratio_psd_charge_eff_bdt.SetLineWidth(3);
        f_ratio_stk_charge_eff_bdt.SetLineWidth(3);

        f_ratio_trigger_eff_het_over_let_bdt.SetLineStyle(2);
        f_ratio_trigger_eff_het_over_unb_bdt.SetLineStyle(2);
        f_ratio_maxrms_eff_bdt.SetLineStyle(2);
        f_ratio_nbarlayer13_eff_bdt.SetLineStyle(2);
        f_ratio_maxrms_and_nbarlayer13_eff_bdt.SetLineStyle(2);
        f_ratio_sumrms_low_energy_eff_bdt.SetLineStyle(2);
        f_ratio_track_selection_eff_bdt.SetLineStyle(2);
        f_ratio_track_selection_eff_within_stk_fvolume_bdt.SetLineStyle(2);
        f_ratio_stk_1_rm_eff_bdt.SetLineStyle(2);
        f_ratio_psd_stk_match_eff_bdt.SetLineStyle(2);
        f_ratio_psd_charge_eff_bdt.SetLineStyle(2);
        f_ratio_stk_charge_eff_bdt.SetLineStyle(2);

        TCanvas c_het_let("c_het_let", "c_het_let", 500, 500);
        c_het_let.SetTicks();

        ratio_trigger_eff_het_over_let_bdt->Draw("AP");
        f_ratio_trigger_eff_het_over_let_bdt.Draw("same");

        gPad->Update(); 
        ratio_trigger_eff_het_over_let_bdt->SetMinimum(0.9);
        ratio_trigger_eff_het_over_let_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_het_let(0.0, 0.95, 0.3, 1, "HET over LET", "tlNDC");
        label_het_let.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_het_unb("c_het_unb", "c_het_unb", 500, 500);
        c_het_unb.SetTicks();

        ratio_trigger_eff_het_over_unb_bdt->Draw("AP");
        f_ratio_trigger_eff_het_over_unb_bdt.Draw("same");

        gPad->Update(); 
        ratio_trigger_eff_het_over_unb_bdt->SetMinimum(0.9);
        ratio_trigger_eff_het_over_unb_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_het_unb(0.0, 0.95, 0.3, 1, "HET over UNB", "tlNDC");
        label_het_unb.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_maxrms("c_maxrms", "c_maxrms", 500, 500);
        c_maxrms.SetTicks();

        ratio_maxrms_eff_bdt->Draw("AP");
        f_ratio_maxrms_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_maxrms_eff_bdt->SetMinimum(0.9);
        ratio_maxrms_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_maxrms(0.0, 0.95, 0.3, 1, "max RMS", "tlNDC");
        label_maxrms.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_nbarlayer13("c_nbarlayer13", "c_nbarlayer13", 500, 500);
        c_nbarlayer13.SetTicks();

        ratio_nbarlayer13_eff_bdt->Draw("AP");
        f_ratio_nbarlayer13_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_nbarlayer13_eff_bdt->SetMinimum(0.9);
        ratio_nbarlayer13_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_nbarlayer13(0.0, 0.95, 0.3, 1, "nbar - last layer", "tlNDC");
        label_nbarlayer13.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_trackselection("c_trackselection", "c_trackselection", 500, 500);
        c_trackselection.SetTicks();

        ratio_track_selection_eff_within_stk_fvolume_bdt->Draw("AP");
        f_ratio_track_selection_eff_within_stk_fvolume_bdt.Draw("same");

        gPad->Update(); 
        ratio_track_selection_eff_within_stk_fvolume_bdt->SetMinimum(0.9);
        ratio_track_selection_eff_within_stk_fvolume_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_trackselection(0.0, 0.95, 0.3, 1, "Track Selection", "tlNDC");
        label_trackselection.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_psdstkmatch("c_psdstkmatch", "c_psdstkmatch", 500, 500);
        c_psdstkmatch.SetTicks();

        ratio_psd_stk_match_eff_bdt->Draw("AP");
        f_ratio_psd_stk_match_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_psd_stk_match_eff_bdt->SetMinimum(0.9);
        ratio_psd_stk_match_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_psdstkmatch(0.0, 0.95, 0.3, 1, "PSD/STK match", "tlNDC");
        label_psdstkmatch.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_psdcharge("c_psdcharge", "c_hetc_psdcharge_let", 500, 500);
        c_psdcharge.SetTicks();

        ratio_psd_charge_eff_bdt->Draw("AP");
        f_ratio_psd_charge_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_psd_charge_eff_bdt->SetMinimum(0.9);
        ratio_psd_charge_eff_bdt->SetMaximum(1.1); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_psdcharge(0.0, 0.95, 0.3, 1, "PSD Charge", "tlNDC");
        label_psdcharge.Draw();
        gStyle->SetOptTitle(0);

        TCanvas c_stkcharge("c_stkcharge", "c_stkcharge", 500, 500);
        c_stkcharge.SetTicks();

        ratio_stk_charge_eff_bdt->Draw("AP");
        f_ratio_stk_charge_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_stk_charge_eff_bdt->SetMinimum(0.9);
        ratio_stk_charge_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        TPaveLabel label_stkcharge(0.0, 0.95, 0.3, 1, "STK Charge", "tlNDC");
        label_stkcharge.Draw();
        gStyle->SetOptTitle(0);

        // Write canvases

        c_het_let.Write();
        c_het_unb.Write();
        c_maxrms.Write();
        c_nbarlayer13.Write();
        c_trackselection.Write();
        c_psdstkmatch.Write();
        c_psdcharge.Write();
        c_stkcharge.Write();

        output_eff_ratio_file.Close();

        // Print summary PDF
        TCanvas print_canvas("print_canvas", "print_canvas");
        print_canvas.SetTicks();

        // Trigger Efficiency - HET over LET
        ratio_trigger_eff_het_over_let_bdt->Draw("AP");
        f_ratio_trigger_eff_het_over_let_bdt.Draw("same");

        gPad->Update(); 
        ratio_trigger_eff_het_over_let_bdt->SetMinimum(0.9);
        ratio_trigger_eff_het_over_let_bdt->SetMaximum(1.1); 
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
        ratio_trigger_eff_het_over_unb_bdt->Draw("AP");
        f_ratio_trigger_eff_het_over_unb_bdt.Draw("same");

        gPad->Update(); 
        ratio_trigger_eff_het_over_unb_bdt->SetMinimum(0.9);
        ratio_trigger_eff_het_over_unb_bdt->SetMaximum(1.1); 
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
        ratio_maxrms_eff_bdt->Draw("AP");
        f_ratio_maxrms_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_maxrms_eff_bdt->SetMinimum(0.9);
        ratio_maxrms_eff_bdt->SetMaximum(1.1); 
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
        ratio_nbarlayer13_eff_bdt->Draw("AP");
        f_ratio_nbarlayer13_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_nbarlayer13_eff_bdt->SetMinimum(0.9);
        ratio_nbarlayer13_eff_bdt->SetMaximum(1.1); 
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
        ratio_maxrms_and_nbarlayer13_eff_bdt->Draw("AP");
        f_ratio_maxrms_and_nbarlayer13_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_maxrms_and_nbarlayer13_eff_bdt->SetMinimum(0.9);
        ratio_maxrms_and_nbarlayer13_eff_bdt->SetMaximum(1.1); 
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
        ratio_sumrms_low_energy_eff_bdt->Draw("AP");
        f_ratio_sumrms_low_energy_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_sumrms_low_energy_eff_bdt->SetMinimum(0.9);
        ratio_sumrms_low_energy_eff_bdt->SetMaximum(1.1); 
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
        ratio_track_selection_eff_bdt->Draw("AP");
        f_ratio_track_selection_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_track_selection_eff_bdt->SetMinimum(0.9);
        ratio_track_selection_eff_bdt->SetMaximum(1.1); 
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
        ratio_track_selection_eff_within_stk_fvolume_bdt->Draw("AP");
        f_ratio_track_selection_eff_within_stk_fvolume_bdt.Draw("same");

        gPad->Update(); 
        ratio_track_selection_eff_within_stk_fvolume_bdt->SetMinimum(0.9);
        ratio_track_selection_eff_within_stk_fvolume_bdt->SetMaximum(1.1); 
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
        ratio_stk_1_rm_eff_bdt->Draw("AP");
        f_ratio_stk_1_rm_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_stk_1_rm_eff_bdt->SetMinimum(0.9);
        ratio_stk_1_rm_eff_bdt->SetMaximum(1.1); 
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
        ratio_psd_stk_match_eff_bdt->Draw("AP");
        f_ratio_psd_stk_match_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_psd_stk_match_eff_bdt->SetMinimum(0.9);
        ratio_psd_stk_match_eff_bdt->SetMaximum(1.1); 
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
        ratio_psd_charge_eff_bdt->Draw("AP");
        f_ratio_psd_charge_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_psd_charge_eff_bdt->SetMinimum(0.9);
        ratio_psd_charge_eff_bdt->SetMaximum(1.1); 
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
        ratio_stk_charge_eff_bdt->Draw("AP");
        f_ratio_stk_charge_eff_bdt.Draw("same");

        gPad->Update(); 
        ratio_stk_charge_eff_bdt->SetMinimum(0.9);
        ratio_stk_charge_eff_bdt->SetMaximum(1.1); 
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