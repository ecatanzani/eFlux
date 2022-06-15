#include <iostream>
#include <string>

#include "TPad.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TPaveLabel.h"
#include "TEfficiency.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"

void produceEfficiencyPlots(
    const char* data_file,
    const char* mc_file,
    const char* output_file = "efficiency_plots_canvases.root",
    const char* data_label="DATA",
    const char* mc_label="electron MC") {

        TFile* datafile = TFile::Open(data_file, "READ");
        if (datafile->IsZombie()) {
            std::cerr << "\n\nError opening input DATA ROOT file [" << data_file << "]\n\n";
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

        data_trigger_eff_het_over_let_bdt                       ->SetLineColor(kRed);
        data_trigger_eff_het_over_unb_bdt                       ->SetLineColor(kRed);
        data_maxrms_eff_bdt                                     ->SetLineColor(kRed);
        data_nbarlayer13_eff_bdt                                ->SetLineColor(kRed);
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetLineColor(kRed);
        data_sumrms_low_energy_eff_bdt                          ->SetLineColor(kRed);
        data_track_selection_eff_bdt                            ->SetLineColor(kRed);
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetLineColor(kRed);
        data_stk_1_rm_eff_bdt                                   ->SetLineColor(kRed);
        data_psd_stk_match_eff_bdt                              ->SetLineColor(kRed);
        data_psd_charge_eff_bdt                                 ->SetLineColor(kRed);
        data_stk_charge_eff_bdt                                 ->SetLineColor(kRed);

        data_trigger_eff_het_over_let_bdt                       ->SetLineWidth(2);
        data_trigger_eff_het_over_unb_bdt                       ->SetLineWidth(2);
        data_maxrms_eff_bdt                                     ->SetLineWidth(2);
        data_nbarlayer13_eff_bdt                                ->SetLineWidth(2);
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetLineWidth(2);
        data_sumrms_low_energy_eff_bdt                          ->SetLineWidth(2);
        data_track_selection_eff_bdt                            ->SetLineWidth(2);
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetLineWidth(2);
        data_stk_1_rm_eff_bdt                                   ->SetLineWidth(2);
        data_psd_stk_match_eff_bdt                              ->SetLineWidth(2);
        data_psd_charge_eff_bdt                                 ->SetLineWidth(2);
        data_stk_charge_eff_bdt                                 ->SetLineWidth(2);

        data_trigger_eff_het_over_let_bdt                       ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_trigger_eff_het_over_unb_bdt                       ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_maxrms_eff_bdt                                     ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_nbarlayer13_eff_bdt                                ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_maxrms_and_nbarlayer13_eff_bdt                     ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_sumrms_low_energy_eff_bdt                          ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_track_selection_eff_bdt                            ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_track_selection_eff_within_stk_fvolume_bdt         ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_stk_1_rm_eff_bdt                                   ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_psd_stk_match_eff_bdt                              ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_psd_charge_eff_bdt                                 ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());
        data_stk_charge_eff_bdt                                 ->SetTitle((std::string(data_label) + std::string(" (BDT)")).c_str());

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

        TFile* mcfile = TFile::Open(mc_file, "READ");
        if (mcfile->IsZombie()) {
            std::cerr << "\n\nError opening input MC ROOT file [" << mc_file << "]\n\n";
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

        mc_trigger_eff_het_over_let_bdt                         ->SetLineColor(kBlue);
        mc_trigger_eff_het_over_unb_bdt                         ->SetLineColor(kBlue);
        mc_maxrms_eff_bdt                                       ->SetLineColor(kBlue);
        mc_nbarlayer13_eff_bdt                                  ->SetLineColor(kBlue);
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetLineColor(kBlue);
        mc_sumrms_low_energy_eff_bdt                            ->SetLineColor(kBlue);
        mc_track_selection_eff_bdt                              ->SetLineColor(kBlue);
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetLineColor(kBlue);
        mc_stk_1_rm_eff_bdt                                     ->SetLineColor(kBlue);
        mc_psd_stk_match_eff_bdt                                ->SetLineColor(kBlue);
        mc_psd_charge_eff_bdt                                   ->SetLineColor(kBlue);
        mc_stk_charge_eff_bdt                                   ->SetLineColor(kBlue);

        mc_trigger_eff_het_over_let_bdt                         ->SetLineWidth(2);
        mc_trigger_eff_het_over_unb_bdt                         ->SetLineWidth(2);
        mc_maxrms_eff_bdt                                       ->SetLineWidth(2);
        mc_nbarlayer13_eff_bdt                                  ->SetLineWidth(2);
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetLineWidth(2);
        mc_sumrms_low_energy_eff_bdt                            ->SetLineWidth(2);
        mc_track_selection_eff_bdt                              ->SetLineWidth(2);
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetLineWidth(2);
        mc_stk_1_rm_eff_bdt                                     ->SetLineWidth(2);
        mc_psd_stk_match_eff_bdt                                ->SetLineWidth(2);
        mc_psd_charge_eff_bdt                                   ->SetLineWidth(2);
        mc_stk_charge_eff_bdt                                   ->SetLineWidth(2);

        mc_trigger_eff_het_over_let_bdt                         ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_trigger_eff_het_over_unb_bdt                         ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_maxrms_eff_bdt                                       ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_nbarlayer13_eff_bdt                                  ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_maxrms_and_nbarlayer13_eff_bdt                       ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_sumrms_low_energy_eff_bdt                            ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_track_selection_eff_bdt                              ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_track_selection_eff_within_stk_fvolume_bdt           ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_stk_1_rm_eff_bdt                                     ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_psd_stk_match_eff_bdt                                ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_psd_charge_eff_bdt                                   ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());
        mc_stk_charge_eff_bdt                                   ->SetTitle((std::string(mc_label) + std::string(" (BDT)")).c_str());

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

        auto remove_x_error_from_TEfficiency = [](TEfficiency* eff) -> std::shared_ptr<TGraphAsymmErrors> {
            auto gr = std::make_shared<TGraphAsymmErrors>(*static_cast<TGraphAsymmErrors*>(eff->CreateGraph()));
            for (int i = 0; i < gr->GetN(); ++i) {
                gr->GetEXlow()[i] = 0;  
                gr->GetEXhigh()[i] = 0;  
            }
            return gr;
        };

        // Print trigger HET over LET efficiency
        TCanvas c_eff_het_let("c_eff_het_let", "c_eff_het_let", 500, 500);
        
        auto gr_data_trigger_eff_het_over_let_bdt = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_let_bdt);
        auto gr_mc_trigger_eff_het_over_let_bdt = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_let_bdt);
        
        gr_data_trigger_eff_het_over_let_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_let_bdt->SetMarkerStyle(kFullDotMedium);

        gr_data_trigger_eff_het_over_let_bdt->SetMarkerColor(kRed);
        gr_mc_trigger_eff_het_over_let_bdt->SetMarkerColor(kBlue);

        gr_data_trigger_eff_het_over_let_bdt->SetMarkerStyle(20);
        gr_mc_trigger_eff_het_over_let_bdt->SetMarkerStyle(20);
        
        gr_data_trigger_eff_het_over_let_bdt->Draw("AP");
        gr_mc_trigger_eff_het_over_let_bdt->Draw("P");

        gPad->Update(); 
        gr_data_trigger_eff_het_over_let_bdt->SetMinimum(0.9);
        gr_data_trigger_eff_het_over_let_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        auto legend = c_eff_het_let.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_het_let(0.0, 0.95, 0.3, 1, "HET over LET", "tlNDC");
        label_c_eff_het_let.Draw();

        // Print trigger HET over UNB efficiency
        TCanvas c_eff_het_unb("c_eff_het_unb", "c_eff_het_unb", 500, 500);
        
        auto gr_data_trigger_eff_het_over_unb_bdt = remove_x_error_from_TEfficiency(data_trigger_eff_het_over_unb_bdt);
        auto gr_mc_trigger_eff_het_over_unb_bdt = remove_x_error_from_TEfficiency(mc_trigger_eff_het_over_unb_bdt);
        
        gr_data_trigger_eff_het_over_unb_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_trigger_eff_het_over_unb_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_trigger_eff_het_over_unb_bdt->SetMarkerColor(kRed);
        gr_mc_trigger_eff_het_over_unb_bdt->SetMarkerColor(kBlue);

        gr_data_trigger_eff_het_over_unb_bdt->SetMarkerStyle(20);
        gr_mc_trigger_eff_het_over_unb_bdt->SetMarkerStyle(20);
        
        gr_data_trigger_eff_het_over_unb_bdt->Draw("AP");
        gr_mc_trigger_eff_het_over_unb_bdt->Draw("P");

        gPad->Update(); 
        gr_data_trigger_eff_het_over_unb_bdt->SetMinimum(0.9);
        gr_data_trigger_eff_het_over_unb_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_het_unb.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_het_unb(0.0, 0.95, 0.3, 1, "HET over UNB", "tlNDC");
        label_c_eff_het_unb.Draw();

        // Print max RMS efficiency
        TCanvas c_eff_max_rms("c_eff_max_rms", "c_eff_max_rms", 500, 500);
        
        auto gr_data_maxrms_eff_bdt = remove_x_error_from_TEfficiency(data_maxrms_eff_bdt);
        auto gr_mc_maxrms_eff_bdt = remove_x_error_from_TEfficiency(mc_maxrms_eff_bdt);
        
        gr_data_maxrms_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_maxrms_eff_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_maxrms_eff_bdt->SetMarkerColor(kRed);
        gr_mc_maxrms_eff_bdt->SetMarkerColor(kBlue);

        gr_data_maxrms_eff_bdt->SetMarkerStyle(20);
        gr_mc_maxrms_eff_bdt->SetMarkerStyle(20);
        
        gr_data_maxrms_eff_bdt->Draw("AP");
        gr_mc_maxrms_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_maxrms_eff_bdt->SetMinimum(0.9);
        gr_data_maxrms_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_max_rms.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_max_rms(0.0, 0.95, 0.3, 1, "max RMS", "tlNDC");
        label_c_eff_max_rms.Draw();

        // Print nbarlayer13 efficiency
        TCanvas c_eff_nbarlayer13("c_eff_nbarlayer13", "c_eff_nbarlayer13", 500, 500);
        
        auto gr_data_nbarlayer13_eff_bdt = remove_x_error_from_TEfficiency(data_nbarlayer13_eff_bdt);
        auto gr_mc_nbarlayer13_eff_bdt = remove_x_error_from_TEfficiency(mc_nbarlayer13_eff_bdt);
        
        gr_data_nbarlayer13_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_nbarlayer13_eff_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_nbarlayer13_eff_bdt->SetMarkerColor(kRed);
        gr_mc_nbarlayer13_eff_bdt->SetMarkerColor(kBlue);

        gr_data_nbarlayer13_eff_bdt->SetMarkerStyle(20);
        gr_mc_nbarlayer13_eff_bdt->SetMarkerStyle(20);
        
        gr_data_nbarlayer13_eff_bdt->Draw("AP");
        gr_mc_nbarlayer13_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_nbarlayer13_eff_bdt->SetMinimum(0.9);
        gr_data_nbarlayer13_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_nbarlayer13.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_nbarlayer13(0.0, 0.95, 0.3, 1, "max RMS", "tlNDC");
        label_c_eff_nbarlayer13.Draw();

        // Print trackselection efficiency within STK fiducial volume
        TCanvas c_eff_trackselection_stk("c_eff_trackselection_stk", "c_eff_trackselection_stk", 500, 500);
        
        auto gr_data_track_selection_eff_within_stk_fvolume_bdt = remove_x_error_from_TEfficiency(data_track_selection_eff_within_stk_fvolume_bdt);
        auto gr_mc_track_selection_eff_within_stk_fvolume_bdt = remove_x_error_from_TEfficiency(mc_track_selection_eff_within_stk_fvolume_bdt);
        
        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMarkerColor(kRed);
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->SetMarkerColor(kBlue);

        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(20);
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->SetMarkerStyle(20);
        
        gr_data_track_selection_eff_within_stk_fvolume_bdt->Draw("AP");
        gr_mc_track_selection_eff_within_stk_fvolume_bdt->Draw("P");

        gPad->Update(); 
        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMinimum(0.9);
        gr_data_track_selection_eff_within_stk_fvolume_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_trackselection_stk.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_trackselection_stk(0.0, 0.95, 0.3, 1, "Track Selection", "tlNDC");
        label_c_eff_trackselection_stk.Draw();

        // Print psd-stk match efficiency
        TCanvas c_eff_psd_stk("c_eff_psd_stk", "c_eff_psd_stk", 500, 500);
        
        auto gr_data_psd_stk_match_eff_bdt = remove_x_error_from_TEfficiency(data_psd_stk_match_eff_bdt);
        auto gr_mc_psd_stk_match_eff_bdt = remove_x_error_from_TEfficiency(mc_psd_stk_match_eff_bdt);
        
        gr_data_psd_stk_match_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_stk_match_eff_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_psd_stk_match_eff_bdt->SetMarkerColor(kRed);
        gr_mc_psd_stk_match_eff_bdt->SetMarkerColor(kBlue);

        gr_data_psd_stk_match_eff_bdt->SetMarkerStyle(20);
        gr_mc_psd_stk_match_eff_bdt->SetMarkerStyle(20);
        
        gr_data_psd_stk_match_eff_bdt->Draw("AP");
        gr_mc_psd_stk_match_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_psd_stk_match_eff_bdt->SetMinimum(0.9);
        gr_data_psd_stk_match_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_psd_stk.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }
        
        TPaveLabel label_c_eff_psd_stk(0.0, 0.95, 0.3, 1, "PSD/STK match", "tlNDC");
        label_c_eff_psd_stk.Draw();

        // Print psd charge efficiency
        TCanvas c_eff_psd_charge("c_eff_psd_charge", "c_eff_psd_charge", 500, 500);
        
        auto gr_data_psd_charge_eff_bdt = remove_x_error_from_TEfficiency(data_psd_charge_eff_bdt);
        auto gr_mc_psd_charge_eff_bdt = remove_x_error_from_TEfficiency(mc_psd_charge_eff_bdt);
        
        gr_data_psd_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_psd_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_psd_charge_eff_bdt->SetMarkerColor(kRed);
        gr_mc_psd_charge_eff_bdt->SetMarkerColor(kBlue);

        gr_data_psd_charge_eff_bdt->SetMarkerStyle(20);
        gr_mc_psd_charge_eff_bdt->SetMarkerStyle(20);
        
        gr_data_psd_charge_eff_bdt->Draw("AP");
        gr_mc_psd_charge_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_psd_charge_eff_bdt->SetMinimum(0.9);
        gr_data_psd_charge_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_psd_charge.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_psd_charge(0.0, 0.95, 0.3, 1, "PSD Charge", "tlNDC");
        label_c_eff_psd_charge.Draw();
        
        // Print stk charge efficiency 
        TCanvas c_eff_stk_charge("c_eff_stk_charge", "c_eff_stk_charge", 500, 500);
        
        auto gr_data_stk_charge_eff_bdt = remove_x_error_from_TEfficiency(data_stk_charge_eff_bdt);
        auto gr_mc_stk_charge_eff_bdt = remove_x_error_from_TEfficiency(mc_stk_charge_eff_bdt);
        
        gr_data_stk_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);
        gr_mc_stk_charge_eff_bdt->SetMarkerStyle(kFullDotMedium);
        
        gr_data_stk_charge_eff_bdt->SetMarkerColor(kRed);
        gr_mc_stk_charge_eff_bdt->SetMarkerColor(kBlue);

        gr_data_stk_charge_eff_bdt->SetMarkerStyle(20);
        gr_mc_psd_charge_eff_bdt->SetMarkerStyle(20);

        gr_data_stk_charge_eff_bdt->Draw("AP");
        gr_mc_stk_charge_eff_bdt->Draw("P");

        gPad->Update(); 
        gr_data_stk_charge_eff_bdt->SetMinimum(0.95);
        gr_data_stk_charge_eff_bdt->SetMaximum(1.01); 
        gPad->Update();

        gPad->SetLogx();
        gPad->SetGrid(1,1);

        legend = c_eff_stk_charge.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("p");
        }

        TPaveLabel label_c_eff_stk_charge(0.0, 0.95, 0.3, 1, "STK Charge", "tlNDC");
        label_c_eff_stk_charge.Draw();

        TFile output(output_file, "RECREATE");
        if (output.IsZombie())
        {
            std::cerr << "Error writing output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        c_eff_het_let.Write();
        c_eff_het_unb.Write();
        c_eff_max_rms.Write();
        c_eff_nbarlayer13.Write();
        c_eff_trackselection_stk.Write();
        c_eff_psd_stk.Write();
        c_eff_psd_charge.Write();
        c_eff_stk_charge.Write();

        output.Close();

    }