#include <iostream>

#include "TPad.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TEfficiency.h"
#include "TApplication.h"

void produceEfficiencyPlots(
    const char* data_file,
    const char* mc_file,
    const bool verbose) {

        TFile* datafile = TFile::Open(data_file, "READ");
        if (datafile->IsZombie()) {
            std::cerr << "\n\nError opening input DATA ROOT file [" << data_file << "]\n\n";
            exit(100);
        }

        auto data_trigger_eff_het_over_let_xtrl_tight       = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_let_xtrl_tight"));
        auto data_trigger_eff_het_over_unb_xtrl_tight       = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_tight"));
        auto data_maxrms_eff_xtrl_tight                     = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_eff_xtrl_tight"));
        auto data_nbarlayer13_eff_xtrl_tight                = static_cast<TEfficiency*>(datafile->Get("efficiencies/nbarlayer13_eff_xtrl_tight"));
        auto data_maxrms_and_nbarlayer13_eff_xtrl_tight     = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_tight"));
        auto data_track_selection_eff_xtrl_tight            = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_xtrl_tight"));
        auto data_psd_stk_match_eff_xtrl_tight              = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_stk_match_eff_xtrl_tight"));
        auto data_psd_charge_eff_xtrl_tight                 = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_charge_eff_xtrl_tight"));
        auto data_trigger_eff_het_over_let_xtrl_loose       = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_let_xtrl_loose"));
        auto data_trigger_eff_het_over_unb_xtrl_loose       = static_cast<TEfficiency*>(datafile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_loose"));
        auto data_maxrms_eff_xtrl_loose                     = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_eff_xtrl_loose"));
        auto data_nbarlayer13_eff_xtrl_loose                = static_cast<TEfficiency*>(datafile->Get("efficiencies/nbarlayer13_eff_xtrl_loose"));
        auto data_maxrms_and_nbarlayer13_eff_xtrl_loose     = static_cast<TEfficiency*>(datafile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_loose"));
        auto data_track_selection_eff_xtrl_loose            = static_cast<TEfficiency*>(datafile->Get("efficiencies/track_selection_eff_xtrl_loose"));
        auto data_psd_stk_match_eff_xtrl_loose              = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_stk_match_eff_xtrl_loose"));
        auto data_psd_charge_eff_xtrl_loose                 = static_cast<TEfficiency*>(datafile->Get("efficiencies/psd_charge_eff_xtrl_loose"));

        data_trigger_eff_het_over_let_xtrl_tight            ->SetLineColor(kRed);
        data_trigger_eff_het_over_unb_xtrl_tight            ->SetLineColor(kRed);
        data_maxrms_eff_xtrl_tight                          ->SetLineColor(kRed);
        data_nbarlayer13_eff_xtrl_tight                     ->SetLineColor(kRed);
        data_maxrms_and_nbarlayer13_eff_xtrl_tight          ->SetLineColor(kRed);
        data_track_selection_eff_xtrl_tight                 ->SetLineColor(kRed);
        data_psd_stk_match_eff_xtrl_tight                   ->SetLineColor(kRed);
        data_psd_charge_eff_xtrl_tight                      ->SetLineColor(kRed);
        data_trigger_eff_het_over_let_xtrl_loose            ->SetLineColor(kRed);
        data_trigger_eff_het_over_unb_xtrl_loose            ->SetLineColor(kRed);
        data_maxrms_eff_xtrl_loose                          ->SetLineColor(kRed);
        data_nbarlayer13_eff_xtrl_loose                     ->SetLineColor(kRed);
        data_maxrms_and_nbarlayer13_eff_xtrl_loose          ->SetLineColor(kRed);
        data_track_selection_eff_xtrl_loose                 ->SetLineColor(kRed);
        data_psd_stk_match_eff_xtrl_loose                   ->SetLineColor(kRed);
        data_psd_charge_eff_xtrl_loose                      ->SetLineColor(kRed);

        data_trigger_eff_het_over_let_xtrl_tight            ->SetDirectory(0);
        data_trigger_eff_het_over_unb_xtrl_tight            ->SetDirectory(0);
        data_maxrms_eff_xtrl_tight                          ->SetDirectory(0);
        data_nbarlayer13_eff_xtrl_tight                     ->SetDirectory(0);
        data_maxrms_and_nbarlayer13_eff_xtrl_tight          ->SetDirectory(0);
        data_track_selection_eff_xtrl_tight                 ->SetDirectory(0);
        data_psd_stk_match_eff_xtrl_tight                   ->SetDirectory(0);
        data_psd_charge_eff_xtrl_tight                      ->SetDirectory(0);
        data_trigger_eff_het_over_let_xtrl_loose            ->SetDirectory(0);
        data_trigger_eff_het_over_unb_xtrl_loose            ->SetDirectory(0);
        data_maxrms_eff_xtrl_loose                          ->SetDirectory(0);
        data_nbarlayer13_eff_xtrl_loose                     ->SetDirectory(0);
        data_maxrms_and_nbarlayer13_eff_xtrl_loose          ->SetDirectory(0);
        data_track_selection_eff_xtrl_loose                 ->SetDirectory(0);
        data_psd_stk_match_eff_xtrl_loose                   ->SetDirectory(0);
        data_psd_charge_eff_xtrl_loose                      ->SetDirectory(0);

        datafile->Close();

        TFile* mcfile = TFile::Open(mc_file, "READ");
        if (mcfile->IsZombie()) {
            std::cerr << "\n\nError opening input MC ROOT file [" << mc_file << "]\n\n";
            exit(100);
        }

        auto mc_trigger_eff_het_over_let_xtrl_tight       = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_let_xtrl_tight"));
        auto mc_trigger_eff_het_over_unb_xtrl_tight       = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_tight"));
        auto mc_maxrms_eff_xtrl_tight                     = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_eff_xtrl_tight"));
        auto mc_nbarlayer13_eff_xtrl_tight                = static_cast<TEfficiency*>(mcfile->Get("efficiencies/nbarlayer13_eff_xtrl_tight"));
        auto mc_maxrms_and_nbarlayer13_eff_xtrl_tight     = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_tight"));
        auto mc_track_selection_eff_xtrl_tight            = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_xtrl_tight"));
        auto mc_psd_stk_match_eff_xtrl_tight              = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_stk_match_eff_xtrl_tight"));
        auto mc_psd_charge_eff_xtrl_tight                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_charge_eff_xtrl_tight"));
        auto mc_trigger_eff_het_over_let_xtrl_loose       = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_let_xtrl_loose"));
        auto mc_trigger_eff_het_over_unb_xtrl_loose       = static_cast<TEfficiency*>(mcfile->Get("efficiencies/trigger_eff_het_over_unb_xtrl_loose"));
        auto mc_maxrms_eff_xtrl_loose                     = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_eff_xtrl_loose"));
        auto mc_nbarlayer13_eff_xtrl_loose                = static_cast<TEfficiency*>(mcfile->Get("efficiencies/nbarlayer13_eff_xtrl_loose"));
        auto mc_maxrms_and_nbarlayer13_eff_xtrl_loose     = static_cast<TEfficiency*>(mcfile->Get("efficiencies/maxrms_and_nbarlayer13_eff_xtrl_loose"));
        auto mc_track_selection_eff_xtrl_loose            = static_cast<TEfficiency*>(mcfile->Get("efficiencies/track_selection_eff_xtrl_loose"));
        auto mc_psd_stk_match_eff_xtrl_loose              = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_stk_match_eff_xtrl_loose"));
        auto mc_psd_charge_eff_xtrl_loose                 = static_cast<TEfficiency*>(mcfile->Get("efficiencies/psd_charge_eff_xtrl_loose"));

        mc_trigger_eff_het_over_let_xtrl_tight            ->SetLineColor(kBlue);
        mc_trigger_eff_het_over_unb_xtrl_tight            ->SetLineColor(kBlue);
        mc_maxrms_eff_xtrl_tight                          ->SetLineColor(kBlue);
        mc_nbarlayer13_eff_xtrl_tight                     ->SetLineColor(kBlue);
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight          ->SetLineColor(kBlue);
        mc_track_selection_eff_xtrl_tight                 ->SetLineColor(kBlue);
        mc_psd_stk_match_eff_xtrl_tight                   ->SetLineColor(kBlue);
        mc_psd_charge_eff_xtrl_tight                      ->SetLineColor(kBlue);
        mc_trigger_eff_het_over_let_xtrl_loose            ->SetLineColor(kBlue);
        mc_trigger_eff_het_over_unb_xtrl_loose            ->SetLineColor(kBlue);
        mc_maxrms_eff_xtrl_loose                          ->SetLineColor(kBlue);
        mc_nbarlayer13_eff_xtrl_loose                     ->SetLineColor(kBlue);
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose          ->SetLineColor(kBlue);
        mc_track_selection_eff_xtrl_loose                 ->SetLineColor(kBlue);
        mc_psd_stk_match_eff_xtrl_loose                   ->SetLineColor(kBlue);
        mc_psd_charge_eff_xtrl_loose                      ->SetLineColor(kBlue);

        mc_trigger_eff_het_over_let_xtrl_tight            ->SetDirectory(0);
        mc_trigger_eff_het_over_unb_xtrl_tight            ->SetDirectory(0);
        mc_maxrms_eff_xtrl_tight                          ->SetDirectory(0);
        mc_nbarlayer13_eff_xtrl_tight                     ->SetDirectory(0);
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight          ->SetDirectory(0);
        mc_track_selection_eff_xtrl_tight                 ->SetDirectory(0);
        mc_psd_stk_match_eff_xtrl_tight                   ->SetDirectory(0);
        mc_psd_charge_eff_xtrl_tight                      ->SetDirectory(0);
        mc_trigger_eff_het_over_let_xtrl_loose            ->SetDirectory(0);
        mc_trigger_eff_het_over_unb_xtrl_loose            ->SetDirectory(0);
        mc_maxrms_eff_xtrl_loose                          ->SetDirectory(0);
        mc_nbarlayer13_eff_xtrl_loose                     ->SetDirectory(0);
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose          ->SetDirectory(0);
        mc_track_selection_eff_xtrl_loose                 ->SetDirectory(0);
        mc_psd_stk_match_eff_xtrl_loose                   ->SetDirectory(0);
        mc_psd_charge_eff_xtrl_loose                      ->SetDirectory(0);

        datafile->Close();

        int argc;
        char** argv;
        TApplication window_application("app", &argc, argv);
        
        TCanvas eff_canvas("eff_canvas", "Efficiency Plots");
        eff_canvas.Divide(4, 4);
        
        eff_canvas.cd(1);
        data_trigger_eff_het_over_let_xtrl_tight->Draw();
        mc_trigger_eff_het_over_let_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(2);
        data_trigger_eff_het_over_unb_xtrl_tight->Draw();
        mc_trigger_eff_het_over_unb_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(3);
        data_maxrms_eff_xtrl_tight->Draw();
        mc_maxrms_eff_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(4);
        data_nbarlayer13_eff_xtrl_tight->Draw();
        mc_nbarlayer13_eff_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(5);
        data_maxrms_and_nbarlayer13_eff_xtrl_tight->Draw();
        mc_maxrms_and_nbarlayer13_eff_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(6);
        data_track_selection_eff_xtrl_tight->Draw();
        mc_track_selection_eff_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(7);
        data_psd_stk_match_eff_xtrl_tight->Draw();
        mc_psd_stk_match_eff_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(8);
        data_psd_charge_eff_xtrl_tight->Draw();
        mc_psd_charge_eff_xtrl_tight->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(9);
        data_trigger_eff_het_over_let_xtrl_loose->Draw();
        mc_trigger_eff_het_over_let_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(10);
        data_trigger_eff_het_over_unb_xtrl_loose->Draw();
        mc_trigger_eff_het_over_unb_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(11);
        data_maxrms_eff_xtrl_loose->Draw();
        mc_maxrms_eff_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(12);
        data_nbarlayer13_eff_xtrl_loose->Draw();
        mc_nbarlayer13_eff_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(13);
        data_maxrms_and_nbarlayer13_eff_xtrl_loose->Draw();
        mc_maxrms_and_nbarlayer13_eff_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(14);
        data_track_selection_eff_xtrl_loose->Draw();
        mc_track_selection_eff_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(15);
        data_psd_stk_match_eff_xtrl_loose->Draw();
        mc_psd_stk_match_eff_xtrl_loose->Draw("same");
        gPad->SetLogx();
        eff_canvas.cd(16);
        data_psd_charge_eff_xtrl_loose->Draw();
        mc_psd_charge_eff_xtrl_loose->Draw("same");
        gPad->SetLogx();

        eff_canvas.Modified(); 
        eff_canvas.Update();

        TRootCanvas* rc = static_cast<TRootCanvas*>(eff_canvas.GetCanvasImp());
        rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        window_application.Run();
    }