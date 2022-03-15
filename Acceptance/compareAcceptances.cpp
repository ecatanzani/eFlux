#include <iostream>

#include "TPad.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TRootCanvas.h"
#include "TApplication.h"

void compareAcceptances(const char* acc_file, const char* acc_from_eff_file) {

    TFile* acc = TFile::Open(acc_file, "READ");
    if (acc->IsZombie()) {
        std::cerr << "\n\nError reading acceptance input file [" << acc_file << "]\n\n";
        exit(100);
    }

    TFile* eff = TFile::Open(acc_from_eff_file, "READ");
    if (eff->IsZombie()) {
        std::cerr << "\n\nError reading acceptance (efficiency) input file [" << acc_from_eff_file << "]\n\n";
        exit(100);
    }

    int argc;
    char** argv;
    TApplication window_application("app", &argc, argv);

    TCanvas c_fiducial("c_fiducial", "c_fiducial", 500, 500);
    acc->cd();
    auto acc_bgo_fiducial_het = static_cast<TGraphErrors*>(acc->Get("gr_acc_bgo_fiducial_het"));
    acc_bgo_fiducial_het->SetLineColor(kGreen+2);
    acc_bgo_fiducial_het->SetMarkerColor(kGreen+2);
    acc_bgo_fiducial_het->Draw();
    eff->cd();
    auto acc_bgo_fiducial_het_fe_xtrl_tight = static_cast<TGraph*>(eff->Get("gr_acceptance_bgo_fiducial_het_fe_xtrl_tight"));
    auto acc_bgo_fiducial_het_fe_xtrl_loose = static_cast<TGraph*>(eff->Get("gr_acceptance_bgo_fiducial_het_fe_xtrl_loose"));
    acc_bgo_fiducial_het_fe_xtrl_tight->SetLineColor(kBlue);
    acc_bgo_fiducial_het_fe_xtrl_tight->SetMarkerColor(kBlue);
    acc_bgo_fiducial_het_fe_xtrl_tight->SetLineStyle(3);
    acc_bgo_fiducial_het_fe_xtrl_tight->SetLineWidth(3);
    acc_bgo_fiducial_het_fe_xtrl_loose->SetLineColor(kMagenta);
    acc_bgo_fiducial_het_fe_xtrl_loose->SetMarkerColor(kMagenta);
    acc_bgo_fiducial_het_fe_xtrl_loose->SetLineStyle(3);
    acc_bgo_fiducial_het_fe_xtrl_loose->SetLineWidth(3);
    acc_bgo_fiducial_het_fe_xtrl_loose->Draw("same");
    gPad->SetLogx();

    c_fiducial.Modified();
    c_fiducial.Update();

    TCanvas c_nbarlayer13("c_nbarlayer13", "c_nbarlayer13", 500, 500);
    acc->cd();
    auto acc_nBarLayer13 = static_cast<TGraphErrors*>(acc->Get("gr_acc_nBarLayer13"));
    acc_nBarLayer13->SetLineColor(kGreen+2);
    acc_nBarLayer13->SetMarkerColor(kGreen+2);
    acc_nBarLayer13->Draw();
    eff->cd();
    auto acc_nbarlayer13_fe_xtrl_tight = static_cast<TGraph*>(eff->Get("gr_acceptance_nbarlayer13_fe_xtrl_tight"));
    auto acc_nbarlayer13_fe_xtrl_loose = static_cast<TGraph*>(eff->Get("gr_acceptance_nbarlayer13_fe_xtrl_loose"));
    acc_nbarlayer13_fe_xtrl_tight->SetLineColor(kBlue);
    acc_nbarlayer13_fe_xtrl_tight->SetMarkerColor(kBlue);
    acc_nbarlayer13_fe_xtrl_tight->SetLineStyle(3);
    acc_nbarlayer13_fe_xtrl_tight->SetLineWidth(3);
    acc_nbarlayer13_fe_xtrl_tight->Draw("same");
    acc_nbarlayer13_fe_xtrl_loose->SetLineColor(kMagenta);
    acc_nbarlayer13_fe_xtrl_loose->SetMarkerColor(kMagenta);
    acc_nbarlayer13_fe_xtrl_loose->SetLineStyle(3);
    acc_nbarlayer13_fe_xtrl_loose->SetLineWidth(3);
    acc_nbarlayer13_fe_xtrl_loose->Draw("same");
    gPad->SetLogx();

    c_nbarlayer13.Modified();
    c_nbarlayer13.Update();

    TCanvas c_maxrms("c_maxrms", "c_maxrms", 500, 500);
    acc->cd();
    auto acc_maxrms = static_cast<TGraphErrors*>(acc->Get("gr_acc_maxrms"));
    acc_maxrms->SetLineColor(kGreen+2);
    acc_maxrms->SetMarkerColor(kGreen+2);
    acc_maxrms->Draw();
    eff->cd();
    auto acc_maxrms_fe_xtrl_tight = static_cast<TGraph*>(eff->Get("gr_accepatance_maxrms_fe_xtrl_tight"));
    auto acc_maxrms_fe_xtrl_loose = static_cast<TGraph*>(eff->Get("gr_accepatance_maxrms_fe_xtrl_loose"));
    acc_maxrms_fe_xtrl_tight->SetLineColor(kBlue);
    acc_maxrms_fe_xtrl_tight->SetMarkerColor(kBlue);
    acc_maxrms_fe_xtrl_tight->SetLineStyle(3);
    acc_maxrms_fe_xtrl_tight->SetLineWidth(3);
    acc_maxrms_fe_xtrl_tight->Draw("same");
    acc_maxrms_fe_xtrl_loose->SetLineColor(kMagenta);
    acc_maxrms_fe_xtrl_loose->SetMarkerColor(kMagenta);
    acc_maxrms_fe_xtrl_loose->SetLineStyle(3);
    acc_maxrms_fe_xtrl_loose->SetLineWidth(3);
    acc_maxrms_fe_xtrl_loose->Draw("same");
    gPad->SetLogx();
    
    c_maxrms.Modified();
    c_maxrms.Update();

    TCanvas c_trackselection("c_trackselection", "c_trackselection", 500, 500);
    acc->cd();
    auto acc_trackselection = static_cast<TGraphErrors*>(acc->Get("gr_acc_trackselection"));
    acc_trackselection->SetLineColor(kGreen+2);
    acc_trackselection->SetMarkerColor(kGreen+2);
    acc_trackselection->Draw();
    eff->cd();
    auto acc_trackselection_fe_xtrl_tight = static_cast<TGraph*>(eff->Get("gr_acceptance_track_selection_fe_xtrl_tight"));
    auto acc_trackselection_fe_xtrl_loose = static_cast<TGraph*>(eff->Get("gr_acceptance_track_selection_fe_xtrl_loose"));
    acc_trackselection_fe_xtrl_tight->SetLineColor(kBlue);
    acc_trackselection_fe_xtrl_tight->SetMarkerColor(kBlue);
    acc_trackselection_fe_xtrl_tight->SetLineStyle(3);
    acc_trackselection_fe_xtrl_tight->SetLineWidth(3);
    acc_trackselection_fe_xtrl_tight->Draw("same");
    acc_trackselection_fe_xtrl_loose->SetLineColor(kMagenta);
    acc_trackselection_fe_xtrl_loose->SetMarkerColor(kMagenta);
    acc_trackselection_fe_xtrl_loose->SetLineStyle(3);
    acc_trackselection_fe_xtrl_loose->SetLineWidth(3);
    acc_trackselection_fe_xtrl_loose->Draw("same");
    gPad->SetLogx();
    
    c_trackselection.Modified();
    c_trackselection.Update();

    TCanvas c_psdstkmatch("c_psdstkmatch", "c_psdstkmatch", 500, 500);
    acc->cd();
    auto acc_psdstkmatch = static_cast<TGraphErrors*>(acc->Get("gr_acc_psdstkmatch"));
    acc_psdstkmatch->SetLineColor(kGreen+2);
    acc_psdstkmatch->SetMarkerColor(kGreen+2);
    acc_psdstkmatch->Draw();
    eff->cd();
    auto acc_psdstkmatch_fe_xtrl_tight = static_cast<TGraph*>(eff->Get("gr_acceptance_psd_stk_match_fe_xtrl_tight"));
    auto acc_psdstkmatch_fe_xtrl_loose = static_cast<TGraph*>(eff->Get("gr_acceptance_psd_stk_match_fe_xtrl_loose"));
    acc_psdstkmatch_fe_xtrl_tight->SetLineColor(kBlue);
    acc_psdstkmatch_fe_xtrl_tight->SetMarkerColor(kBlue);
    acc_psdstkmatch_fe_xtrl_tight->SetLineStyle(3);
    acc_psdstkmatch_fe_xtrl_tight->SetLineWidth(3);
    acc_psdstkmatch_fe_xtrl_tight->Draw("same");
    acc_psdstkmatch_fe_xtrl_loose->SetLineColor(kMagenta);
    acc_psdstkmatch_fe_xtrl_loose->SetMarkerColor(kMagenta);
    acc_psdstkmatch_fe_xtrl_loose->SetLineStyle(3);
    acc_psdstkmatch_fe_xtrl_loose->SetLineWidth(3);
    acc_psdstkmatch_fe_xtrl_loose->Draw("same");
    gPad->SetLogx();
    
    c_psdstkmatch.Modified();
    c_psdstkmatch.Update();

    TCanvas c_psdcharge("c_psdcharge", "c_psdcharge", 500, 500);
    acc->cd();
    auto acc_psdcharge = static_cast<TGraphErrors*>(acc->Get("gr_acc_all_cut"));
    acc_psdcharge->SetLineColor(kGreen+2);
    acc_psdcharge->SetMarkerColor(kGreen+2);
    acc_psdcharge->Draw();
    eff->cd();
    auto acc_psdcharge_fe_xtrl_tight = static_cast<TGraph*>(eff->Get("gr_acceptance_psd_charge_fe_xtrl_tight"));
    auto acc_psdcharge_fe_xtrl_loose = static_cast<TGraph*>(eff->Get("gr_acceptance_psd_charge_fe_xtrl_loose"));
    acc_psdcharge_fe_xtrl_tight->SetLineColor(kBlue);
    acc_psdcharge_fe_xtrl_tight->SetMarkerColor(kBlue);
    acc_psdcharge_fe_xtrl_tight->SetLineStyle(3);
    acc_psdcharge_fe_xtrl_tight->SetLineWidth(3);
    acc_psdcharge_fe_xtrl_tight->Draw("same");
    acc_psdcharge_fe_xtrl_loose->SetLineColor(kMagenta);
    acc_psdcharge_fe_xtrl_loose->SetMarkerColor(kMagenta);
    acc_psdcharge_fe_xtrl_loose->SetLineStyle(3);
    acc_psdcharge_fe_xtrl_loose->SetLineWidth(3);
    acc_psdcharge_fe_xtrl_loose->Draw("same");
    gPad->SetLogx();
    
    c_psdcharge.Modified();
    c_psdcharge.Update();

    TRootCanvas* rc_bgo_fiducial = static_cast<TRootCanvas*>(c_fiducial.GetCanvasImp());
    rc_bgo_fiducial->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TRootCanvas* rc_nbarlayer13 = static_cast<TRootCanvas*>(c_nbarlayer13.GetCanvasImp());
    rc_nbarlayer13->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TRootCanvas* rc_maxrms = static_cast<TRootCanvas*>(c_maxrms.GetCanvasImp());
    rc_maxrms->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TRootCanvas* rc_trackselection = static_cast<TRootCanvas*>(c_trackselection.GetCanvasImp());
    rc_trackselection->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TRootCanvas* rc_psdstkmatch = static_cast<TRootCanvas*>(c_psdstkmatch.GetCanvasImp());
    rc_psdstkmatch->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TRootCanvas* rc_psdcharge = static_cast<TRootCanvas*>(c_psdcharge.GetCanvasImp());
    rc_psdcharge->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    window_application.Run();
}