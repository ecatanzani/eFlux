#include <iostream>

#include "TPad.h"
#include "TPDF.h"
#include "TFile.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveLabel.h"

void buildROCPlots(
    const char* roc_10_100,
    const char* roc_100_250,
    const char* roc_250_500,
    const char* roc_500_1000,
    const char* roc_1000_3000,
    const char* roc_3000)
    {
        TCanvas print_canvas("print_canvas", "print_canvas");
        print_canvas.SetTicks();
        
        print_canvas.Divide(2, 1);

        TFile input_file_10_100(roc_10_100, "READ");
        if (input_file_10_100.IsZombie())
        {
            std::cerr << "\n\nError opening input file [" << roc_10_100 << "]\n\n";
            exit(100);
        }

        auto gr_roc_10_100_xtrl = static_cast<TGraph*>(input_file_10_100.Get("gr_xtrl_roc"));
        auto gr_roc_10_100_bdt = static_cast<TGraph*>(input_file_10_100.Get("gr_bdt_roc"));
        auto gr_rpower_10_100_xtrl = static_cast<TGraph*>(input_file_10_100.Get("gr_xtrl_proton_rejection_power"));
        auto gr_rpower_10_100_bdt = static_cast<TGraph*>(input_file_10_100.Get("gr_bdt_proton_rejection_power"));

        gr_roc_10_100_xtrl->SetLineWidth(2);
        gr_roc_10_100_xtrl->SetLineColor(kBlue);
        gr_roc_10_100_bdt->SetLineWidth(2);
        gr_roc_10_100_bdt->SetLineColor(kRed);
        gr_rpower_10_100_xtrl->SetLineWidth(2);
        gr_rpower_10_100_xtrl->SetLineColor(kBlue);
        gr_rpower_10_100_bdt->SetLineWidth(2);
        gr_rpower_10_100_bdt->SetLineColor(kRed);

        print_canvas.cd(1);
        gr_roc_10_100_xtrl->GetXaxis()->SetRangeUser(0, 1.05);
        gr_roc_10_100_xtrl->GetYaxis()->SetRangeUser(0.98, 1.001);
        gr_roc_10_100_xtrl->Draw();
        gr_roc_10_100_bdt->Draw("same");
        gPad->SetGrid(1,1);
        print_canvas.cd(2);
        gr_rpower_10_100_xtrl->GetXaxis()->SetRangeUser(0.6, 1.01);
        gr_rpower_10_100_xtrl->GetYaxis()->SetRangeUser(10, 1e+5);
        gr_rpower_10_100_xtrl->Draw();
        gr_rpower_10_100_bdt->Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogy();

        print_canvas.cd(0);
        TPaveLabel label(0.0, 0.95, 0.3, 1, "BDT performances - 10 GeV - 100 GeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_performances.pdf(","Title:BDT performances - 10 GeV - 100 GeV");


        print_canvas.Clear();
        print_canvas.Divide(2, 1);

        TFile input_file_100_250(roc_100_250, "READ");
        if (input_file_100_250.IsZombie())
        {
            std::cerr << "\n\nError opening input file [" << roc_100_250 << "]\n\n";
            exit(100);
        }

        auto gr_roc_100_250_xtrl = static_cast<TGraph*>(input_file_100_250.Get("gr_xtrl_roc"));
        auto gr_roc_100_250_bdt = static_cast<TGraph*>(input_file_100_250.Get("gr_bdt_roc"));
        auto gr_rpower_100_250_xtrl = static_cast<TGraph*>(input_file_100_250.Get("gr_xtrl_proton_rejection_power"));
        auto gr_rpower_100_250_bdt = static_cast<TGraph*>(input_file_100_250.Get("gr_bdt_proton_rejection_power"));

        gr_roc_100_250_xtrl->SetLineWidth(2);
        gr_roc_100_250_xtrl->SetLineColor(kBlue);
        gr_roc_100_250_bdt->SetLineWidth(2);
        gr_roc_100_250_bdt->SetLineColor(kRed);
        gr_rpower_100_250_xtrl->SetLineWidth(2);
        gr_rpower_100_250_xtrl->SetLineColor(kBlue);
        gr_rpower_100_250_bdt->SetLineWidth(2);
        gr_rpower_100_250_bdt->SetLineColor(kRed);
        
        print_canvas.cd(1);
        gr_roc_100_250_xtrl->GetXaxis()->SetRangeUser(0, 1.05);
        gr_roc_100_250_xtrl->GetYaxis()->SetRangeUser(0.99, 1.001);
        gr_roc_100_250_xtrl->Draw();
        gr_roc_100_250_bdt->Draw("same");
        gPad->SetGrid(1,1);
        print_canvas.cd(2);
        gr_rpower_100_250_xtrl->GetXaxis()->SetRangeUser(0.6, 1.01);
        gr_rpower_100_250_xtrl->GetYaxis()->SetRangeUser(10, 1e+5);
        gr_rpower_100_250_xtrl->Draw();
        gr_rpower_100_250_bdt->Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogy();

        print_canvas.cd(0);
        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT performances - 100 GeV - 250 GeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_performances.pdf","Title:BDT performances - 100 GeV - 250 GeV");

        print_canvas.Clear();
        print_canvas.Divide(2, 1);

        TFile input_file_250_500(roc_250_500, "READ");
        if (input_file_250_500.IsZombie())
        {
            std::cerr << "\n\nError opening input file [" << roc_250_500 << "]\n\n";
            exit(100);
        }

        auto gr_roc_250_500_xtrl = static_cast<TGraph*>(input_file_250_500.Get("gr_xtrl_roc"));
        auto gr_roc_250_500_bdt = static_cast<TGraph*>(input_file_250_500.Get("gr_bdt_roc"));
        auto gr_rpower_250_500_xtrl = static_cast<TGraph*>(input_file_250_500.Get("gr_xtrl_proton_rejection_power"));
        auto gr_rpower_250_500_bdt = static_cast<TGraph*>(input_file_250_500.Get("gr_bdt_proton_rejection_power"));

        gr_roc_250_500_xtrl->SetLineWidth(2);
        gr_roc_250_500_xtrl->SetLineColor(kBlue);
        gr_roc_250_500_bdt->SetLineWidth(2);
        gr_roc_250_500_bdt->SetLineColor(kRed);
        gr_rpower_250_500_xtrl->SetLineWidth(2);
        gr_rpower_250_500_xtrl->SetLineColor(kBlue);
        gr_rpower_250_500_bdt->SetLineWidth(2);
        gr_rpower_250_500_bdt->SetLineColor(kRed);
        
        print_canvas.cd(1);
        gr_roc_250_500_xtrl->GetXaxis()->SetRangeUser(0, 1.05);
        gr_roc_250_500_xtrl->GetYaxis()->SetRangeUser(0.99, 1.001);
        gr_roc_250_500_xtrl->Draw();
        gr_roc_250_500_bdt->Draw("same");
        gPad->SetGrid(1,1);
        print_canvas.cd(2);
        gr_rpower_250_500_xtrl->GetXaxis()->SetRangeUser(0.6, 1.01);
        gr_rpower_250_500_xtrl->GetYaxis()->SetRangeUser(10, 1e+5);
        gr_rpower_250_500_xtrl->Draw();
        gr_rpower_250_500_bdt->Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogy();

        print_canvas.cd(0);
        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT performances - 250 GeV - 500 GeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_performances.pdf","Title:BDT performances - 250 GeV - 500 GeV");

        print_canvas.Clear();
        print_canvas.Divide(2, 1);

        TFile input_file_500_1000(roc_500_1000, "READ");
        if (input_file_500_1000.IsZombie())
        {
            std::cerr << "\n\nError opening input file [" << roc_500_1000 << "]\n\n";
            exit(100);
        }

        auto gr_roc_500_1000_xtrl = static_cast<TGraph*>(input_file_500_1000.Get("gr_xtrl_roc"));
        auto gr_roc_500_1000_bdt = static_cast<TGraph*>(input_file_500_1000.Get("gr_bdt_roc"));
        auto gr_rpower_500_1000_xtrl = static_cast<TGraph*>(input_file_500_1000.Get("gr_xtrl_proton_rejection_power"));
        auto gr_rpower_500_1000_bdt = static_cast<TGraph*>(input_file_500_1000.Get("gr_bdt_proton_rejection_power"));

        gr_roc_500_1000_xtrl->SetLineWidth(2);
        gr_roc_500_1000_xtrl->SetLineColor(kBlue);
        gr_roc_500_1000_bdt->SetLineWidth(2);
        gr_roc_500_1000_bdt->SetLineColor(kRed);
        gr_rpower_500_1000_xtrl->SetLineWidth(2);
        gr_rpower_500_1000_xtrl->SetLineColor(kBlue);
        gr_rpower_500_1000_bdt->SetLineWidth(2);
        gr_rpower_500_1000_bdt->SetLineColor(kRed);
        
        print_canvas.cd(1);
        gr_roc_500_1000_xtrl->GetXaxis()->SetRangeUser(0, 1.05);
        gr_roc_500_1000_xtrl->GetYaxis()->SetRangeUser(0.99, 1.001);
        gr_roc_500_1000_xtrl->Draw();
        gr_roc_500_1000_bdt->Draw("same");
        gPad->SetGrid(1,1);
        print_canvas.cd(2);
        gr_rpower_500_1000_xtrl->GetXaxis()->SetRangeUser(0.6, 1.01);
        gr_rpower_500_1000_xtrl->GetYaxis()->SetRangeUser(10, 1e+5);
        gr_rpower_500_1000_xtrl->Draw();
        gr_rpower_500_1000_bdt->Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogy();

        print_canvas.cd(0);
        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT performances - 500 GeV - 1 TeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_performances.pdf","Title:BDT performances - 500 GeV - 1 TeV");

        print_canvas.Clear();
        print_canvas.Divide(2, 1);

        TFile input_file_1000_3000(roc_1000_3000, "READ");
        if (input_file_1000_3000.IsZombie())
        {
            std::cerr << "\n\nError opening input file [" << roc_1000_3000 << "]\n\n";
            exit(100);
        }

        auto gr_roc_1000_3000_xtrl = static_cast<TGraph*>(input_file_1000_3000.Get("gr_xtrl_roc"));
        auto gr_roc_1000_3000_bdt = static_cast<TGraph*>(input_file_1000_3000.Get("gr_bdt_roc"));
        auto gr_rpower_1000_3000_xtrl = static_cast<TGraph*>(input_file_1000_3000.Get("gr_xtrl_proton_rejection_power"));
        auto gr_rpower_1000_3000_bdt = static_cast<TGraph*>(input_file_1000_3000.Get("gr_bdt_proton_rejection_power"));

        gr_roc_1000_3000_xtrl->SetLineWidth(2);
        gr_roc_1000_3000_xtrl->SetLineColor(kBlue);
        gr_roc_1000_3000_bdt->SetLineWidth(2);
        gr_roc_1000_3000_bdt->SetLineColor(kRed);
        gr_rpower_1000_3000_xtrl->SetLineWidth(2);
        gr_rpower_1000_3000_xtrl->SetLineColor(kBlue);
        gr_rpower_1000_3000_bdt->SetLineWidth(2);
        gr_rpower_1000_3000_bdt->SetLineColor(kRed);
        
        print_canvas.cd(1);
        gr_roc_1000_3000_xtrl->GetXaxis()->SetRangeUser(0, 1.05);
        gr_roc_1000_3000_xtrl->GetYaxis()->SetRangeUser(0.99, 1.001);
        gr_roc_1000_3000_xtrl->Draw();
        gr_roc_1000_3000_bdt->Draw("same");
        gPad->SetGrid(1,1);
        print_canvas.cd(2);
        gr_rpower_1000_3000_xtrl->GetXaxis()->SetRangeUser(0.6, 1.01);
        gr_rpower_1000_3000_xtrl->GetYaxis()->SetRangeUser(10, 1e+5);
        gr_rpower_1000_3000_xtrl->Draw();
        gr_rpower_1000_3000_bdt->Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogy();

        print_canvas.cd(0);
        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT performances - 1 TeV - 3 TeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_performances.pdf","Title:BDT performances - 1 TeV - 3 TeV");

        print_canvas.Clear();
        print_canvas.Divide(2, 1);

        TFile input_file_3000(roc_3000, "READ");
        if (input_file_3000.IsZombie())
        {
            std::cerr << "\n\nError opening input file [" << roc_3000 << "]\n\n";
            exit(100);
        }

        auto gr_roc_3000_xtrl = static_cast<TGraph*>(input_file_3000.Get("gr_xtrl_roc"));
        auto gr_roc_3000_bdt = static_cast<TGraph*>(input_file_3000.Get("gr_bdt_roc"));
        auto gr_rpower_3000_xtrl = static_cast<TGraph*>(input_file_3000.Get("gr_xtrl_proton_rejection_power"));
        auto gr_rpower_3000_bdt = static_cast<TGraph*>(input_file_3000.Get("gr_bdt_proton_rejection_power"));

        gr_roc_3000_xtrl->SetLineWidth(2);
        gr_roc_3000_xtrl->SetLineColor(kBlue);
        gr_roc_3000_bdt->SetLineWidth(2);
        gr_roc_3000_bdt->SetLineColor(kRed);
        gr_rpower_3000_xtrl->SetLineWidth(2);
        gr_rpower_3000_xtrl->SetLineColor(kBlue);
        gr_rpower_3000_bdt->SetLineWidth(2);
        gr_rpower_3000_bdt->SetLineColor(kRed);
        
        print_canvas.cd(1);
        gr_roc_3000_xtrl->GetXaxis()->SetRangeUser(0, 1.05);
        gr_roc_3000_xtrl->GetYaxis()->SetRangeUser(0.99, 1.001);
        gr_roc_3000_xtrl->Draw();
        gr_roc_3000_bdt->Draw("same");
        gPad->SetGrid(1,1);
        print_canvas.cd(2);
        gr_rpower_3000_xtrl->GetXaxis()->SetRangeUser(0.6, 1.01);
        gr_rpower_3000_xtrl->GetYaxis()->SetRangeUser(10, 1e+5);
        gr_rpower_3000_xtrl->Draw();
        gr_rpower_3000_bdt->Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogy();

        print_canvas.cd(0);
        label = TPaveLabel(0.0, 0.95, 0.3, 1, "BDT performances - 3 TeV - 10 TeV", "tlNDC");
        label.Draw();
        gStyle->SetOptTitle(0);
        print_canvas.Print("bdt_performances.pdf)","Title:BDT performances - 3 TeV - 10 TeV");
        

    }