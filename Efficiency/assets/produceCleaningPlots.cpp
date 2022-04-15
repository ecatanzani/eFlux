#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

#include "TF1.h"
#include "TH2D.h"
#include "TPDF.h"
#include "TPad.h"
#include "TLine.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TEfficiency.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"


void produceCleaningPlots(const char* input_file_name) {

    auto slicenorm = [](TH2D* histo) -> TH2D* {
        for (int bidx=1; bidx<=histo->GetNbinsX(); ++bidx)
        {
            auto scale_factor = histo->Integral(bidx, bidx, 1, histo->GetNbinsY()) ? 1/histo->Integral(bidx, bidx, 1, histo->GetNbinsY()) : 1;
            if (scale_factor!=1) {
                for (int bidy=1; bidy<=histo->GetNbinsY(); ++bidy) {
                    histo->SetBinContent(bidx, bidy, histo->GetBinContent(bidx, bidy)*scale_factor);
                    histo->SetBinError(bidx, bidy, histo->GetBinError(bidx, bidy)*scale_factor);
                }
            }
        }
        return histo; 
    };

    TFile *input_file = TFile::Open(input_file_name, "READ");
    if (input_file->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
        exit(100);
    }

    auto h_stk_cleaning_20_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_20_100"));
    auto h_stk_cleaning_100_250 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_100_250"));
    auto h_stk_cleaning_250_500 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_250_500"));
    auto h_stk_cleaning_500_1000 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_500_1000"));
    auto h_stk_cleaning_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000"));
    auto h_stk_cleaning_3000 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_3000"));

    auto h_stk_cleaning_20_100_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12"));
    auto h_stk_cleaning_100_250_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12"));
    auto h_stk_cleaning_250_500_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12"));
    auto h_stk_cleaning_500_1000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12"));
    auto h_stk_cleaning_1000_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12"));
    auto h_stk_cleaning_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12"));

    auto h_stk_cleaning_20_100_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12_100"));
    auto h_stk_cleaning_100_250_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12_100"));
    auto h_stk_cleaning_250_500_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12_100"));
    auto h_stk_cleaning_500_1000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12_100"));
    auto h_stk_cleaning_1000_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12_100"));
    auto h_stk_cleaning_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12_100"));

    auto h_rvalue = slicenorm(static_cast<TH2D*>(input_file->Get("rvalue/h_rvalue")));
    auto h_rvalue_xtrl_12 = slicenorm(static_cast<TH2D*>(input_file->Get("rvalue/h_rvalue_xtrl_12")));
    auto h_rvalue_xtrl_12_100 = slicenorm(static_cast<TH2D*>(input_file->Get("rvalue/h_rvalue_xtrl_12_100")));

    auto h_lvalue = slicenorm(static_cast<TH2D*>(input_file->Get("lvalue/h_lvalue")));
    auto h_lvalue_xtrl_12 = slicenorm(static_cast<TH2D*>(input_file->Get("lvalue/h_lvalue_xtrl_12")));
    auto h_lvalue_xtrl_12_100 = slicenorm(static_cast<TH2D*>(input_file->Get("lvalue/h_lvalue_xtrl_12_100")));

    h_stk_cleaning_20_100->SetDirectory(0);
    h_stk_cleaning_100_250->SetDirectory(0);
    h_stk_cleaning_250_500->SetDirectory(0);
    h_stk_cleaning_500_1000->SetDirectory(0);
    h_stk_cleaning_1000_3000->SetDirectory(0);
    h_stk_cleaning_3000->SetDirectory(0);

    h_stk_cleaning_20_100_xtrl_12->SetDirectory(0);
    h_stk_cleaning_100_250_xtrl_12->SetDirectory(0);
    h_stk_cleaning_250_500_xtrl_12->SetDirectory(0);
    h_stk_cleaning_500_1000_xtrl_12->SetDirectory(0);
    h_stk_cleaning_1000_3000_xtrl_12->SetDirectory(0);
    h_stk_cleaning_3000_xtrl_12->SetDirectory(0);

    h_stk_cleaning_20_100_xtrl_12_100->SetDirectory(0);
    h_stk_cleaning_100_250_xtrl_12_100->SetDirectory(0);
    h_stk_cleaning_250_500_xtrl_12_100->SetDirectory(0);
    h_stk_cleaning_500_1000_xtrl_12_100->SetDirectory(0);
    h_stk_cleaning_1000_3000_xtrl_12_100->SetDirectory(0);
    h_stk_cleaning_3000_xtrl_12_100->SetDirectory(0);

    h_rvalue->SetDirectory(0);
    h_rvalue_xtrl_12->SetDirectory(0);
    h_rvalue_xtrl_12_100->SetDirectory(0);

    h_lvalue->SetDirectory(0);
    h_lvalue_xtrl_12->SetDirectory(0);
    h_lvalue_xtrl_12_100->SetDirectory(0);

    input_file->Close(); 

    TCanvas print_canvas("print_canvas", "print_canvas");
    print_canvas.Divide(2, 3);

    print_canvas.SetTicks();
    
    TPaveLabel label(0.0, 0.97, 0.3, 1, "STK cleaning cuts", "tlNDC");
    label.Draw();

    print_canvas.cd(1);

    h_stk_cleaning_20_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(2);

    h_stk_cleaning_100_250->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(3);

    h_stk_cleaning_250_500->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(4);

    h_stk_cleaning_500_1000->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(5);

    h_stk_cleaning_1000_3000->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(6);

    h_stk_cleaning_3000->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.Print("stk_cleaning_cuts.pdf(","Title:STK cleaning cuts");

    label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts (xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);

    h_stk_cleaning_20_100_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(2);

    h_stk_cleaning_100_250_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(3);

    h_stk_cleaning_250_500_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(4);

    h_stk_cleaning_500_1000_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(5);

    h_stk_cleaning_1000_3000_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(6);

    h_stk_cleaning_3000_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.Print("stk_cleaning_cuts.pdf","Title:STK cleaning cuts (xtrl < 12)");

    label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts (12 < xtrl < 100)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);

    h_stk_cleaning_20_100_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(2);

    h_stk_cleaning_100_250_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(3);

    h_stk_cleaning_250_500_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(4);

    h_stk_cleaning_500_1000_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(5);

    h_stk_cleaning_1000_3000_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(6);

    h_stk_cleaning_3000_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.Print("stk_cleaning_cuts.pdf)","Title:STK cleaning cuts (12 < xtrl < 100)");

    TLine l1(0, 1, 1e+4, 1);
    TLine l2(0, 6, 1e+4, 6);

    l1.SetLineWidth(2);
    l1.SetLineStyle(2);
    l1.SetLineColor(kMagenta);

    l2.SetLineWidth(2);
    l2.SetLineStyle(2);
    l2.SetLineColor(kMagenta);

    TCanvas print_canvas_lr ("print_canvas_lr", "print_canvas_lr");

    h_rvalue->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    l1.Draw("same");
    l2.Draw("same");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "R Value", "tlNDC");
    label.Draw();

    print_canvas_lr.Print("rl_cleaning_cuts.pdf(","Title:RValue");

    h_rvalue_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    l1.Draw("same");
    l2.Draw("same");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "R Value (xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas_lr.Print("rl_cleaning_cuts.pdf","Title:RValue (xtrl < 12)");
    
    h_rvalue_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    l1.Draw("same");
    l2.Draw("same");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "R Value (12 < xtrl < 100)", "tlNDC");
    label.Draw();

    print_canvas_lr.Print("rl_cleaning_cuts.pdf","Title:RValue (12 < xtrl < 100)");

    TF1 l_lower_limit("l_lower_limit", "0.1", 0, 1e+4);
    TF1 l_upper_limit_under_1TeV("l_upper_limit_under_1TeV","0.059 * (log10(x) - 3) + 0.5", 10, 1e+3);
    TF1 l_upper_limit_above_1TeV("l_upper_limit_under_1TeV","0.038 * (log10(x) - 3) + 0.5", 1e+3, 1e+4);

    l_lower_limit.SetLineWidth(2);
    l_lower_limit.SetLineStyle(2);
    l_lower_limit.SetLineColor(kMagenta);

    l_upper_limit_under_1TeV.SetLineWidth(2);
    l_upper_limit_under_1TeV.SetLineStyle(2);
    l_upper_limit_under_1TeV.SetLineColor(kMagenta);

    l_upper_limit_above_1TeV.SetLineWidth(2);
    l_upper_limit_above_1TeV.SetLineStyle(2);
    l_upper_limit_above_1TeV.SetLineColor(kMagenta);

    h_lvalue->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    l_lower_limit.Draw("same");
    l_upper_limit_under_1TeV.Draw("same");
    l_upper_limit_above_1TeV.Draw("same");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "L Value", "tlNDC");
    label.Draw();

    print_canvas_lr.Print("rl_cleaning_cuts.pdf(","Title:RValue");

    h_lvalue_xtrl_12->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    l_lower_limit.Draw("same");
    l_upper_limit_under_1TeV.Draw("same");
    l_upper_limit_above_1TeV.Draw("same");

    label = TPaveLabel(0.0, 0.97, 0.3, 1, "L Value (xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas_lr.Print("rl_cleaning_cuts.pdf","Title:RValue (xtrl < 12)");

    h_lvalue_xtrl_12_100->Draw("colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    l_lower_limit.Draw("same");
    l_upper_limit_under_1TeV.Draw("same");
    l_upper_limit_above_1TeV.Draw("same");

    label = TPaveLabel(0.0, 0.97, 0.3, 1, "L Value (12 < xtrl < 100)", "tlNDC");
    label.Draw();

    print_canvas_lr.Print("rl_cleaning_cuts.pdf)","Title:RValue (12 < xtrl < 100)");

}

void produceStkCleaningPlotsProjections(
    const char* data_file,
    const char* mc_file,
    const char* data_label="DATA",
    const char* mc_label="electron MC") {

        TFile* datafile = TFile::Open(data_file, "READ");
        if (datafile->IsZombie()) {
            std::cerr << "\n\nError opening input DATA ROOT file [" << data_file << "]\n\n";
            exit(100);
        }

        auto data_h_stk_cleaning_20_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_20_100"));
        auto data_h_stk_cleaning_100_250 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_100_250"));
        auto data_h_stk_cleaning_250_500 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_250_500"));
        auto data_h_stk_cleaning_500_1000 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000"));
        auto data_h_stk_cleaning_1000_3000 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000"));
        auto data_h_stk_cleaning_3000 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_3000"));

        auto data_h_stk_cleaning_20_100_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12"));
        auto data_h_stk_cleaning_100_250_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12"));
        auto data_h_stk_cleaning_250_500_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12"));
        auto data_h_stk_cleaning_500_1000_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12"));
        auto data_h_stk_cleaning_1000_3000_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12"));
        auto data_h_stk_cleaning_3000_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12"));

        auto data_h_stk_cleaning_20_100_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12_100"));
        auto data_h_stk_cleaning_100_250_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12_100"));
        auto data_h_stk_cleaning_250_500_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12_100"));
        auto data_h_stk_cleaning_500_1000_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12_100"));
        auto data_h_stk_cleaning_1000_3000_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12_100"));
        auto data_h_stk_cleaning_3000_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12_100"));

        data_h_stk_cleaning_20_100->SetTitle(data_label);
        data_h_stk_cleaning_100_250->SetTitle(data_label);
        data_h_stk_cleaning_250_500->SetTitle(data_label);
        data_h_stk_cleaning_500_1000->SetTitle(data_label);
        data_h_stk_cleaning_1000_3000->SetTitle(data_label);
        data_h_stk_cleaning_3000->SetTitle(data_label);

        data_h_stk_cleaning_20_100_xtrl_12->SetTitle(data_label);
        data_h_stk_cleaning_100_250_xtrl_12->SetTitle(data_label);
        data_h_stk_cleaning_250_500_xtrl_12->SetTitle(data_label);
        data_h_stk_cleaning_500_1000_xtrl_12->SetTitle(data_label);
        data_h_stk_cleaning_1000_3000_xtrl_12->SetTitle(data_label);
        data_h_stk_cleaning_3000_xtrl_12->SetTitle(data_label);

        data_h_stk_cleaning_20_100_xtrl_12_100->SetTitle(data_label);
        data_h_stk_cleaning_100_250_xtrl_12_100->SetTitle(data_label);
        data_h_stk_cleaning_250_500_xtrl_12_100->SetTitle(data_label);
        data_h_stk_cleaning_500_1000_xtrl_12_100->SetTitle(data_label);
        data_h_stk_cleaning_1000_3000_xtrl_12_100->SetTitle(data_label);
        data_h_stk_cleaning_3000_xtrl_12_100->SetTitle(data_label);

        data_h_stk_cleaning_20_100->SetDirectory(0);
        data_h_stk_cleaning_100_250->SetDirectory(0);
        data_h_stk_cleaning_250_500->SetDirectory(0);
        data_h_stk_cleaning_500_1000->SetDirectory(0);
        data_h_stk_cleaning_1000_3000->SetDirectory(0);
        data_h_stk_cleaning_3000->SetDirectory(0);

        data_h_stk_cleaning_20_100_xtrl_12->SetDirectory(0);
        data_h_stk_cleaning_100_250_xtrl_12->SetDirectory(0);
        data_h_stk_cleaning_250_500_xtrl_12->SetDirectory(0);
        data_h_stk_cleaning_500_1000_xtrl_12->SetDirectory(0);
        data_h_stk_cleaning_1000_3000_xtrl_12->SetDirectory(0);
        data_h_stk_cleaning_3000_xtrl_12->SetDirectory(0);

        data_h_stk_cleaning_20_100_xtrl_12_100->SetDirectory(0);
        data_h_stk_cleaning_100_250_xtrl_12_100->SetDirectory(0);
        data_h_stk_cleaning_250_500_xtrl_12_100->SetDirectory(0);
        data_h_stk_cleaning_500_1000_xtrl_12_100->SetDirectory(0);
        data_h_stk_cleaning_1000_3000_xtrl_12_100->SetDirectory(0);
        data_h_stk_cleaning_3000_xtrl_12_100->SetDirectory(0);

        datafile->Close();

        TFile* mcfile = TFile::Open(mc_file, "READ");
        if (mcfile->IsZombie()) {
            std::cerr << "\n\nError opening input MC ROOT file [" << mc_file << "]\n\n";
            exit(100);
        }

        auto mc_h_stk_cleaning_20_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_20_100"));
        auto mc_h_stk_cleaning_100_250 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_100_250"));
        auto mc_h_stk_cleaning_250_500 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_250_500"));
        auto mc_h_stk_cleaning_500_1000 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000"));
        auto mc_h_stk_cleaning_1000_3000 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000"));
        auto mc_h_stk_cleaning_3000 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_3000"));

        auto mc_h_stk_cleaning_20_100_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12"));
        auto mc_h_stk_cleaning_100_250_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12"));
        auto mc_h_stk_cleaning_250_500_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12"));
        auto mc_h_stk_cleaning_500_1000_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12"));
        auto mc_h_stk_cleaning_1000_3000_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12"));
        auto mc_h_stk_cleaning_3000_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12"));

        auto mc_h_stk_cleaning_20_100_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12_100"));
        auto mc_h_stk_cleaning_100_250_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12_100"));
        auto mc_h_stk_cleaning_250_500_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12_100"));
        auto mc_h_stk_cleaning_500_1000_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12_100"));
        auto mc_h_stk_cleaning_1000_3000_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12_100"));
        auto mc_h_stk_cleaning_3000_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12_100"));

        mc_h_stk_cleaning_20_100->SetTitle(mc_label);
        mc_h_stk_cleaning_100_250->SetTitle(mc_label);
        mc_h_stk_cleaning_250_500->SetTitle(mc_label);
        mc_h_stk_cleaning_500_1000->SetTitle(mc_label);
        mc_h_stk_cleaning_1000_3000->SetTitle(mc_label);
        mc_h_stk_cleaning_3000->SetTitle(mc_label);

        mc_h_stk_cleaning_20_100_xtrl_12->SetTitle(mc_label);
        mc_h_stk_cleaning_100_250_xtrl_12->SetTitle(mc_label);
        mc_h_stk_cleaning_250_500_xtrl_12->SetTitle(mc_label);
        mc_h_stk_cleaning_500_1000_xtrl_12->SetTitle(mc_label);
        mc_h_stk_cleaning_1000_3000_xtrl_12->SetTitle(mc_label);
        mc_h_stk_cleaning_3000_xtrl_12->SetTitle(mc_label);

        mc_h_stk_cleaning_20_100_xtrl_12_100->SetTitle(mc_label);
        mc_h_stk_cleaning_100_250_xtrl_12_100->SetTitle(mc_label);
        mc_h_stk_cleaning_250_500_xtrl_12_100->SetTitle(mc_label);
        mc_h_stk_cleaning_500_1000_xtrl_12_100->SetTitle(mc_label);
        mc_h_stk_cleaning_1000_3000_xtrl_12_100->SetTitle(mc_label);
        mc_h_stk_cleaning_3000_xtrl_12_100->SetTitle(mc_label);

        mc_h_stk_cleaning_20_100->SetDirectory(0);
        mc_h_stk_cleaning_100_250->SetDirectory(0);
        mc_h_stk_cleaning_250_500->SetDirectory(0);
        mc_h_stk_cleaning_500_1000->SetDirectory(0);
        mc_h_stk_cleaning_1000_3000->SetDirectory(0);
        mc_h_stk_cleaning_3000->SetDirectory(0);

        mc_h_stk_cleaning_20_100_xtrl_12->SetDirectory(0);
        mc_h_stk_cleaning_100_250_xtrl_12->SetDirectory(0);
        mc_h_stk_cleaning_250_500_xtrl_12->SetDirectory(0);
        mc_h_stk_cleaning_500_1000_xtrl_12->SetDirectory(0);
        mc_h_stk_cleaning_1000_3000_xtrl_12->SetDirectory(0);
        mc_h_stk_cleaning_3000_xtrl_12->SetDirectory(0);

        mc_h_stk_cleaning_20_100_xtrl_12_100->SetDirectory(0);
        mc_h_stk_cleaning_100_250_xtrl_12_100->SetDirectory(0);
        mc_h_stk_cleaning_250_500_xtrl_12_100->SetDirectory(0);
        mc_h_stk_cleaning_500_1000_xtrl_12_100->SetDirectory(0);
        mc_h_stk_cleaning_1000_3000_xtrl_12_100->SetDirectory(0);
        mc_h_stk_cleaning_3000_xtrl_12_100->SetDirectory(0);

        mcfile->Close();

        int nbinsX {data_h_stk_cleaning_20_100->GetNbinsX()};

        std::vector<double> nstk_cuts_20_100;
        std::vector<double> nstk_cuts_100_250;
        std::vector<double> nstk_cuts_250_500;
        std::vector<double> nstk_cuts_500_1000;
        std::vector<double> nstk_cuts_1000_3000;
        std::vector<double> nstk_cuts_3000;

        std::vector<double> nstk_20_100;
        std::vector<double> nstk_100_250;
        std::vector<double> nstk_250_500;
        std::vector<double> nstk_500_1000;
        std::vector<double> nstk_1000_3000;
        std::vector<double> nstk_3000;

        std::vector<double> cat_record (nbinsX, -999);

        // Projections of the STK cleaning cuts for the energy range 20-100 GeV
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_20_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_20_100->ProjectionY("mc_proj", bidx, bidx));

            bool max_found {false};
            if (mc_proj->Integral())
            {
                double event_threshold {mc_proj->Integral()*(1-1e-3)};
                for (int idx=mc_proj->GetNbinsX(); idx>=1; --idx)
                {
                    if (mc_proj->Integral(1, idx) <= event_threshold)
                    {
                        nstk_cuts_20_100.push_back(mc_proj->GetBinCenter(idx+2));
                        nstk_20_100.push_back(data_h_stk_cleaning_20_100->GetXaxis()->GetBinCenter(bidx));
                        max_found = true;
                        break;
                    }
                }
            }

            if (max_found)
                cat_record[bidx] = nstk_cuts_20_100.back();

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");
            
            std::unique_ptr<TLine> line_cut;
            if (max_found)
            {
                line_cut = std::make_unique<TLine>(nstk_cuts_20_100.back(), 0, nstk_cuts_20_100.back(), max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters bin ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_20_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_20_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_20_100.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 20-100 GeV (xtrl < 12)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_20_100_xtrl_12->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_20_100_xtrl_12->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_20_100_xtrl_12.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_20_100_xtrl_12.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_20_100_xtrl_12.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 20-100 GeV ( 12 < xtrl < 100)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_20_100_xtrl_12_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_20_100_xtrl_12_100->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);
            
            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_20_100_xtrl_12_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_20_100_xtrl_12_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_20_100_xtrl_12_100.pdf)","Title:STK cleaning cuts");
        }

        cat_record = std::vector<double> (nbinsX, -999);

        // Projections of the STK cleaning cuts for the energy range 100-250 GeV
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_100_250->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_100_250->ProjectionY("mc_proj", bidx, bidx));

            bool max_found {false};
            if (mc_proj->Integral())
            {
                double event_threshold {mc_proj->Integral()*(1-1e-3)};
                for (int idx=mc_proj->GetNbinsX(); idx>=1; --idx)
                {
                    if (mc_proj->Integral(1, idx) <= event_threshold)
                    {
                        nstk_cuts_100_250.push_back(mc_proj->GetBinCenter(idx+2));
                        nstk_100_250.push_back(data_h_stk_cleaning_20_100->GetXaxis()->GetBinCenter(bidx));
                        max_found = true;
                        break;
                    }
                }
            }

            if (max_found)
                cat_record[bidx] = nstk_cuts_100_250.back();

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (max_found)
            {
                line_cut = std::make_unique<TLine>(nstk_cuts_100_250.back(), 0, nstk_cuts_100_250.back(), max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_100_250.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_100_250.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_100_250.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 100-250 GeV (xtrl < 12)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_100_250_xtrl_12->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_100_250_xtrl_12->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_100_250_xtrl_12.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_100_250_xtrl_12.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_100_250_xtrl_12.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 100-250 GeV (12 < xtrl < 100)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_100_250_xtrl_12_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_100_250_xtrl_12_100->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_100_250_xtrl_12_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_100_250_xtrl_12_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_100_250_xtrl_12_100.pdf)","Title:STK cleaning cuts");
        }

        cat_record = std::vector<double> (nbinsX, -999);

        // Projections of the STK cleaning cuts for the energy range 250-500 GeV
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_250_500->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_250_500->ProjectionY("mc_proj", bidx, bidx));

            bool max_found {false};
            if (mc_proj->Integral())
            {
                double event_threshold {mc_proj->Integral()*(1-1e-3)};
                for (int idx=mc_proj->GetNbinsX(); idx>=1; --idx)
                {
                    if (mc_proj->Integral(1, idx) <= event_threshold)
                    {
                        nstk_cuts_250_500.push_back(mc_proj->GetBinCenter(idx+2));
                        nstk_250_500.push_back(data_h_stk_cleaning_20_100->GetXaxis()->GetBinCenter(bidx));
                        max_found = true;
                        break;
                    }
                }
            }

            if (max_found)
                cat_record[bidx] = nstk_cuts_250_500.back();

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (max_found)
            {
                line_cut = std::make_unique<TLine>(nstk_cuts_250_500.back(), 0, nstk_cuts_250_500.back(), max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_250_500.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_250_500.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_250_500.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 250-500 GeV (xtrl < 12)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_250_500_xtrl_12->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_250_500_xtrl_12->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_250_500_xtrl_12.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_250_500_xtrl_12.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_250_500_xtrl_12.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 250-500 GeV (12 < xtrl < 100)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_250_500_xtrl_12_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_250_500_xtrl_12_100->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_250_500_xtrl_12_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_250_500_xtrl_12_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_250_500_xtrl_12_100.pdf)","Title:STK cleaning cuts");
        }

        cat_record = std::vector<double> (nbinsX, -999);

        // Projections of the STK cleaning cuts for the energy range 500 GeV -1 TeV
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_500_1000->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_500_1000->ProjectionY("mc_proj", bidx, bidx));

            bool max_found {false};
            if (mc_proj->Integral())
            {
                double event_threshold {mc_proj->Integral()*(1-1e-3)};
                for (int idx=mc_proj->GetNbinsX(); idx>=1; --idx)
                {
                    if (mc_proj->Integral(1, idx) <= event_threshold)
                    {
                        nstk_cuts_500_1000.push_back(mc_proj->GetBinCenter(idx+2));
                        nstk_500_1000.push_back(data_h_stk_cleaning_20_100->GetXaxis()->GetBinCenter(bidx));
                        max_found = true;
                        break;
                    }
                }
            }

            if (max_found)
                cat_record[bidx] = nstk_cuts_500_1000.back();

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (max_found)
            {
                line_cut = std::make_unique<TLine>(nstk_cuts_500_1000.back(), 0, nstk_cuts_500_1000.back(), max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 500 GeV -1 TeV (xtrl < 12)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_500_1000_xtrl_12->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_500_1000_xtrl_12->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000_xtrl_12.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000_xtrl_12.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000_xtrl_12.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 500 GeV -1 TeV (12 < xtrl < 100)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_500_1000_xtrl_12_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_500_1000_xtrl_12_100->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000_xtrl_12_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000_xtrl_12_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_500_1000_xtrl_12_100.pdf)","Title:STK cleaning cuts");
        }
        
        cat_record = std::vector<double> (nbinsX, -999);

        // Projections of the STK cleaning cuts for the energy range 1 TeV - 3 TeV
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_1000_3000->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_1000_3000->ProjectionY("mc_proj", bidx, bidx));

            bool max_found {false};
            if (mc_proj->Integral())
            {
                double event_threshold {mc_proj->Integral()*(1-1e-3)};
                for (int idx=mc_proj->GetNbinsX(); idx>=1; --idx)
                {
                    if (mc_proj->Integral(1, idx) <= event_threshold)
                    {
                        nstk_cuts_1000_3000.push_back(mc_proj->GetBinCenter(idx+4));
                        nstk_1000_3000.push_back(data_h_stk_cleaning_20_100->GetXaxis()->GetBinCenter(bidx));
                        max_found = true;
                        break;
                    }
                }
            }

            if (max_found)
                cat_record[bidx] = nstk_cuts_1000_3000.back();

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (max_found)
            {
                line_cut = std::make_unique<TLine>(nstk_cuts_1000_3000.back(), 0, nstk_cuts_1000_3000.back(), max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 1 TeV - 3 TeV (xtrl < 12)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_1000_3000_xtrl_12->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_1000_3000_xtrl_12->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000_xtrl_12.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000_xtrl_12.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000_xtrl_12.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range 1 TeV - 3 TeV (12 < xtrl < 100)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_1000_3000_xtrl_12_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_1000_3000_xtrl_12_100->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000_xtrl_12_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000_xtrl_12_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_1000_3000_xtrl_12_100.pdf)","Title:STK cleaning cuts");
        }

        cat_record = std::vector<double> (nbinsX, -999);

        // Projections of the STK cleaning cuts for the energy range > 3 TeV
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_3000->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_3000->ProjectionY("mc_proj", bidx, bidx));

            bool max_found {false};
            if (mc_proj->Integral())
            {
                double event_threshold {mc_proj->Integral()*(1-1e-3)};
                for (int idx=mc_proj->GetNbinsX(); idx>=1; --idx)
                {
                    if (mc_proj->Integral(1, idx) <= event_threshold)
                    {
                        nstk_cuts_3000.push_back(mc_proj->GetBinCenter(idx+2));
                        nstk_3000.push_back(data_h_stk_cleaning_20_100->GetXaxis()->GetBinCenter(bidx));
                        max_found = true;
                        break;
                    }
                }
            }

            if (max_found)
                cat_record[bidx] = nstk_cuts_3000.back();

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (max_found)
            {
                line_cut = std::make_unique<TLine>(nstk_cuts_3000.back(), 0, nstk_cuts_3000.back(), max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_3000.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_3000.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_3000.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range > 3 TeV (xtrl < 12)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_3000_xtrl_12->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_3000_xtrl_12->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_3000_xtrl_12.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_3000_xtrl_12.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_3000_xtrl_12.pdf)","Title:STK cleaning cuts");
        }

        // Projections of the STK cleaning cuts for the energy range > 3 TeV (12 < xtrl < 100)
        for (int bidx=1; bidx<=nbinsX; ++bidx)
        {
            TCanvas print_canvas("print_canvas", "print_canvas");

            auto data_proj = static_cast<TH1D*>(data_h_stk_cleaning_3000_xtrl_12_100->ProjectionY("data_proj", bidx, bidx));
            auto mc_proj = static_cast<TH1D*>(mc_h_stk_cleaning_3000_xtrl_12_100->ProjectionY("mc_proj", bidx, bidx));

            if (data_proj->GetEntries())
                data_proj->Scale(1./data_proj->GetEntries());
            
            if (mc_proj->GetEntries())
                mc_proj->Scale(1./mc_proj->GetEntries());

            std::vector<double> data_bin_content (data_proj->GetNbinsX());
            std::vector<double> mc_bin_content (mc_proj->GetNbinsX());

            for (int bidx=1; bidx<=data_proj->GetNbinsX(); ++bidx)
            {
                data_bin_content[bidx-1] = static_cast<double>(data_proj->GetBinContent(bidx));
                mc_bin_content[bidx-1] = static_cast<double>(mc_proj->GetBinContent(bidx));
            }

            auto max_pad_data  = max_element(std::begin(data_bin_content), std::end(data_bin_content));
            auto max_pad_mc = max_element(std::begin(mc_bin_content), std::end(mc_bin_content));

            auto max_pad = std::max(*max_pad_data, *max_pad_mc);
            max_pad += 0.1*max_pad;

            data_proj->SetLineColor(kRed);
            mc_proj->SetLineColor(kBlue);

            data_proj->SetLineWidth(2);
            mc_proj->SetLineWidth(2);
            
            data_proj->SetStats(0);
            mc_proj->SetStats(0);

            data_proj->GetYaxis()->SetRangeUser(0, max_pad);
            mc_proj->GetYaxis()->SetRangeUser(0, max_pad);

            data_proj->Draw("hist");
            mc_proj->Draw("hist, same");

            std::unique_ptr<TLine> line_cut;
            if (cat_record[bidx]!=-999)
            {
                line_cut = std::make_unique<TLine>(cat_record[bidx], 0, cat_record[bidx], max_pad);
                line_cut->SetLineColor(kMagenta);
                line_cut->SetLineWidth(2);
                line_cut->SetLineStyle(2);
                line_cut->Draw("same");
            }

            gPad->SetLogx();

            gStyle->SetOptTitle(0);
            
            auto label_title = std::string("STK cleaning cuts - projection nSTK clusters: ") + std::to_string(bidx);
            TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
            label.Draw();
            
            if (mc_proj->GetEntries() || data_proj->GetEntries())
            {
                auto legend = print_canvas.BuildLegend();
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
                {
                    auto primitive = (TLegendEntry*)primitiveObj;
                    primitive->SetOption("l");
                }
            }

            gPad->SetGrid(1,1);

            if (bidx==1)
                print_canvas.Print("stk_cleaning_cuts_projections_3000_xtrl_12_100.pdf(","Title:STK cleaning cuts");
            else if (bidx<nbinsX)
                print_canvas.Print("stk_cleaning_cuts_projections_3000_xtrl_12_100.pdf","Title:STK cleaning cuts");
            else
                print_canvas.Print("stk_cleaning_cuts_projections_3000_xtrl_12_100.pdf)","Title:STK cleaning cuts");
        }


        // nSTK cluster cuts
        TGraph gr_cuts_20_100 (static_cast<int>(nstk_20_100.size()), &nstk_20_100[0], &nstk_cuts_20_100[0]);
        TGraph gr_cuts_100_250 (static_cast<int>(nstk_100_250.size()), &nstk_100_250[0], &nstk_cuts_100_250[0]);
        TGraph gr_cuts_250_500 (static_cast<int>(nstk_250_500.size()), &nstk_250_500[0], &nstk_cuts_250_500[0]);
        TGraph gr_cuts_500_1000 (static_cast<int>(nstk_500_1000.size()), &nstk_500_1000[0], &nstk_cuts_500_1000[0]);
        TGraph gr_cuts_1000_3000 (static_cast<int>(nstk_1000_3000.size()), &nstk_1000_3000[0], &nstk_cuts_1000_3000[0]);
        //TGraph gr_cuts_3000 (static_cast<int>(nstk_3000.size()), &nstk_3000[0], &nstk_cuts_3000[0]);
        
        TF1 f_cuts_20_100 ("f_cuts_20_100", "[0] + [1]*x", 5, 300);
        TF1 f_cuts_100_250 ("f_cuts_100_250", "[0] + [1]*x", 5, 300);
        TF1 f_cuts_250_500 ("f_cuts_250_500", "[0] + [1]*x", 5, 300);
        TF1 f_cuts_500_1000 ("f_cuts_500_1000", "[0] + [1]*x", 5, 300);
        TF1 f_cuts_1000_3000 ("f_cuts_1000_3000", "[0] + [1]*x", 5, 300);
        TF1 f_cuts_3000 ("f_cuts_3000", "[0] + [1]*x", 5, 300);

        /*
        TF1 f_cuts_20_100 ("f_cuts_20_100", "[0] + [1]*x", nstk_20_100.front(), nstk_20_100.back());
        TF1 f_cuts_100_250 ("f_cuts_100_250", "[0] + [1]*x", nstk_100_250.front(), nstk_100_250.back());
        TF1 f_cuts_250_500 ("f_cuts_250_500", "[0] + [1]*x", nstk_250_500.front(), nstk_250_500.back());
        TF1 f_cuts_500_1000 ("f_cuts_500_1000", "[0] + [1]*x", nstk_500_1000.front(), nstk_500_1000.back());
        TF1 f_cuts_1000_3000 ("f_cuts_1000_3000", "[0] + [1]*x", nstk_1000_3000.front(), nstk_1000_3000.back());
        */

        //TF1 f_cuts_3000 ("f_cuts_3000", "TMath::Log10([0] + [1]*TMath::Log10(x+1))", 0, 10);

        gr_cuts_20_100.Fit(&f_cuts_20_100, "Q", "", 10, 100);
        gr_cuts_100_250.Fit(&f_cuts_100_250, "Q", "", 10, 100);
        gr_cuts_250_500.Fit(&f_cuts_250_500, "Q", "", 10, 100);
        gr_cuts_500_1000.Fit(&f_cuts_500_1000, "Q", "", 10, 100);
        gr_cuts_1000_3000.Fit(&f_cuts_1000_3000, "Q", "", 5, nstk_1000_3000.back());
        //gr_cuts_3000.Fit(&f_cuts_3000, "R");

        gr_cuts_20_100.SetLineWidth(2);
        gr_cuts_100_250.SetLineWidth(2);
        gr_cuts_250_500.SetLineWidth(2);
        gr_cuts_500_1000.SetLineWidth(2);
        gr_cuts_1000_3000.SetLineWidth(2);
        //gr_cuts_3000.SetLineWidth(2);

        gr_cuts_20_100.SetMarkerStyle(105);
        gr_cuts_100_250.SetMarkerStyle(105);
        gr_cuts_250_500.SetMarkerStyle(105);
        gr_cuts_500_1000.SetMarkerStyle(105);
        gr_cuts_1000_3000.SetMarkerStyle(105);
        //gr_cuts_3000.SetMarkerStyle(105);

        gr_cuts_20_100.SetName("gr_cuts_20_100");
        gr_cuts_100_250.SetName("gr_cuts_100_250");
        gr_cuts_250_500.SetName("gr_cuts_250_500");
        gr_cuts_500_1000.SetName("gr_cuts_500_1000");
        gr_cuts_1000_3000.SetName("gr_cuts_1000_3000");
        //gr_cuts_3000.SetName("gr_cuts_3000");

        gr_cuts_20_100.SetTitle("nSTK cuts - 20 - 100 GeV; nStkClu1Rm; StkEcore1m");
        gr_cuts_100_250.SetTitle("nSTK cuts - 100 - 250 GeV; nStkClu1Rm; StkEcore1m");
        gr_cuts_250_500.SetTitle("nSTK cuts - 250 - 500 GeV; nStkClu1Rm; StkEcore1m");
        gr_cuts_500_1000.SetTitle("nSTK cuts - 500 GeV - 1 TeV; nStkClu1Rm; StkEcore1m");
        gr_cuts_1000_3000.SetTitle("nSTK cuts - 1 - 3 TeV; nStkClu1Rm; StkEcore1m");
        //gr_cuts_3000.SetTitle("nSTK cuts - > 3 TeV; nStkClu1Rm; StkEcore1m");

        TCanvas print_canvas("print_canvas", "print_canvas");

        gr_cuts_20_100.Draw();
        gPad->SetLogx();
        gPad->SetLogy();
        auto label_title = std::string("STK cleaning cuts - 20 GeV - 100 GeV cleaning cut");
        TPaveLabel label(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf(","Title:STK cleaning cuts");

        print_canvas.Clear();
        gr_cuts_100_250.Draw();
        gPad->SetLogx();
        gPad->SetLogy();
        label_title = std::string("STK cleaning cuts - 100 GeV - 250 GeV cleaning cut");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf","Title:STK cleaning cuts");

        print_canvas.Clear();
        gr_cuts_250_500.Draw();
        gPad->SetLogx();
        gPad->SetLogy();
        label_title = std::string("STK cleaning cuts - 250 GeV - 500 GeV cleaning cut");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf","Title:STK cleaning cuts");

        print_canvas.Clear();
        gr_cuts_500_1000.Draw();
        gPad->SetLogx();
        gPad->SetLogy();
        label_title = std::string("STK cleaning cuts - 500 GeV - 1 TeV cleaning cut");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf","Title:STK cleaning cuts");

        print_canvas.Clear();
        gr_cuts_1000_3000.Draw();
        gPad->SetLogx();
        gPad->SetLogy();
        label_title = std::string("STK cleaning cuts - 1 TeV - 3 TeV cleaning cut");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf","Title:STK cleaning cuts");

        /*
        print_canvas.Clear();
        gr_cuts_3000.Draw();
        gPad->SetLogx();
        gPad->SetLogy();
        label_title = std::string("STK cleaning cuts - > 3 TeV cleaning cut");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf)","Title:STK cleaning cuts");
        */

        print_canvas.Clear();

        f_cuts_20_100.SetTitle("20 GeV - 100 GeV; nStkClu1Rm; StkEcore1m");
        f_cuts_20_100.SetLineWidth(2);
        f_cuts_20_100.SetLineColor(kBlue-7);

        f_cuts_100_250.SetTitle("100 GeV - 250 GeV; nStkClu1Rm; StkEcore1m");
        f_cuts_100_250.SetLineWidth(2);
        f_cuts_100_250.SetLineColor(kBlue-9);

        f_cuts_250_500.SetTitle("250 GeV - 500 GeV; nStkClu1Rm; StkEcore1m");
        f_cuts_250_500.SetLineWidth(2);
        f_cuts_250_500.SetLineColor(kBlue-10);

        f_cuts_500_1000.SetTitle("500 GeV - 1 TeV; nStkClu1Rm; StkEcore1m");
        f_cuts_500_1000.SetLineWidth(2);
        f_cuts_500_1000.SetLineColor(kRed-10);

        f_cuts_1000_3000.SetTitle("1 TeV - 3 TeV; nStkClu1Rm; StkEcore1m");
        f_cuts_1000_3000.SetLineWidth(2);
        f_cuts_1000_3000.SetLineColor(kRed-9);

        print_canvas.Clear();

        f_cuts_20_100.Draw();
        f_cuts_100_250.Draw("same");
        f_cuts_250_500.Draw("same");
        f_cuts_500_1000.Draw("same");
        f_cuts_1000_3000.Draw("same");

        gPad->SetLogx();
        gPad->SetLogy();

        auto legend = print_canvas.BuildLegend();
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (auto primitiveObj :  *(legend->GetListOfPrimitives()))
        {
            auto primitive = (TLegendEntry*)primitiveObj;
            primitive->SetOption("l");
        }

        //f_cuts_3000.Draw("same");
        label_title = std::string("STK cleaning cuts - Cuts parametrizations");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();
        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf)","Title:STK cleaning cuts");

        /*
        print_canvas.Clear();

        TF1 f_cuts_energy ("f_cuts_energy", "[0] + [1]*x + [2]*pow(x, 2) + [3]*pow(x, 3)", 1, 1e+4);
        gr_energy_parametrization.Fit(&f_cuts_energy, "Q");
        gr_energy_parametrization.Draw();

        gr_energy_parametrization.GetXaxis()->SetRangeUser(10, 1e+4);
        gr_energy_parametrization.GetYaxis()->SetRangeUser(1e+2, 1e+6);

        label_title = std::string("STK cleaning cuts - energy parametrization");
        label = TPaveLabel(0.0, 0.97, 0.3, 1, label_title.c_str(), "tlNDC");
        label.Draw();

        print_canvas.Print("stk_cleaning_cuts_parametrization.pdf)","Title:STK cleaning cuts");
        */
       
        // Saving output functions to file
        TFile *func_out_file = TFile::Open("stk_cleaning_cuts_parametrization.root", "RECREATE");
        if (func_out_file->IsZombie())
        {
            std::cout << "Error writing output ROOT file [stk_cleaning_cuts_parametrization.root]" << std::endl;
            exit(100);
        }
        f_cuts_20_100.Write();
        f_cuts_100_250.Write();
        f_cuts_250_500.Write();
        f_cuts_500_1000.Write();
        f_cuts_1000_3000.Write();
        f_cuts_3000.Write();

        func_out_file->Close();

        print_canvas.Clear();


        print_canvas.Clear();
        print_canvas.Divide(2, 3);

        print_canvas.SetTicks();
        
        label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts", "tlNDC");
        label.Draw();

        print_canvas.cd(1);

        data_h_stk_cleaning_20_100->Draw();
        mc_h_stk_cleaning_20_100->Draw("same, colz");
        f_cuts_20_100.SetLineStyle(2);
        f_cuts_20_100.SetLineColor(kMagenta);
        f_cuts_20_100.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(2);

        data_h_stk_cleaning_100_250->Draw();
        mc_h_stk_cleaning_100_250->Draw("same, colz");
        f_cuts_100_250.SetLineStyle(2);
        f_cuts_100_250.SetLineColor(kMagenta);
        f_cuts_100_250.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(3);

        data_h_stk_cleaning_250_500->Draw();
        mc_h_stk_cleaning_250_500->Draw("same, colz");
        f_cuts_250_500.SetLineStyle(2);
        f_cuts_250_500.SetLineColor(kMagenta);
        f_cuts_250_500.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(4);

        data_h_stk_cleaning_500_1000->Draw();
        mc_h_stk_cleaning_500_1000->Draw("same, colz");
        f_cuts_500_1000.SetLineStyle(2);
        f_cuts_500_1000.SetLineColor(kMagenta);
        f_cuts_500_1000.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(5);

        data_h_stk_cleaning_1000_3000->Draw();
        mc_h_stk_cleaning_1000_3000->Draw("same, colz");
        f_cuts_1000_3000.SetLineStyle(2);
        f_cuts_1000_3000.SetLineColor(kMagenta);
        f_cuts_1000_3000.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(6);

        data_h_stk_cleaning_3000->Draw();
        mc_h_stk_cleaning_3000->Draw("same, colz");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.Print("stk_cleaning_cuts.pdf(","Title:STK cleaning cuts");

        label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts (xtrl < 12)", "tlNDC");
        label.Draw();

        print_canvas.cd(1);

        data_h_stk_cleaning_20_100_xtrl_12->Draw();
        mc_h_stk_cleaning_20_100_xtrl_12->Draw("same, colz");
        f_cuts_20_100.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(2);

        data_h_stk_cleaning_100_250_xtrl_12->Draw();
        mc_h_stk_cleaning_100_250_xtrl_12->Draw("same, colz");
        f_cuts_100_250.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(3);
        
        data_h_stk_cleaning_250_500_xtrl_12->Draw();
        mc_h_stk_cleaning_250_500_xtrl_12->Draw("same, colz");
        f_cuts_250_500.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(4);

        data_h_stk_cleaning_500_1000_xtrl_12->Draw();
        mc_h_stk_cleaning_500_1000_xtrl_12->Draw("same, colz");
        f_cuts_500_1000.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(5);

        data_h_stk_cleaning_1000_3000_xtrl_12->Draw();
        mc_h_stk_cleaning_1000_3000_xtrl_12->Draw("same, colz");
        f_cuts_1000_3000.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(6);

        data_h_stk_cleaning_3000_xtrl_12->Draw();
        mc_h_stk_cleaning_3000_xtrl_12->Draw("same, colz");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.Print("stk_cleaning_cuts.pdf","Title:STK cleaning cuts (xtrl < 12)");

        label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts (12 < xtrl < 100)", "tlNDC");
        label.Draw();

        print_canvas.cd(1);

        data_h_stk_cleaning_20_100_xtrl_12_100->Draw();
        mc_h_stk_cleaning_20_100_xtrl_12_100->Draw("same, colz");
        f_cuts_20_100.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(2);

        data_h_stk_cleaning_100_250_xtrl_12_100->Draw();
        mc_h_stk_cleaning_100_250_xtrl_12_100->Draw("same, colz");
        f_cuts_100_250.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(3);

        data_h_stk_cleaning_250_500_xtrl_12_100->Draw();
        mc_h_stk_cleaning_250_500_xtrl_12_100->Draw("same, colz");
        f_cuts_250_500.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(4);

        data_h_stk_cleaning_500_1000_xtrl_12_100->Draw();
        mc_h_stk_cleaning_500_1000_xtrl_12_100->Draw("same, colz");
        f_cuts_500_1000.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(5);

        data_h_stk_cleaning_1000_3000_xtrl_12_100->Draw();
        mc_h_stk_cleaning_1000_3000_xtrl_12_100->Draw("same, colz");
        f_cuts_1000_3000.Draw("same");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.cd(6);

        data_h_stk_cleaning_3000_xtrl_12_100->Draw();
        mc_h_stk_cleaning_3000_xtrl_12_100->Draw("same, colz");
        gPad->SetGrid(1,1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gStyle->SetOptStat(0);

        print_canvas.Print("stk_cleaning_cuts.pdf)","Title:STK cleaning cuts (12 < xtrl < 100)");

    }


void produceStkCleaningPlotsProjections_highlightSignal(const char* data_file, const char* mc_file) {

    TFile* datafile = TFile::Open(data_file, "READ");
    if (datafile->IsZombie()) {
        std::cerr << "\n\nError opening input DATA ROOT file [" << data_file << "]\n\n";
        exit(100);
    }

    auto data_h_stk_cleaning_20_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_20_100"));
    auto data_h_stk_cleaning_100_250 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_100_250"));
    auto data_h_stk_cleaning_250_500 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_250_500"));
    auto data_h_stk_cleaning_500_1000 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000"));
    auto data_h_stk_cleaning_1000_3000 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000"));
    auto data_h_stk_cleaning_3000 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_3000"));

    auto data_h_stk_cleaning_20_100_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12"));
    auto data_h_stk_cleaning_100_250_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12"));
    auto data_h_stk_cleaning_250_500_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12"));
    auto data_h_stk_cleaning_500_1000_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12"));
    auto data_h_stk_cleaning_1000_3000_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12"));
    auto data_h_stk_cleaning_3000_xtrl_12 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12"));

    auto data_h_stk_cleaning_20_100_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12_100"));
    auto data_h_stk_cleaning_100_250_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12_100"));
    auto data_h_stk_cleaning_250_500_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12_100"));
    auto data_h_stk_cleaning_500_1000_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12_100"));
    auto data_h_stk_cleaning_1000_3000_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12_100"));
    auto data_h_stk_cleaning_3000_xtrl_12_100 = static_cast<TH2D*>(datafile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12_100"));

    data_h_stk_cleaning_20_100->SetDirectory(0);
    data_h_stk_cleaning_100_250->SetDirectory(0);
    data_h_stk_cleaning_250_500->SetDirectory(0);
    data_h_stk_cleaning_500_1000->SetDirectory(0);
    data_h_stk_cleaning_1000_3000->SetDirectory(0);
    data_h_stk_cleaning_3000->SetDirectory(0);

    data_h_stk_cleaning_20_100_xtrl_12->SetDirectory(0);
    data_h_stk_cleaning_100_250_xtrl_12->SetDirectory(0);
    data_h_stk_cleaning_250_500_xtrl_12->SetDirectory(0);
    data_h_stk_cleaning_500_1000_xtrl_12->SetDirectory(0);
    data_h_stk_cleaning_1000_3000_xtrl_12->SetDirectory(0);
    data_h_stk_cleaning_3000_xtrl_12->SetDirectory(0);

    data_h_stk_cleaning_20_100_xtrl_12_100->SetDirectory(0);
    data_h_stk_cleaning_100_250_xtrl_12_100->SetDirectory(0);
    data_h_stk_cleaning_250_500_xtrl_12_100->SetDirectory(0);
    data_h_stk_cleaning_500_1000_xtrl_12_100->SetDirectory(0);
    data_h_stk_cleaning_1000_3000_xtrl_12_100->SetDirectory(0);
    data_h_stk_cleaning_3000_xtrl_12_100->SetDirectory(0);

    datafile->Close();

    TFile* mcfile = TFile::Open(mc_file, "READ");
    if (mcfile->IsZombie()) {
        std::cerr << "\n\nError opening input MC ROOT file [" << mc_file << "]\n\n";
        exit(100);
    }

    auto mc_h_stk_cleaning_20_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_20_100"));
    auto mc_h_stk_cleaning_100_250 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_100_250"));
    auto mc_h_stk_cleaning_250_500 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_250_500"));
    auto mc_h_stk_cleaning_500_1000 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000"));
    auto mc_h_stk_cleaning_1000_3000 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000"));
    auto mc_h_stk_cleaning_3000 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_3000"));

    auto mc_h_stk_cleaning_20_100_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12"));
    auto mc_h_stk_cleaning_100_250_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12"));
    auto mc_h_stk_cleaning_250_500_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12"));
    auto mc_h_stk_cleaning_500_1000_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12"));
    auto mc_h_stk_cleaning_1000_3000_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12"));
    auto mc_h_stk_cleaning_3000_xtrl_12 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12"));

    auto mc_h_stk_cleaning_20_100_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_20_100_xtrl_12_100"));
    auto mc_h_stk_cleaning_100_250_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_100_250_xtrl_12_100"));
    auto mc_h_stk_cleaning_250_500_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_250_500_xtrl_12_100"));
    auto mc_h_stk_cleaning_500_1000_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_500_1000_xtrl_12_100"));
    auto mc_h_stk_cleaning_1000_3000_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_1000_3000_xtrl_12_100"));
    auto mc_h_stk_cleaning_3000_xtrl_12_100 = static_cast<TH2D*>(mcfile->Get("stk_cleaning_cut/h_stk_cleaning_3000_xtrl_12_100"));

    mc_h_stk_cleaning_20_100->SetDirectory(0);
    mc_h_stk_cleaning_100_250->SetDirectory(0);
    mc_h_stk_cleaning_250_500->SetDirectory(0);
    mc_h_stk_cleaning_500_1000->SetDirectory(0);
    mc_h_stk_cleaning_1000_3000->SetDirectory(0);
    mc_h_stk_cleaning_3000->SetDirectory(0);

    mc_h_stk_cleaning_20_100_xtrl_12->SetDirectory(0);
    mc_h_stk_cleaning_100_250_xtrl_12->SetDirectory(0);
    mc_h_stk_cleaning_250_500_xtrl_12->SetDirectory(0);
    mc_h_stk_cleaning_500_1000_xtrl_12->SetDirectory(0);
    mc_h_stk_cleaning_1000_3000_xtrl_12->SetDirectory(0);
    mc_h_stk_cleaning_3000_xtrl_12->SetDirectory(0);

    mc_h_stk_cleaning_20_100_xtrl_12_100->SetDirectory(0);
    mc_h_stk_cleaning_100_250_xtrl_12_100->SetDirectory(0);
    mc_h_stk_cleaning_250_500_xtrl_12_100->SetDirectory(0);
    mc_h_stk_cleaning_500_1000_xtrl_12_100->SetDirectory(0);
    mc_h_stk_cleaning_1000_3000_xtrl_12_100->SetDirectory(0);
    mc_h_stk_cleaning_3000_xtrl_12_100->SetDirectory(0);

    mcfile->Close();

    TCanvas print_canvas("print_canvas", "print_canvas");
    print_canvas.Divide(2, 3);

    print_canvas.SetTicks();
    
    TPaveLabel label(0.0, 0.97, 0.3, 1, "STK cleaning cuts", "tlNDC");
    label.Draw();

    print_canvas.cd(1);

    data_h_stk_cleaning_20_100->Draw();
    mc_h_stk_cleaning_20_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(2);

    data_h_stk_cleaning_100_250->Draw();
    mc_h_stk_cleaning_100_250->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(3);

    data_h_stk_cleaning_250_500->Draw();
    mc_h_stk_cleaning_250_500->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(4);

    data_h_stk_cleaning_500_1000->Draw();
    mc_h_stk_cleaning_500_1000->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(5);

    data_h_stk_cleaning_1000_3000->Draw();
    mc_h_stk_cleaning_1000_3000->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(6);

    data_h_stk_cleaning_3000->Draw();
    mc_h_stk_cleaning_3000->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.Print("stk_cleaning_cuts.pdf(","Title:STK cleaning cuts");

    label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts (xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);

    data_h_stk_cleaning_20_100_xtrl_12->Draw();
    mc_h_stk_cleaning_20_100_xtrl_12->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(2);

    data_h_stk_cleaning_100_250_xtrl_12->Draw();
    mc_h_stk_cleaning_100_250_xtrl_12->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(3);
    
    data_h_stk_cleaning_250_500_xtrl_12->Draw();
    mc_h_stk_cleaning_250_500_xtrl_12->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(4);

    data_h_stk_cleaning_500_1000_xtrl_12->Draw();
    mc_h_stk_cleaning_500_1000_xtrl_12->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(5);

    data_h_stk_cleaning_1000_3000_xtrl_12->Draw();
    mc_h_stk_cleaning_1000_3000_xtrl_12->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(6);

    data_h_stk_cleaning_3000_xtrl_12->Draw();
    mc_h_stk_cleaning_3000_xtrl_12->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.Print("stk_cleaning_cuts.pdf","Title:STK cleaning cuts (xtrl < 12)");

    label = TPaveLabel(0.0, 0.97, 0.3, 1, "STK cleaning cuts (12 < xtrl < 100)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);

    data_h_stk_cleaning_20_100_xtrl_12_100->Draw();
    mc_h_stk_cleaning_20_100_xtrl_12_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(2);

    data_h_stk_cleaning_100_250_xtrl_12_100->Draw();
    mc_h_stk_cleaning_100_250_xtrl_12_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(3);

    data_h_stk_cleaning_250_500_xtrl_12_100->Draw();
    mc_h_stk_cleaning_250_500_xtrl_12_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(4);

    data_h_stk_cleaning_500_1000_xtrl_12_100->Draw();
    mc_h_stk_cleaning_500_1000_xtrl_12_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(5);

    data_h_stk_cleaning_1000_3000_xtrl_12_100->Draw();
    mc_h_stk_cleaning_1000_3000_xtrl_12_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.cd(6);

    data_h_stk_cleaning_3000_xtrl_12_100->Draw();
    mc_h_stk_cleaning_3000_xtrl_12_100->Draw("same, colz");
    gPad->SetGrid(1,1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    print_canvas.Print("stk_cleaning_cuts.pdf)","Title:STK cleaning cuts (12 < xtrl < 100)");

}