#include <iostream>

#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"

void psd_study_plots(const char* input_file_name, const char* output_file_name) {

    TFile *input_file = TFile::Open(input_file_name, "READ");
    if (input_file->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
        exit(100);
    }

    auto h_psd_stk_match_distance_x_20_100 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_20_100"));
    auto h_psd_stk_match_distance_x_100_250 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_100_250"));
    auto h_psd_stk_match_distance_x_250_500 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_250_500"));
    auto h_psd_stk_match_distance_x_500_1000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_500_1000"));
    auto h_psd_stk_match_distance_x_1000_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_1000_3000"));
    auto h_psd_stk_match_distance_x_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_3000"));
    auto h_psd_stk_match_distance_y_20_100 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_20_100"));
    auto h_psd_stk_match_distance_y_100_250 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_100_250"));
    auto h_psd_stk_match_distance_y_250_500 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_250_500"));
    auto h_psd_stk_match_distance_y_500_1000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_500_1000"));
    auto h_psd_stk_match_distance_y_1000_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_1000_3000"));
    auto h_psd_stk_match_distance_y_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_3000"));

    auto h_psd_stk_match_distance_x_within_psd_fvolume_20_100 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_within_psd_fvolume_20_100"));
    auto h_psd_stk_match_distance_x_within_psd_fvolume_100_250 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_within_psd_fvolume_100_250"));
    auto h_psd_stk_match_distance_x_within_psd_fvolume_250_500 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_within_psd_fvolume_250_500"));
    auto h_psd_stk_match_distance_x_within_psd_fvolume_500_1000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_within_psd_fvolume_500_1000"));
    auto h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000"));
    auto h_psd_stk_match_distance_x_within_psd_fvolume_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_within_psd_fvolume_3000"));
    auto h_psd_stk_match_distance_y_within_psd_fvolume_20_100 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_within_psd_fvolume_20_100"));
    auto h_psd_stk_match_distance_y_within_psd_fvolume_100_250 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_within_psd_fvolume_100_250"));
    auto h_psd_stk_match_distance_y_within_psd_fvolume_250_500 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_within_psd_fvolume_250_500"));
    auto h_psd_stk_match_distance_y_within_psd_fvolume_500_1000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_within_psd_fvolume_500_1000"));
    auto h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000"));
    auto h_psd_stk_match_distance_y_within_psd_fvolume_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_within_psd_fvolume_3000"));

    auto h_psd_stk_match_distance_x_outside_psd_fvolume_20_100 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_outside_psd_fvolume_20_100"));
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_100_250 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_outside_psd_fvolume_100_250"));
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_250_500 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_outside_psd_fvolume_250_500"));
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000"));
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000"));
    auto h_psd_stk_match_distance_x_outside_psd_fvolume_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_x_outside_psd_fvolume_3000"));
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_20_100 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_outside_psd_fvolume_20_100"));
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_100_250 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_outside_psd_fvolume_100_250"));
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_250_500 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_outside_psd_fvolume_250_500"));
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000"));
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000"));
    auto h_psd_stk_match_distance_y_outside_psd_fvolume_3000 = static_cast<TH1D*>(input_file->Get("psd_matching_distance/h_psd_stk_match_distance_y_outside_psd_fvolume_3000"));

    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_20_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_100_250"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_250_500"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_500_1000"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_1000_3000"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_3000"));
    
    auto h_stk_charge_psd_fvolume_psd_cut_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_20_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_100_250"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_250_500"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_500_1000"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_1000_3000"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_3000"));

    auto h_stk_charge_no_psd_cut_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_no_psd_cut_20_100"));
    auto h_stk_charge_no_psd_cut_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_no_psd_cut_100_250"));
    auto h_stk_charge_no_psd_cut_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_no_psd_cut_250_500"));
    auto h_stk_charge_no_psd_cut_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_no_psd_cut_500_1000"));
    auto h_stk_charge_no_psd_cut_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_no_psd_cut_1000_3000"));
    auto h_stk_charge_no_psd_cut_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_no_psd_cut_3000"));
    
    auto h_stk_charge_psd_cut_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_cut_20_100"));
    auto h_stk_charge_psd_cut_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_cut_100_250"));
    auto h_stk_charge_psd_cut_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_cut_250_500"));
    auto h_stk_charge_psd_cut_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_cut_500_1000"));
    auto h_stk_charge_psd_cut_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_cut_1000_3000"));
    auto h_stk_charge_psd_cut_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_cut_3000"));


    h_psd_stk_match_distance_x_20_100->SetDirectory(0);
    h_psd_stk_match_distance_x_100_250->SetDirectory(0);
    h_psd_stk_match_distance_x_250_500->SetDirectory(0);
    h_psd_stk_match_distance_x_500_1000->SetDirectory(0);
    h_psd_stk_match_distance_x_1000_3000->SetDirectory(0);
    h_psd_stk_match_distance_x_3000->SetDirectory(0);
    h_psd_stk_match_distance_y_20_100->SetDirectory(0);
    h_psd_stk_match_distance_y_100_250->SetDirectory(0);
    h_psd_stk_match_distance_y_250_500->SetDirectory(0);
    h_psd_stk_match_distance_y_500_1000->SetDirectory(0);
    h_psd_stk_match_distance_y_1000_3000->SetDirectory(0);
    h_psd_stk_match_distance_y_3000->SetDirectory(0);

    h_psd_stk_match_distance_x_within_psd_fvolume_20_100->SetDirectory(0);
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetDirectory(0);
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetDirectory(0);
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetDirectory(0);
    h_psd_stk_match_distance_x_within_psd_fvolume_1000_3000->SetDirectory(0);
    h_psd_stk_match_distance_x_within_psd_fvolume_3000->SetDirectory(0);
    h_psd_stk_match_distance_y_within_psd_fvolume_20_100->SetDirectory(0);
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetDirectory(0);
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetDirectory(0);
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetDirectory(0);
    h_psd_stk_match_distance_y_within_psd_fvolume_1000_3000->SetDirectory(0);
    h_psd_stk_match_distance_y_within_psd_fvolume_3000->SetDirectory(0);

    h_psd_stk_match_distance_x_outside_psd_fvolume_20_100->SetDirectory(0);
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetDirectory(0);
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetDirectory(0);
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetDirectory(0);
    h_psd_stk_match_distance_x_outside_psd_fvolume_1000_3000->SetDirectory(0);
    h_psd_stk_match_distance_x_outside_psd_fvolume_3000->SetDirectory(0);
    h_psd_stk_match_distance_y_outside_psd_fvolume_20_100->SetDirectory(0);
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetDirectory(0);
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetDirectory(0);
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetDirectory(0);
    h_psd_stk_match_distance_y_outside_psd_fvolume_1000_3000->SetDirectory(0);
    h_psd_stk_match_distance_y_outside_psd_fvolume_3000->SetDirectory(0);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_3000->SetDirectory(0);

    h_stk_charge_psd_fvolume_psd_cut_20_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_100_250->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_250_500->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_500_1000->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_3000->SetDirectory(0);

    h_stk_charge_no_psd_cut_20_100->SetDirectory(0);
    h_stk_charge_no_psd_cut_100_250->SetDirectory(0);
    h_stk_charge_no_psd_cut_250_500->SetDirectory(0);
    h_stk_charge_no_psd_cut_500_1000->SetDirectory(0);
    h_stk_charge_no_psd_cut_1000_3000->SetDirectory(0);
    h_stk_charge_no_psd_cut_3000->SetDirectory(0);

    h_stk_charge_psd_cut_20_100->SetDirectory(0);
    h_stk_charge_psd_cut_100_250->SetDirectory(0);
    h_stk_charge_psd_cut_250_500->SetDirectory(0);
    h_stk_charge_psd_cut_500_1000->SetDirectory(0);
    h_stk_charge_psd_cut_1000_3000->SetDirectory(0);
    h_stk_charge_psd_cut_3000->SetDirectory(0);

    input_file->Close();

    TFile* output_file = TFile::Open(output_file_name, "RECREATE");
    if (output_file->IsZombie()) {
        std::cerr << "\n\nError writing output ROOT file [" << output_file_name << "]\n\n";
        exit(100);
    }

    // Build X PSD distance canvas
    TCanvas c_psd_x_distance("c_psd_x_distance", "c_psd_x_distance", 500, 500);
    
    h_psd_stk_match_distance_x_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_x_100_250->SetTitle("100-250 GeV");
    h_psd_stk_match_distance_x_100_250->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_x_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_x_250_500->SetTitle("250-500 GeV");
    h_psd_stk_match_distance_x_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_x_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_x_500_1000->SetTitle("0.5-1 TeV");
    h_psd_stk_match_distance_x_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_x_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_x_100_250->Draw();
    h_psd_stk_match_distance_x_250_500->Draw("same");
    h_psd_stk_match_distance_x_500_1000->Draw("same");
    
    c_psd_x_distance.SetLogy();
    c_psd_x_distance.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_x_distance_legend = c_psd_x_distance.BuildLegend();
    c_psd_x_distance_legend->SetBorderSize(0);
    c_psd_x_distance_legend->SetFillStyle(0);
    auto c_psd_x_distance_primitives = c_psd_x_distance_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_x_distance_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_x_distance.Write();

    // Build X PSD distance canvas - within PSD fiducial volume
    TCanvas c_psd_x_distance_psd_fvolume ("c_psd_x_distance_psd_fvolume", "c_psd_x_distance_psd_fvolume", 500, 500);
    
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetTitle("100-250 GeV");
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetTitle("250-500 GeV");
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetTitle("0.5-1 TeV");
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->Draw();
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->Draw("same");
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->Draw("same");
    
    c_psd_x_distance_psd_fvolume.SetLogy();
    c_psd_x_distance_psd_fvolume.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_x_distance_psd_fvolume_legend = c_psd_x_distance_psd_fvolume.BuildLegend();
    c_psd_x_distance_psd_fvolume_legend->SetBorderSize(0);
    c_psd_x_distance_psd_fvolume_legend->SetFillStyle(0);
    auto c_psd_x_distance_psd_fvolume_primitives = c_psd_x_distance_psd_fvolume_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_x_distance_psd_fvolume_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_x_distance_psd_fvolume.Write();

    // Build X PSD distance canvas - outside PSD fiducial volume
    TCanvas c_psd_x_distance_outside_psd_fvolume ("c_psd_x_distance_outside_psd_fvolume", "c_psd_x_distance_outside_psd_fvolume", 500, 500);
    
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetTitle("100-250 GeV");
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetTitle("250-500 GeV");
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetTitle("0.5-1 TeV");
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->Draw();
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->Draw("same");
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->Draw("same");
    
    c_psd_x_distance_outside_psd_fvolume.SetLogy();
    c_psd_x_distance_outside_psd_fvolume.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_x_distance_outside_psd_fvolume_legend = c_psd_x_distance_outside_psd_fvolume.BuildLegend();
    c_psd_x_distance_outside_psd_fvolume_legend->SetBorderSize(0);
    c_psd_x_distance_outside_psd_fvolume_legend->SetFillStyle(0);
    auto c_psd_x_distance_outside_psd_fvolume_primitives = c_psd_x_distance_outside_psd_fvolume_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_x_distance_outside_psd_fvolume_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_x_distance_outside_psd_fvolume.Write();

    // Build X PSD distance canvas - 100 - 250 GeV focus
    TCanvas c_psd_x_distance_focus_100_250 ("c_psd_x_distance_focus_100_250", "c_psd_x_distance_focus_100_250", 500, 500);
    
    h_psd_stk_match_distance_x_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_x_100_250->SetTitle("100-250 GeV - total");
    h_psd_stk_match_distance_x_100_250->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetTitle("100-250 GeV - within PSD fiducial volume");
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetTitle("100-250 GeV - outside PSD fiducial volume");
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->SetLineWidth(2);
    
    h_psd_stk_match_distance_x_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_x_100_250->Draw();
    h_psd_stk_match_distance_x_within_psd_fvolume_100_250->Draw("same");
    h_psd_stk_match_distance_x_outside_psd_fvolume_100_250->Draw("same");
    
    c_psd_x_distance_focus_100_250.SetLogy();
    c_psd_x_distance_focus_100_250.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_x_distance_focus_100_250_legend = c_psd_x_distance_focus_100_250.BuildLegend();
    c_psd_x_distance_focus_100_250_legend->SetBorderSize(0);
    c_psd_x_distance_focus_100_250_legend->SetFillStyle(0);
    auto c_psd_x_distance_focus_100_250_primitives = c_psd_x_distance_focus_100_250_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_x_distance_focus_100_250_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_x_distance_focus_100_250.Write();

    // Build X PSD distance canvas - 250 - 500 GeV focus
    TCanvas c_psd_x_distance_focus_250_500 ("c_psd_x_distance_focus_250_500", "c_psd_x_distance_focus_250_500", 500, 500);
    
    h_psd_stk_match_distance_x_250_500->SetLineColor(kRed+2);
    h_psd_stk_match_distance_x_250_500->SetTitle("250-500 GeV - total");
    h_psd_stk_match_distance_x_250_500->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_250_500->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetTitle("250-500 GeV - within PSD fiducial volume");
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetTitle("250-500 GeV - outside PSD fiducial volume");
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->SetLineWidth(2);
    
    h_psd_stk_match_distance_x_250_500->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_x_250_500->Draw();
    h_psd_stk_match_distance_x_within_psd_fvolume_250_500->Draw("same");
    h_psd_stk_match_distance_x_outside_psd_fvolume_250_500->Draw("same");
    
    c_psd_x_distance_focus_250_500.SetLogy();
    c_psd_x_distance_focus_250_500.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_x_distance_focus_250_500_legend = c_psd_x_distance_focus_250_500.BuildLegend();
    c_psd_x_distance_focus_250_500_legend->SetBorderSize(0);
    c_psd_x_distance_focus_250_500_legend->SetFillStyle(0);
    auto c_psd_x_distance_focus_250_500_primitives = c_psd_x_distance_focus_250_500_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_x_distance_focus_250_500_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_x_distance_focus_250_500.Write();

    // Build X PSD distance canvas - 500 - 1000 GeV focus
    TCanvas c_psd_x_distance_focus_500_1000 ("c_psd_x_distance_focus_500_1000", "c_psd_x_distance_focus_500_1000", 500, 500);
    
    h_psd_stk_match_distance_x_500_1000->SetLineColor(kRed+2);
    h_psd_stk_match_distance_x_500_1000->SetTitle("250-500 GeV - total");
    h_psd_stk_match_distance_x_500_1000->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_500_1000->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_500_1000->SetLineWidth(2);

    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetTitle("250-500 GeV - within PSD fiducial volume");
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->GetXaxis()->SetTitle("#DeltaX (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->SetLineWidth(2);

    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetTitle("250-500 GeV - outside PSD fiducial volume");
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_x_500_1000->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_x_500_1000->Draw();
    h_psd_stk_match_distance_x_within_psd_fvolume_500_1000->Draw("same");
    h_psd_stk_match_distance_x_outside_psd_fvolume_500_1000->Draw("same");
    
    c_psd_x_distance_focus_500_1000.SetLogy();
    c_psd_x_distance_focus_500_1000.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_x_distance_focus_500_1000_legend = c_psd_x_distance_focus_500_1000.BuildLegend();
    c_psd_x_distance_focus_500_1000_legend->SetBorderSize(0);
    c_psd_x_distance_focus_500_1000_legend->SetFillStyle(0);
    auto c_psd_x_distance_focus_500_1000_primitives = c_psd_x_distance_focus_500_1000_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_x_distance_focus_500_1000_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_x_distance_focus_500_1000.Write();

    // Build Y PSD distance canvas
    TCanvas c_psd_y_distance("c_psd_y_distance", "c_psd_y_distance", 500, 500);
    
    h_psd_stk_match_distance_y_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_y_100_250->SetTitle("100-250 GeV");
    h_psd_stk_match_distance_y_100_250->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_y_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_y_250_500->SetTitle("250-500 GeV");
    h_psd_stk_match_distance_y_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_y_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_y_500_1000->SetTitle("0.5-1 TeV");
    h_psd_stk_match_distance_y_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_y_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_y_100_250->Draw();
    h_psd_stk_match_distance_y_250_500->Draw("same");
    h_psd_stk_match_distance_y_500_1000->Draw("same");
    
    c_psd_y_distance.SetLogy();
    c_psd_y_distance.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_y_distance_legend = c_psd_y_distance.BuildLegend();
    c_psd_y_distance_legend->SetBorderSize(0);
    c_psd_y_distance_legend->SetFillStyle(0);
    auto c_psd_y_distance_primitives = c_psd_y_distance_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_y_distance_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_y_distance.Write();

    // Build Y PSD distance canvas - within PSD fiducial volume
    TCanvas c_psd_y_distance_psd_fvolume ("c_psd_y_distance_psd_fvolume", "c_psd_y_distance_psd_fvolume", 500, 500);
    
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetTitle("100-250 GeV");
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetTitle("250-500 GeV");
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetTitle("0.5-1 TeV");
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->Draw();
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->Draw("same");
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->Draw("same");
    
    c_psd_y_distance_psd_fvolume.SetLogy();
    c_psd_y_distance_psd_fvolume.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_y_distance_psd_fvolume_legend = c_psd_y_distance_psd_fvolume.BuildLegend();
    c_psd_y_distance_psd_fvolume_legend->SetBorderSize(0);
    c_psd_y_distance_psd_fvolume_legend->SetFillStyle(0);
    auto c_psd_y_distance_psd_fvolume_primitives = c_psd_y_distance_psd_fvolume_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_y_distance_psd_fvolume_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_y_distance_psd_fvolume.Write();

    // Build Y PSD distance canvas - outside PSD fiducial volume
    TCanvas c_psd_y_distance_outside_psd_fvolume ("c_psd_y_distance_outside_psd_fvolume", "c_psd_y_distance_outside_psd_fvolume", 500, 500);
    
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetTitle("100-250 GeV");
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetTitle("250-500 GeV");
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetTitle("0.5-1 TeV");
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->Draw();
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->Draw("same");
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->Draw("same");
    
    c_psd_y_distance_outside_psd_fvolume.SetLogy();
    c_psd_y_distance_outside_psd_fvolume.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_y_distance_outside_psd_fvolume_legend = c_psd_y_distance_outside_psd_fvolume.BuildLegend();
    c_psd_y_distance_outside_psd_fvolume_legend->SetBorderSize(0);
    c_psd_y_distance_outside_psd_fvolume_legend->SetFillStyle(0);
    auto c_psd_y_distance_outside_psd_fvolume_primitives = c_psd_y_distance_outside_psd_fvolume_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_y_distance_outside_psd_fvolume_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_y_distance_outside_psd_fvolume.Write();

    // Build Y PSD distance canvas - 100 - 250 GeV focus
    TCanvas c_psd_y_distance_focus_100_250 ("c_psd_y_distance_focus_100_250", "c_psd_y_distance_focus_100_250", 500, 500);
    
    h_psd_stk_match_distance_y_100_250->SetLineColor(kRed+2);
    h_psd_stk_match_distance_y_100_250->SetTitle("100-250 GeV - total");
    h_psd_stk_match_distance_y_100_250->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetTitle("100-250 GeV - within PSD fiducial volume");
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->SetLineWidth(2);

    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetTitle("100-250 GeV - outside PSD fiducial volume");
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->SetLineWidth(2);
    
    h_psd_stk_match_distance_y_100_250->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_y_100_250->Draw();
    h_psd_stk_match_distance_y_within_psd_fvolume_100_250->Draw("same");
    h_psd_stk_match_distance_y_outside_psd_fvolume_100_250->Draw("same");
    
    c_psd_y_distance_focus_100_250.SetLogy();
    c_psd_y_distance_focus_100_250.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_y_distance_focus_100_250_legend = c_psd_y_distance_focus_100_250.BuildLegend();
    c_psd_y_distance_focus_100_250_legend->SetBorderSize(0);
    c_psd_y_distance_focus_100_250_legend->SetFillStyle(0);
    auto c_psd_y_distance_focus_100_250_primitives = c_psd_y_distance_focus_100_250_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_y_distance_focus_100_250_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_y_distance_focus_100_250.Write();

    // Build Y PSD distance canvas - 250 - 500 GeV focus
    TCanvas c_psd_y_distance_focus_250_500 ("c_psd_y_distance_focus_250_500", "c_psd_y_distance_focus_250_500", 500, 500);
    
    h_psd_stk_match_distance_y_250_500->SetLineColor(kRed+2);
    h_psd_stk_match_distance_y_250_500->SetTitle("250-500 GeV - total");
    h_psd_stk_match_distance_y_250_500->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_250_500->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetTitle("250-500 GeV - within PSD fiducial volume");
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->SetLineWidth(2);

    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetTitle("250-500 GeV - outside PSD fiducial volume");
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->SetLineWidth(2);
    
    h_psd_stk_match_distance_y_250_500->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_y_250_500->Draw();
    h_psd_stk_match_distance_y_within_psd_fvolume_250_500->Draw("same");
    h_psd_stk_match_distance_y_outside_psd_fvolume_250_500->Draw("same");
    
    c_psd_y_distance_focus_250_500.SetLogy();
    c_psd_y_distance_focus_250_500.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_y_distance_focus_250_500_legend = c_psd_y_distance_focus_250_500.BuildLegend();
    c_psd_y_distance_focus_250_500_legend->SetBorderSize(0);
    c_psd_y_distance_focus_250_500_legend->SetFillStyle(0);
    auto c_psd_y_distance_focus_250_500_primitives = c_psd_y_distance_focus_250_500_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_y_distance_focus_250_500_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_y_distance_focus_250_500.Write();
    
    // Build Y PSD distance canvas - 500 - 1000 GeV focus
    TCanvas c_psd_y_distance_focus_500_1000 ("c_psd_y_distance_focus_500_1000", "c_psd_y_distance_focus_500_1000", 500, 500);
    
    h_psd_stk_match_distance_y_500_1000->SetLineColor(kRed+2);
    h_psd_stk_match_distance_y_500_1000->SetTitle("250-500 GeV - total");
    h_psd_stk_match_distance_y_500_1000->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_500_1000->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_500_1000->SetLineWidth(2);

    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetLineColor(kGreen+2);
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetTitle("250-500 GeV - within PSD fiducial volume");
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->GetXaxis()->SetTitle("#DeltaY (Track-closest PSD hit) [mm]");
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->GetYaxis()->SetTitle("entries");
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->SetLineWidth(2);

    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetLineColor(kBlue+2);
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetTitle("250-500 GeV - outside PSD fiducial volume");
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->SetLineWidth(2);
    
    h_psd_stk_match_distance_y_500_1000->GetYaxis()->SetRangeUser(1, 1e+6);
    h_psd_stk_match_distance_y_500_1000->Draw();
    h_psd_stk_match_distance_y_within_psd_fvolume_500_1000->Draw("same");
    h_psd_stk_match_distance_y_outside_psd_fvolume_500_1000->Draw("same");
    
    c_psd_y_distance_focus_500_1000.SetLogy();
    c_psd_y_distance_focus_500_1000.SetTicks();
    gPad->SetGrid(1,1);

    auto c_psd_y_distance_focus_500_1000_legend = c_psd_y_distance_focus_500_1000.BuildLegend();
    c_psd_y_distance_focus_500_1000_legend->SetBorderSize(0);
    c_psd_y_distance_focus_500_1000_legend->SetFillStyle(0);
    auto c_psd_y_distance_focus_500_1000_primitives = c_psd_y_distance_focus_500_1000_legend->GetListOfPrimitives();
    for (auto primitiveObj :  *c_psd_y_distance_focus_500_1000_primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("l");
    }

    c_psd_y_distance_focus_500_1000.Write();

    // STK charge within PSD fiducial volume - no PSD charge cut
    TCanvas stk_charge_within_psd_fvolume_no_psd_cut ("stk_charge_within_psd_fvolume_no_psd_cut", "stk_charge_within_psd_fvolume_no_psd_cut");
    stk_charge_within_psd_fvolume_no_psd_cut.Divide(3, 1);

    stk_charge_within_psd_fvolume_no_psd_cut.cd(1);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250->SetTitle("100-250 GeV");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250->Draw("colz");
    stk_charge_within_psd_fvolume_no_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_within_psd_fvolume_no_psd_cut.cd(2);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500->SetTitle("250-500 GeV");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500->Draw("colz");
    stk_charge_within_psd_fvolume_no_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_within_psd_fvolume_no_psd_cut.cd(3);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000->SetTitle("500-1000 GeV");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000->Draw("colz");
    stk_charge_within_psd_fvolume_no_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_within_psd_fvolume_no_psd_cut.Write();

    // STK charge within PSD fiducial volume - PSD charge cut
    TCanvas stk_charge_within_psd_fvolume_psd_cut ("stk_charge_within_psd_fvolume_psd_cut", "stk_charge_within_psd_fvolume_psd_cut");
    stk_charge_within_psd_fvolume_psd_cut.Divide(3, 1);

    stk_charge_within_psd_fvolume_psd_cut.cd(1);
    h_stk_charge_psd_fvolume_psd_cut_100_250->SetTitle("100-250 GeV");
    h_stk_charge_psd_fvolume_psd_cut_100_250->Draw("colz");
    stk_charge_within_psd_fvolume_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_within_psd_fvolume_psd_cut.cd(2);
    h_stk_charge_psd_fvolume_psd_cut_250_500->SetTitle("250-500 GeV");
    h_stk_charge_psd_fvolume_psd_cut_250_500->Draw("colz");
    stk_charge_within_psd_fvolume_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_within_psd_fvolume_psd_cut.cd(3);
    h_stk_charge_psd_fvolume_psd_cut_500_1000->SetTitle("500-1000 GeV");
    h_stk_charge_psd_fvolume_psd_cut_500_1000->Draw("colz");
    stk_charge_within_psd_fvolume_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_within_psd_fvolume_psd_cut.Write();

    // STK charge outside PSD fiducial volume - no PSD charge cut
    TCanvas stk_charge_outside_psd_fvolume_no_psd_cut ("stk_charge_outside_psd_fvolume_no_psd_cut", "stk_charge_outside_psd_fvolume_no_psd_cut");
    stk_charge_outside_psd_fvolume_no_psd_cut.Divide(3, 1);

    stk_charge_outside_psd_fvolume_no_psd_cut.cd(1);
    h_stk_charge_no_psd_cut_100_250->SetTitle("100-250 GeV");
    h_stk_charge_no_psd_cut_100_250->Draw("colz");
    stk_charge_outside_psd_fvolume_no_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_outside_psd_fvolume_no_psd_cut.cd(2);
    h_stk_charge_no_psd_cut_250_500->SetTitle("250-500 GeV");
    h_stk_charge_no_psd_cut_250_500->Draw("colz");
    stk_charge_outside_psd_fvolume_no_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_outside_psd_fvolume_no_psd_cut.cd(3);
    h_stk_charge_no_psd_cut_500_1000->SetTitle("500-1000 GeV");
    h_stk_charge_no_psd_cut_500_1000->Draw("colz");
    stk_charge_outside_psd_fvolume_no_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_outside_psd_fvolume_no_psd_cut.Write();

    // STK charge within PSD fiducial volume - PSD charge cut
    TCanvas stk_charge_outside_psd_fvolume_psd_cut ("stk_charge_outside_psd_fvolume_psd_cut", "stk_charge_outside_psd_fvolume_psd_cut");
    stk_charge_outside_psd_fvolume_psd_cut.Divide(3, 1);

    stk_charge_outside_psd_fvolume_psd_cut.cd(1);
    h_stk_charge_psd_cut_100_250->SetTitle("100-250 GeV");
    h_stk_charge_psd_cut_100_250->Draw("colz");
    stk_charge_outside_psd_fvolume_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_outside_psd_fvolume_psd_cut.cd(2);
    h_stk_charge_psd_cut_250_500->SetTitle("250-500 GeV");
    h_stk_charge_psd_cut_250_500->Draw("colz");
    stk_charge_outside_psd_fvolume_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);

    stk_charge_outside_psd_fvolume_psd_cut.cd(3);
    h_stk_charge_psd_cut_500_1000->SetTitle("500-1000 GeV");
    h_stk_charge_psd_cut_500_1000->Draw("colz");
    stk_charge_outside_psd_fvolume_psd_cut.SetTicks();
    gPad->SetLogz();
    gPad->SetGrid(1,1);
    
    stk_charge_outside_psd_fvolume_psd_cut.Write();
    
    output_file->Close();
}