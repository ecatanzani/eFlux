#include <iostream>

#include "TH2D.h"
#include "TPDF.h"
#include "TPad.h"
#include "TLine.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"

void stk_charge(const char* input_file_name) {
    
    TFile *input_file = TFile::Open(input_file_name, "READ");
    if (input_file->IsZombie()) {
        std::cerr << "\n\nError opening input file [" << input_file << "]\n\n";
        exit(100);
    }

    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_20_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_100_250"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_250_500"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_500_1000"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_1000_3000"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_3000"));

    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12"));

    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100"));

    auto h_stk_charge_psd_fvolume_psd_cut_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_20_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_100_250"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_250_500"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_500_1000"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_1000_3000"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_3000"));
    
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12"));

    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100"));

    auto h_stk_charge_20_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_20_100"));
    auto h_stk_charge_100_250 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_100_250"));
    auto h_stk_charge_250_500 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_250_500"));
    auto h_stk_charge_500_1000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_500_1000"));
    auto h_stk_charge_1000_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_1000_3000"));
    auto h_stk_charge_3000 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_cut_3000"));

    auto h_stk_charge_20_100_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_20_100_xtrl_12"));
    auto h_stk_charge_100_250_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_100_250_xtrl_12"));
    auto h_stk_charge_250_500_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_250_500_xtrl_12"));
    auto h_stk_charge_500_1000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_500_1000_xtrl_12"));
    auto h_stk_charge_1000_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_1000_3000_xtrl_12"));
    auto h_stk_charge_3000_xtrl_12 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_3000_xtrl_12"));
    
    auto h_stk_charge_20_100_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_20_100_xtrl_12_100"));
    auto h_stk_charge_100_250_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_100_250_xtrl_12_100"));
    auto h_stk_charge_250_500_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_250_500_xtrl_12_100"));
    auto h_stk_charge_500_1000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_500_1000_xtrl_12_100"));
    auto h_stk_charge_1000_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_1000_3000_xtrl_12_100"));
    auto h_stk_charge_3000_xtrl_12_100 = static_cast<TH2D*>(input_file->Get("stk_charge/h_stk_charge_3000_xtrl_12_100"));

    h_stk_charge_psd_fvolume_no_psd_cut_20_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_3000->SetDirectory(0);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12->SetDirectory(0);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100->SetDirectory(0);

    h_stk_charge_psd_fvolume_psd_cut_20_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_100_250->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_250_500->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_500_1000->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_3000->SetDirectory(0);
    
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12->SetDirectory(0);

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100->SetDirectory(0);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100->SetDirectory(0);

    h_stk_charge_20_100->SetDirectory(0);
    h_stk_charge_100_250->SetDirectory(0);
    h_stk_charge_250_500->SetDirectory(0);
    h_stk_charge_500_1000->SetDirectory(0);
    h_stk_charge_1000_3000->SetDirectory(0);
    h_stk_charge_3000->SetDirectory(0);

    h_stk_charge_20_100_xtrl_12->SetDirectory(0);
    h_stk_charge_100_250_xtrl_12->SetDirectory(0);
    h_stk_charge_250_500_xtrl_12->SetDirectory(0);
    h_stk_charge_500_1000_xtrl_12->SetDirectory(0);
    h_stk_charge_1000_3000_xtrl_12->SetDirectory(0);
    h_stk_charge_3000_xtrl_12->SetDirectory(0);
    
    h_stk_charge_20_100_xtrl_12_100->SetDirectory(0);
    h_stk_charge_100_250_xtrl_12_100->SetDirectory(0);
    h_stk_charge_250_500_xtrl_12_100->SetDirectory(0);
    h_stk_charge_500_1000_xtrl_12_100->SetDirectory(0);
    h_stk_charge_1000_3000_xtrl_12_100->SetDirectory(0);
    h_stk_charge_3000_xtrl_12_100->SetDirectory(0);

    input_file->Close();

    TLine l1(0, 0, 0, 25);
    TLine l2(0, 25, 10, 25);
    TLine l3(10, 25, 10, 40);
    TLine l4(10, 40, 40, 40);
    TLine l5(40, 40, 40, 10);
    TLine l6(40, 10, 25, 10);
    TLine l7(25, 10, 25, 0);
    TLine l8(25, 0, 0, 0);

    l1.SetLineWidth(2);
    l1.SetLineStyle(2);
    l1.SetLineColor(kMagenta);

    l2.SetLineWidth(2);
    l2.SetLineStyle(2);
    l2.SetLineColor(kMagenta);

    l3.SetLineWidth(2);
    l3.SetLineStyle(2);
    l3.SetLineColor(kMagenta);

    l4.SetLineWidth(2);
    l4.SetLineStyle(2);
    l4.SetLineColor(kMagenta);

    l5.SetLineWidth(2);
    l5.SetLineStyle(2);
    l5.SetLineColor(kMagenta);

    l6.SetLineWidth(2);
    l6.SetLineStyle(2);
    l6.SetLineColor(kMagenta);

    l7.SetLineWidth(2);
    l7.SetLineStyle(2);
    l7.SetLineColor(kMagenta);

    l8.SetLineWidth(2);
    l8.SetLineStyle(2);
    l8.SetLineColor(kMagenta);


    TCanvas print_canvas("print_canvas", "print_canvas");
    print_canvas.Divide(2, 3);

    print_canvas.SetTicks();
    print_canvas.SetLogy();

    TPaveLabel label(0.0, 0.95, 0.3, 1, "STK charge within PSD fiducial volume (no PSD charge cut)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");

    gPad->SetGrid(1,1);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_no_psd_cut_3000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);

    print_canvas.Print("stk_charge.pdf(","Title:STK charges within PSD fiducial volume (no PSD charge cut)");

    label.Draw();
    
    print_canvas.cd(1);
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_20_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_20_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();

    TLegend legend(0.4,0.7,0.7,0.9);
    legend.AddEntry(h_stk_charge_psd_fvolume_no_psd_cut_20_100_chx,"STK Charge X","l");
    legend.AddEntry(h_stk_charge_psd_fvolume_no_psd_cut_20_100_chy,"STK Charge Y","l");
    legend.Draw();
    
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250->Draw("colz");

    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_100_250->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_100_250->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_100_250_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500->Draw("colz");

    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_250_500->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_250_500->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_250_500_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000->Draw("colz");

    
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_500_1000->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_500_1000->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000->Draw("colz");

    
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_1000_3000->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_1000_3000->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_no_psd_cut_3000->Draw("colz");

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_3000->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_3000_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_3000->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_3000_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_3000_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_3000_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_3000_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_3000_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_3000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_3000_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_3000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut)");


    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge within PSD fiducial volume (no PSD charge cut - xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut - xtrl < 12)");

    print_canvas.cd(1);
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);

    
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut - xtrl < 12)");


    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge within PSD fiducial volume (no PSD charge cut - xtrl < 12 < 100)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut - xtrl < 12 < 100)");

    print_canvas.cd(1);
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_20_100_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_100_250_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_250_500_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_500_1000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);

    
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_1000_3000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_no_psd_cut_3000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut - xtrl < 12 < 100)");
    
    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge within PSD fiducial volume (PSD charge cut)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_psd_fvolume_psd_cut_20_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_psd_cut_100_250->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_psd_cut_250_500->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_psd_cut_500_1000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_psd_cut_3000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (PSD charge cut)");

    print_canvas.cd(1);
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_20_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_20_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_20_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_20_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_20_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_20_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_20_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_20_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_20_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_20_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_20_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_20_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_20_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_psd_fvolume_psd_cut_100_250_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_100_250->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_100_250_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_100_250->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_100_250_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_100_250_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_100_250_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_100_250_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_100_250_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_100_250_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_100_250_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_100_250_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_100_250_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_100_250_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_psd_fvolume_psd_cut_250_500_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_250_500->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_250_500_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_250_500->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_250_500_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_250_500_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_250_500_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_250_500_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_250_500_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_250_500_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_250_500_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_250_500_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_250_500_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_250_500_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_500_1000->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_500_1000_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_500_1000->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_500_1000_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_500_1000_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_500_1000_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_500_1000_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_500_1000_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_500_1000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_500_1000_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_500_1000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);

    
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_1000_3000->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_1000_3000->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_1000_3000_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_psd_fvolume_psd_cut_3000_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_3000->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_3000_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_3000->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_3000_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_3000_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_3000_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_3000_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_3000_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_3000_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_3000_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_3000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_3000_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_3000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (PSD charge cut)");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge within PSD fiducial volume (PSD charge cut - xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (PSD charge cut - xtrl < 12)");

    print_canvas.cd(1);
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);
    
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (PSD charge cut - xtrl < 12)");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge within PSD fiducial volume (PSD charge cut - xtrl < 12 < 100)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (PSD charge cut - xtrl < 12 < 100)");

    print_canvas.cd(1);
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_20_100_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_100_250_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_250_500_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_500_1000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);
    
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_1000_3000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100->ProjectionX("h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx"));
    auto h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100->ProjectionY("h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chy"));
    
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chx->Draw();
    h_stk_charge_psd_fvolume_psd_cut_3000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (PSD charge cut - xtrl < 12 < 100)");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge outside PSD fiducial volume (no PSD/STK match)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_20_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_100_250->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_250_500->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_500_1000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_1000_3000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_3000->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges outside PSD fiducial volume (no PSD/STK match)");

    print_canvas.cd(1);
    auto h_stk_charge_20_100_chx = static_cast<TH1D*>(h_stk_charge_20_100->ProjectionX("h_stk_charge_20_100_chx"));
    auto h_stk_charge_20_100_chy = static_cast<TH1D*>(h_stk_charge_20_100->ProjectionY("h_stk_charge_20_100_chy"));
    
    h_stk_charge_20_100_chx->SetLineColor(kBlue);
    h_stk_charge_20_100_chy->SetLineColor(kRed);

    h_stk_charge_20_100_chx->SetLineWidth(2);
    h_stk_charge_20_100_chy->SetLineWidth(2);

    h_stk_charge_20_100_chx->SetTitle("Charge X");
    h_stk_charge_20_100_chy->SetTitle("Charge Y");

    h_stk_charge_20_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_20_100_chx->Draw();
    h_stk_charge_20_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_100_250_chx = static_cast<TH1D*>(h_stk_charge_100_250->ProjectionX("h_stk_charge_100_250_chx"));
    auto h_stk_charge_100_250_chy = static_cast<TH1D*>(h_stk_charge_100_250->ProjectionY("h_stk_charge_100_250_chy"));
    
    h_stk_charge_100_250_chx->SetLineColor(kBlue);
    h_stk_charge_100_250_chy->SetLineColor(kRed);

    h_stk_charge_100_250_chx->SetLineWidth(2);
    h_stk_charge_100_250_chy->SetLineWidth(2);

    h_stk_charge_100_250_chx->SetTitle("Charge X");
    h_stk_charge_100_250_chy->SetTitle("Charge Y");

    h_stk_charge_100_250_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_100_250_chx->Draw();
    h_stk_charge_100_250_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_250_500_chx = static_cast<TH1D*>(h_stk_charge_250_500->ProjectionX("h_stk_charge_250_500_chx"));
    auto h_stk_charge_250_500_chy = static_cast<TH1D*>(h_stk_charge_250_500->ProjectionY("h_stk_charge_250_500_chy"));
    
    h_stk_charge_250_500_chx->SetLineColor(kBlue);
    h_stk_charge_250_500_chy->SetLineColor(kRed);

    h_stk_charge_250_500_chx->SetLineWidth(2);
    h_stk_charge_250_500_chy->SetLineWidth(2);

    h_stk_charge_250_500_chx->SetTitle("Charge X");
    h_stk_charge_250_500_chy->SetTitle("Charge Y");

    h_stk_charge_250_500_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_250_500_chx->Draw();
    h_stk_charge_250_500_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_500_1000_chx = static_cast<TH1D*>(h_stk_charge_500_1000->ProjectionX("h_stk_charge_500_1000_chx"));
    auto h_stk_charge_500_1000_chy = static_cast<TH1D*>(h_stk_charge_500_1000->ProjectionY("h_stk_charge_500_1000_chy"));
    
    h_stk_charge_500_1000_chx->SetLineColor(kBlue);
    h_stk_charge_500_1000_chy->SetLineColor(kRed);

    h_stk_charge_500_1000_chx->SetLineWidth(2);
    h_stk_charge_500_1000_chy->SetLineWidth(2);

    h_stk_charge_500_1000_chx->SetTitle("Charge X");
    h_stk_charge_500_1000_chy->SetTitle("Charge Y");

    h_stk_charge_500_1000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_500_1000_chx->Draw();
    h_stk_charge_500_1000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);

    
    auto h_stk_charge_1000_3000_chx = static_cast<TH1D*>(h_stk_charge_1000_3000->ProjectionX("h_stk_charge_1000_3000_chx"));
    auto h_stk_charge_1000_3000_chy = static_cast<TH1D*>(h_stk_charge_1000_3000->ProjectionY("h_stk_charge_1000_3000_chy"));
    
    h_stk_charge_1000_3000_chx->SetLineColor(kBlue);
    h_stk_charge_1000_3000_chy->SetLineColor(kRed);

    h_stk_charge_1000_3000_chx->SetLineWidth(2);
    h_stk_charge_1000_3000_chy->SetLineWidth(2);

    h_stk_charge_1000_3000_chx->SetTitle("Charge X");
    h_stk_charge_1000_3000_chy->SetTitle("Charge Y");

    h_stk_charge_1000_3000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_1000_3000_chx->Draw();
    h_stk_charge_1000_3000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_3000_chx = static_cast<TH1D*>(h_stk_charge_3000->ProjectionX("h_stk_charge_3000_chx"));
    auto h_stk_charge_3000_chy = static_cast<TH1D*>(h_stk_charge_3000->ProjectionY("h_stk_charge_3000_chy"));
    
    h_stk_charge_3000_chx->SetLineColor(kBlue);
    h_stk_charge_3000_chy->SetLineColor(kRed);

    h_stk_charge_3000_chx->SetLineWidth(2);
    h_stk_charge_3000_chy->SetLineWidth(2);

    h_stk_charge_3000_chx->SetTitle("Charge X");
    h_stk_charge_3000_chy->SetTitle("Charge Y");

    h_stk_charge_3000_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_3000_chx->Draw();
    h_stk_charge_3000_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut)");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge outside PSD fiducial volume (no PSD/STK match - xtrl < 12)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_20_100_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_100_250_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_250_500_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_500_1000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_1000_3000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_3000_xtrl_12->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges outside PSD fiducial volume (no PSD/STK match - xtrl < 12)");

    print_canvas.cd(1);
    auto h_stk_charge_20_100_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_20_100_xtrl_12->ProjectionX("h_stk_charge_20_100_xtrl_12_chx"));
    auto h_stk_charge_20_100_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_20_100_xtrl_12->ProjectionY("h_stk_charge_20_100_xtrl_12_chy"));
    
    h_stk_charge_20_100_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_20_100_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_20_100_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_20_100_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_20_100_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_20_100_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_20_100_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_20_100_xtrl_12_chx->Draw();
    h_stk_charge_20_100_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_100_250_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_100_250_xtrl_12->ProjectionX("h_stk_charge_100_250_xtrl_12_chx"));
    auto h_stk_charge_100_250_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_100_250_xtrl_12->ProjectionY("h_stk_charge_100_250_xtrl_12_chy"));
    
    h_stk_charge_100_250_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_100_250_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_100_250_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_100_250_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_100_250_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_100_250_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_100_250_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_100_250_xtrl_12_chx->Draw();
    h_stk_charge_100_250_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_250_500_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_250_500_xtrl_12->ProjectionX("h_stk_charge_250_500_xtrl_12_chx"));
    auto h_stk_charge_250_500_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_250_500_xtrl_12->ProjectionY("h_stk_charge_250_500_xtrl_12_chy"));
    
    h_stk_charge_250_500_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_250_500_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_250_500_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_250_500_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_250_500_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_250_500_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_250_500_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_250_500_xtrl_12_chx->Draw();
    h_stk_charge_250_500_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_500_1000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_500_1000_xtrl_12->ProjectionX("h_stk_charge_500_1000_xtrl_12_chx"));
    auto h_stk_charge_500_1000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_500_1000_xtrl_12->ProjectionY("h_stk_charge_500_1000_xtrl_12_chx"));
    
    h_stk_charge_500_1000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_500_1000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_500_1000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_500_1000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_500_1000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_500_1000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_500_1000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_500_1000_xtrl_12_chx->Draw();
    h_stk_charge_500_1000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);

    
    auto h_stk_charge_1000_3000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_1000_3000_xtrl_12->ProjectionX("h_stk_charge_1000_3000_xtrl_12_chx"));
    auto h_stk_charge_1000_3000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_1000_3000_xtrl_12->ProjectionY("h_stk_charge_1000_3000_xtrl_12_chy"));
    
    h_stk_charge_1000_3000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_1000_3000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_1000_3000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_1000_3000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_1000_3000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_1000_3000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_1000_3000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_1000_3000_xtrl_12_chx->Draw();
    h_stk_charge_1000_3000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_3000_xtrl_12_chx = static_cast<TH1D*>(h_stk_charge_3000_xtrl_12->ProjectionX("h_stk_charge_3000_xtrl_12_chx"));
    auto h_stk_charge_3000_xtrl_12_chy = static_cast<TH1D*>(h_stk_charge_3000_xtrl_12->ProjectionY("h_stk_charge_3000_xtrl_12_chy"));
    
    h_stk_charge_3000_xtrl_12_chx->SetLineColor(kBlue);
    h_stk_charge_3000_xtrl_12_chy->SetLineColor(kRed);

    h_stk_charge_3000_xtrl_12_chx->SetLineWidth(2);
    h_stk_charge_3000_xtrl_12_chy->SetLineWidth(2);

    h_stk_charge_3000_xtrl_12_chx->SetTitle("Charge X");
    h_stk_charge_3000_xtrl_12_chy->SetTitle("Charge Y");

    h_stk_charge_3000_xtrl_12_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_3000_xtrl_12_chx->Draw();
    h_stk_charge_3000_xtrl_12_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf","Title:STK charges within PSD fiducial volume (no PSD charge cut - xtrl < 12 < 100)");

    label = TPaveLabel(0.0, 0.95, 0.3, 1, "STK charge outside PSD fiducial volume (no PSD/STK match - xtrl < 12 < 100)", "tlNDC");
    label.Draw();

    print_canvas.cd(1);
    h_stk_charge_20_100_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(2);
    h_stk_charge_100_250_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(3);
    h_stk_charge_250_500_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(4);
    h_stk_charge_500_1000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(5);
    h_stk_charge_1000_3000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();
    print_canvas.cd(6);
    h_stk_charge_3000_xtrl_12_100->Draw("colz");
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    l5.Draw("same");
    l6.Draw("same");
    l7.Draw("same");
    l8.Draw("same");
    gPad->SetGrid(1,1);
    gPad->SetLogy(0);
    gPad->SetLogz();

    print_canvas.Print("stk_charge.pdf","Title:STK charges outside PSD fiducial volume (no PSD/STK match - xtrl < 12 < 100)");

    print_canvas.cd(1);
    auto h_stk_charge_20_100_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_20_100_xtrl_12_100->ProjectionX("h_stk_charge_20_100_xtrl_12_100_chx"));
    auto h_stk_charge_20_100_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_20_100_xtrl_12_100->ProjectionY("h_stk_charge_20_100_xtrl_12_100_chy"));
    
    h_stk_charge_20_100_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_20_100_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_20_100_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_20_100_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_20_100_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_20_100_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_20_100_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_20_100_xtrl_12_100_chx->Draw();
    h_stk_charge_20_100_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();
    
    print_canvas.cd(2);

    auto h_stk_charge_100_250_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_100_250_xtrl_12_100->ProjectionX("h_stk_charge_100_250_xtrl_12_100_chx"));
    auto h_stk_charge_100_250_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_100_250_xtrl_12_100->ProjectionY("h_stk_charge_100_250_xtrl_12_100_chy"));
    
    h_stk_charge_100_250_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_100_250_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_100_250_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_100_250_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_100_250_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_100_250_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_100_250_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_100_250_xtrl_12_100_chx->Draw();
    h_stk_charge_100_250_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(3);

    auto h_stk_charge_250_500_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_250_500_xtrl_12_100->ProjectionX("h_stk_charge_250_500_xtrl_12_100_chx"));
    auto h_stk_charge_250_500_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_250_500_xtrl_12_100->ProjectionY("h_stk_charge_250_500_xtrl_12_100_chy"));
    
    h_stk_charge_250_500_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_250_500_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_250_500_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_250_500_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_250_500_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_250_500_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_250_500_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_250_500_xtrl_12_100_chx->Draw();
    h_stk_charge_250_500_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(4);
    
    auto h_stk_charge_500_1000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_500_1000_xtrl_12_100->ProjectionX("h_stk_charge_500_1000_xtrl_12_100_chx"));
    auto h_stk_charge_500_1000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_500_1000_xtrl_12_100->ProjectionY("h_stk_charge_500_1000_xtrl_12_100_chy"));
    
    h_stk_charge_500_1000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_500_1000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_500_1000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_500_1000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_500_1000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_500_1000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_500_1000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_500_1000_xtrl_12_100_chx->Draw();
    h_stk_charge_500_1000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.cd(5);

    
    auto h_stk_charge_1000_3000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_1000_3000_xtrl_12->ProjectionX("h_stk_charge_1000_3000_xtrl_12_100_chx"));
    auto h_stk_charge_1000_3000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_1000_3000_xtrl_12->ProjectionY("h_stk_charge_1000_3000_xtrl_12_100_chy"));
    
    h_stk_charge_1000_3000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_1000_3000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_1000_3000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_1000_3000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_1000_3000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_1000_3000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_1000_3000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_1000_3000_xtrl_12_100_chx->Draw();
    h_stk_charge_1000_3000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
   
    legend.Draw();

    print_canvas.cd(6);

    auto h_stk_charge_3000_xtrl_12_100_chx = static_cast<TH1D*>(h_stk_charge_3000_xtrl_12_100->ProjectionX("h_stk_charge_3000_xtrl_12_100_chx"));
    auto h_stk_charge_3000_xtrl_12_100_chy = static_cast<TH1D*>(h_stk_charge_3000_xtrl_12_100->ProjectionY("h_stk_charge_3000_xtrl_12_100_chy"));
    
    h_stk_charge_3000_xtrl_12_100_chx->SetLineColor(kBlue);
    h_stk_charge_3000_xtrl_12_100_chy->SetLineColor(kRed);

    h_stk_charge_3000_xtrl_12_100_chx->SetLineWidth(2);
    h_stk_charge_3000_xtrl_12_100_chy->SetLineWidth(2);

    h_stk_charge_3000_xtrl_12_100_chx->SetTitle("Charge X");
    h_stk_charge_3000_xtrl_12_100_chy->SetTitle("Charge Y");

    h_stk_charge_3000_xtrl_12_100_chx->GetXaxis()->SetTitle("STK charge");
    h_stk_charge_3000_xtrl_12_100_chx->Draw();
    h_stk_charge_3000_xtrl_12_100_chy->Draw("same");
    
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    
    legend.Draw();

    print_canvas.Print("stk_charge.pdf)","Title:STK charges within PSD fiducial volume (no PSD charge cut - xtrl < 12 < 100)");
}