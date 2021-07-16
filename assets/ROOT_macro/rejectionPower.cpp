#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>

void rejectionPower(const char* electronMC, const char* protonMC, const char* outfile)
{
    TFile electron_file(electronMC, "READ");
    if (!electron_file.IsOpen())
    {
        std::cerr << "\n\nError reading input file [" << electronMC << "]\n\n";
        exit(100);
    }
    auto h_electrons = static_cast<TH1D*>(electron_file.Get("h_incoming"));
    h_electrons->SetName("h_incoming_electrons");
    auto h_electrons_sel = static_cast<TH1D*>(electron_file.Get("h_all_cut"));
    h_electrons_sel->SetName("h_all_cut_electrons");
    h_electrons->SetDirectory(0);
    h_electrons_sel->SetDirectory(0);
    electron_file.Close();
    
    TFile proton_file(protonMC, "READ");
    if (!proton_file.IsOpen())
    {
        std::cerr << "\n\nError reading input file [" << protonMC << "]\n\n";
        exit(100);
    }
    auto h_protons = static_cast<TH1D*>(proton_file.Get("h_incoming"));
    h_protons->SetName("h_incoming_protons");
    auto h_protons_sel = static_cast<TH1D*>(proton_file.Get("h_all_cut"));
    h_protons_sel->SetName("h_all_cut_protons");
    h_protons->SetDirectory(0);
    h_protons_sel->SetDirectory(0);
    proton_file.Close();
    
    auto electron_eff = static_cast<TH1D*>(h_electrons_sel->Clone("electron_eff"));
    electron_eff->SetName("electron_eff");
    electron_eff->Divide(h_electrons);
    auto proton_eff = static_cast<TH1D*>(h_protons_sel->Clone("proton_eff"));
    proton_eff->SetName("proton_eff");
    proton_eff->Divide(h_protons);

    auto rpower = static_cast<TH1D*>(electron_eff->Clone("rpower"));
    rpower->SetName("rpower");
    rpower->Divide(proton_eff);

    TFile out_file(outfile, "RECREATE");
    if (!out_file.IsOpen())
    {
        std::cerr << "\n\nError reading input file [" << outfile << "]\n\n";
        exit(100);
    }

    h_electrons->Write();
    h_electrons_sel->Write();
    h_protons->Write();
    h_protons_sel->Write();
    electron_eff->Write();
    proton_eff->Write();
    rpower->Write();

    TCanvas power_c("power_c", "power_c",11,49,700,502);
	gStyle->SetOptFit(1);
  	power_c.Range(-1.008394,-80,3.707364,453.3333);
  	power_c.SetFillColor(0);
  	power_c.SetBorderMode(0);
  	power_c.SetBorderSize(2);
  	power_c.SetLogx();
    power_c.SetLogy();
	power_c.SetLeftMargin(0.15);
  	power_c.SetRightMargin(0.15);
  	power_c.SetBottomMargin(0.15);
  	power_c.SetFrameLineWidth(2);
  	power_c.SetFrameBorderMode(0);
  	power_c.SetFrameBorderSize(2);
  	power_c.SetFrameLineWidth(2);
  	power_c.SetFrameBorderMode(0);
  	power_c.SetFrameBorderSize(2);
    rpower->Draw();
    power_c.Write();

    TCanvas ele_eff("ele_eff", "ele_eff",11,49,700,502);
	gStyle->SetOptFit(1);
  	ele_eff.Range(-1.008394,-80,3.707364,453.3333);
  	ele_eff.SetFillColor(0);
  	ele_eff.SetBorderMode(0);
  	ele_eff.SetBorderSize(2);
  	ele_eff.SetLogx();
	ele_eff.SetLeftMargin(0.15);
  	ele_eff.SetRightMargin(0.15);
  	ele_eff.SetBottomMargin(0.15);
  	ele_eff.SetFrameLineWidth(2);
  	ele_eff.SetFrameBorderMode(0);
  	ele_eff.SetFrameBorderSize(2);
  	ele_eff.SetFrameLineWidth(2);
  	ele_eff.SetFrameBorderMode(0);
  	ele_eff.SetFrameBorderSize(2);
    electron_eff->Draw();
    ele_eff.Write();

    TCanvas proton_eff_c("proton_eff_c", "proton_eff_c",11,49,700,502);
	gStyle->SetOptFit(1);
  	proton_eff_c.Range(-1.008394,-80,3.707364,453.3333);
  	proton_eff_c.SetFillColor(0);
  	proton_eff_c.SetBorderMode(0);
  	proton_eff_c.SetBorderSize(2);
  	proton_eff_c.SetLogx();
	proton_eff_c.SetLeftMargin(0.15);
  	proton_eff_c.SetRightMargin(0.15);
  	proton_eff_c.SetBottomMargin(0.15);
  	proton_eff_c.SetFrameLineWidth(2);
  	proton_eff_c.SetFrameBorderMode(0);
  	proton_eff_c.SetFrameBorderSize(2);
  	proton_eff_c.SetFrameLineWidth(2);
  	proton_eff_c.SetFrameBorderMode(0);
  	proton_eff_c.SetFrameBorderSize(2);
    proton_eff->Draw();
    proton_eff_c.Write();
    
    out_file.Close();
}