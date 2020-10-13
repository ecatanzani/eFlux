#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <memory>

#include "TFile.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TStyle.h"

void getEff(const char* input)
{
	TFile *infile = TFile::Open(input, "READ");
	auto h_geo = static_cast<TH1D*>(infile->Get("h_geometric_cut"));
	auto h_bgo_fid = static_cast<TH1D*>(infile->Get("h_geometric_BGO_fiducial_cut"));
	auto h_maxrms = static_cast<TH1D*>(infile->Get("h_BGOfiducial_maxRms_cut"));
	auto h_track = static_cast<TH1D*>(infile->Get("h_BGOfiducial_track_selection_cut"));
	auto h_psd_stk_match = static_cast<TH1D*>(infile->Get("h_BGOfiducial_psd_stk_match_cut"));
	auto h_psd_charge = static_cast<TH1D*>(infile->Get("h_BGOfiducial_psd_charge_cut"));
	auto h_stk_charge = static_cast<TH1D*>(infile->Get("h_BGOfiducial_stk_charge_cut"));	

	h_geo->SetDirectory(0);
	h_bgo_fid->SetDirectory(0);
	h_maxrms->SetDirectory(0);
	h_track->SetDirectory(0);
	h_psd_stk_match->SetDirectory(0);
	h_psd_charge->SetDirectory(0);
	h_stk_charge->SetDirectory(0);

	h_geo->GetYaxis()->SetTitle("efficiency");
	h_bgo_fid->GetYaxis()->SetTitle("efficiency");
	h_maxrms->GetYaxis()->SetTitle("efficiency");
	h_track->GetYaxis()->SetTitle("efficiency");
	h_psd_stk_match->GetYaxis()->SetTitle("efficiency");
	h_psd_charge->GetYaxis()->SetTitle("efficiency");
	h_stk_charge->GetYaxis()->SetTitle("efficiency");


	infile->Close();
	
	TEfficiency* bgo_fid_eff = nullptr;
	TEfficiency* tr_eff = nullptr;
	TEfficiency* psd_charge_eff = nullptr;
	TEfficiency* stk_charge_eff = nullptr;

	if (TEfficiency::CheckConsistency(*h_bgo_fid, *h_geo))
                bgo_fid_eff = new TEfficiency(*h_bgo_fid, *h_geo);
	if (TEfficiency::CheckConsistency(*h_track, *h_maxrms))
		tr_eff = new TEfficiency(*h_track, *h_maxrms);
	if (TEfficiency::CheckConsistency(*h_psd_charge, *h_psd_stk_match))
                psd_charge_eff = new TEfficiency(*h_psd_charge, *h_psd_stk_match);
	if (TEfficiency::CheckConsistency(*h_stk_charge, *h_psd_charge))
                stk_charge_eff = new TEfficiency(*h_stk_charge, *h_psd_charge);

	bgo_fid_eff->SetName("bgo_fid_eff");
	bgo_fid_eff->SetTitle("BGO fiducial");
	tr_eff->SetName("tr_eff");
	tr_eff->SetTitle("Track selection");
	psd_charge_eff->SetName("psd_charge_eff");
	psd_charge_eff->SetTitle("PSD charge");
	stk_charge_eff->SetName("stk_charge_eff");
	stk_charge_eff->SetTitle("STK charge");

	TFile* outfile = TFile::Open("efficiencies.root", "RECREATE");
	bgo_fid_eff->Write();
	tr_eff->Write();
	psd_charge_eff->Write();
	stk_charge_eff->Write();
	
	TCanvas ceff_all("ceff_all", "ceff_all",11,49,700,502);
	gStyle->SetOptFit(1);
  	ceff_all.Range(-1.008394,-80,3.707364,453.3333);
  	ceff_all.SetFillColor(0);
  	ceff_all.SetBorderMode(0);
  	ceff_all.SetBorderSize(2);
  	ceff_all.SetLogx();
	ceff_all.SetLeftMargin(0.15);
  	ceff_all.SetRightMargin(0.15);
  	ceff_all.SetBottomMargin(0.15);
  	ceff_all.SetFrameLineWidth(2);
  	ceff_all.SetFrameBorderMode(0);
  	ceff_all.SetFrameBorderSize(2);
  	ceff_all.SetFrameLineWidth(2);
  	ceff_all.SetFrameBorderMode(0);
  	ceff_all.SetFrameBorderSize(2);

	bgo_fid_eff->Draw();
	tr_eff->Draw("same");
	psd_charge_eff->Draw("same");
	stk_charge_eff->Draw("same");
	ceff_all.Write();

	outfile->Close();
}

void getAcc(const char* input)
{

	TFile *file = TFile::Open(input, "READ");
	auto acceptance = static_cast<TH1D*>(file->Get("acceptance"));
	acceptance->SetDirectory(0);
	file->Close();

	acceptance->GetXaxis()->SetTitle("Energy (GeV)");
	acceptance->GetYaxis()->SetTitle("acceptance (m^{2} sr)");
	TFile* outfile = TFile::Open("effeptance.root", "RECREATE");
	TCanvas ceff_all("ceff_all", "ceff_all",11,49,700,502);
        gStyle->SetOptFit(1);
        ceff_all.Range(-1.008394,-80,3.707364,453.3333);
        ceff_all.SetFillColor(0);
        ceff_all.SetBorderMode(0);
        ceff_all.SetBorderSize(2);
        ceff_all.SetLogx();
        ceff_all.SetLeftMargin(0.15);
        ceff_all.SetRightMargin(0.15);
        ceff_all.SetBottomMargin(0.15);
        ceff_all.SetFrameLineWidth(2);
        ceff_all.SetFrameBorderMode(0);
        ceff_all.SetFrameBorderSize(2);
        ceff_all.SetFrameLineWidth(2);
        ceff_all.SetFrameBorderMode(0);
        ceff_all.SetFrameBorderSize(2);
	acceptance->Draw();
	ceff_all.Write();
	outfile->Close();
}
