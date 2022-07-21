#include <memory>
#include <vector>
#include <string>
#include <iostream>

#include "TROOT.h"
#include "TRint.h"
#include "TAxis.h"
#include "TKey.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TLegendEntry.h"
#include "TApplication.h"
#include "TGraphAsymmErrors.h"

const bool AMS01_allelectron = true;
const bool AMS02_allelectron = true;
const bool DAMPE_allelectron = true;
const bool FERMILAT_allelectron = true;
const bool CALET_allelectron = true;
const bool CAPRICE_allelectron = true;
const bool PAMELA_allelectron = true;
const bool HEAT_allelectron = true;
const bool HESS_allelectron = true;

TColor *kAMS01_allelectron = gROOT->GetColor(kOrange+1);
TColor *kAMS02_allelectron = gROOT->GetColor(kOrange+8);
TColor *kDAMPE_allelectron = gROOT->GetColor(kRed);
TColor *kFERMILAT_allelectron = gROOT->GetColor(kGreen+2);
TColor *kCALET_allelectron = gROOT->GetColor(kBlue+1);
TColor *kCAPRICE_allelectron = gROOT->GetColor(kBlue-9);
TColor *kPAMELA_allelectron = gROOT->GetColor(kMagenta-4);
TColor *kHEAT_allelectron = gROOT->GetColor(kCyan-6);
TColor *kHESS_allelectron = gROOT->GetColor(kGray+1);
TColor *kHESS_allelectron_le = gROOT->GetColor(kGray+1);

int ams01_marker_style = 21;
int ams02_marker_style = 22;
int fermi_marker_style = 23;
int calet_marker_style = 29;
int dampe_marker_style = 20;
int caprice_marker_style = 33;
int pamela_marker_style = 41;
int heat_marker_style = 34;
int heat_hess_style = 45;

const char* ams01_allelectron_legendentry = "AMS01 (2007)";
const char* ams02_allelectron_legendentry = "AMS02 (2021)";
const char* dampe_allelectron_legendentry = "DAMPE (2017)";
const char* fermilat_allelectron_legendentry = "FERMI-LAT (2017)";
const char* calet_allelectron_legendentry = "CALET (20118)";
const char* caprice_allelectron_legendentry = "CAPRICE (2000)";
const char* pamela_allelectron_legendentry = "PAMELA (2013)";
const char* heat_allelectron_legendentry = "HEAT (2001)";
const char* hess_allelectron_legendentry = "HESS (2008-2009)";

/*
const char* hess_allelectron_le_legendentry = "HESS (2009)";
const char* hess_allelectron_legendentry = "HESS (2008)";
*/

const char* ams01_allelectron_path = "data/crdb_AMS01_alle_1998_E3.root";
const char* ams02_allelectron_path = "data/crdb_AMS02_alle_2011_2018_E3.root";
const char* dampe_allelectron_path = "data/crdb_DAMPE_alle_2015_2017_E3.root";
const char* fermilat_allelectron_path = "data/ssdc_FERMI_LAT_alle_2008_2015_E3.root";
const char* calet_allelectron_path = "data/crdb_CALET_alle_2015_2017_E3.root";
const char* caprice_allelectron_path = "data/crdb_CAPRICE_alle_1994_E3.root";
const char* pamela_allelectron_path = "data/crdb_PAMELA_alle_2006_2009_E3.root";
const char* heat_allelectron_path = "data/ssdc_HEAT_alle_1994_1995_E3.root";
const char* hess_allelectron_path = "data/hess_he_syst_allele.root";
const char* hess_allelectron_path_le = "data/hess_le_syst_allele.root";

const double emin = 0.5;
const double emax = 6000;
const double phi_min = 0;
const double phi_max = 300;

void add_to_mg(
    const char* path, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3,
    const bool draw_syst_band,
    const char* provider,
    const std::vector<std::string> &filter=std::vector<std::string>());
std::shared_ptr<TMultiGraph> buildflux(const char* name, const bool hess=false);
void write(const char* outfile, std::shared_ptr<TMultiGraph> mgE3, std::shared_ptr<TMultiGraph> mgE3_hess);
void drawflux(std::shared_ptr<TMultiGraph> mgr, std::shared_ptr<TMultiGraph> mgr_hess);
//void drawflux(const char* infile);

void draw(
    const char* outfile = "flux.root",
    const bool plot = true)
{
    auto mgE3 = buildflux("flux_mg_E3");
    auto mgE3_hess = buildflux("flux_mg_E_hess", true);
    write(outfile, mgE3, mgE3_hess);
    if (plot) drawflux(mgE3, mgE3_hess);
}

std::shared_ptr<TMultiGraph> buildflux(const char* name, const bool hess)
{
    std::shared_ptr<TMultiGraph> mgE3 = std::make_shared<TMultiGraph>();
    mgE3->SetName(name);
    
    bool draw_band {false};

    if (!hess) {
        if (HESS_allelectron) {
            add_to_mg(hess_allelectron_path_le, kHESS_allelectron, heat_hess_style, 1, hess_allelectron_legendentry, mgE3, draw_band, "other");
            add_to_mg(hess_allelectron_path, kHESS_allelectron, heat_hess_style, 1, hess_allelectron_legendentry, mgE3, draw_band, "other");
        }
        if (AMS01_allelectron) add_to_mg(ams01_allelectron_path, kAMS01_allelectron, ams01_marker_style, 1, ams01_allelectron_legendentry, mgE3, draw_band, "crdb");
        if (AMS02_allelectron) add_to_mg(ams02_allelectron_path, kAMS02_allelectron, ams02_marker_style, 1, ams02_allelectron_legendentry, mgE3, draw_band, "crdb");
        if (FERMILAT_allelectron) add_to_mg(fermilat_allelectron_path, kFERMILAT_allelectron, fermi_marker_style, 1, fermilat_allelectron_legendentry, mgE3, draw_band, "ssdc");
        if (CALET_allelectron) add_to_mg(calet_allelectron_path, kCALET_allelectron, calet_marker_style, 1, calet_allelectron_legendentry, mgE3, draw_band, "crdb");
        if (CAPRICE_allelectron) add_to_mg(caprice_allelectron_path, kCAPRICE_allelectron, caprice_marker_style, 1, caprice_allelectron_legendentry, mgE3, draw_band, "crdb");
        if (PAMELA_allelectron) add_to_mg(pamela_allelectron_path, kPAMELA_allelectron, pamela_marker_style, 1, pamela_allelectron_legendentry, mgE3, draw_band, "crdb");
        if (HEAT_allelectron) add_to_mg(heat_allelectron_path, kHEAT_allelectron, heat_marker_style, 1, heat_allelectron_legendentry, mgE3, draw_band, "ssdc");
        if (DAMPE_allelectron) add_to_mg(dampe_allelectron_path, kDAMPE_allelectron, dampe_marker_style, 1, dampe_allelectron_legendentry, mgE3, draw_band, "crdb");
        
    }
    else {
        if (HESS_allelectron) {
            draw_band = true;
            add_to_mg(hess_allelectron_path_le, kHESS_allelectron, heat_hess_style, 1, hess_allelectron_legendentry, mgE3, draw_band, "other");
            add_to_mg(hess_allelectron_path, kHESS_allelectron, heat_hess_style, 1, hess_allelectron_legendentry, mgE3, draw_band, "other");
        }
    }

    mgE3->GetXaxis()->SetTitle("Energy [GeV]");
    mgE3->GetYaxis()->SetTitle("E^{3} #times #Phi_{e^{-}+e^{+}} (GeV^{2} [m^{2} sr s]^{-1})");
    
    return mgE3;
}

void add_to_mg(
    const char* path, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3,
    const bool draw_syst_band,
    const char* provider,
    const std::vector<std::string> &filter)
{
    TFile* inFile = TFile::Open(path, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << path << "]\n";
        exit(100);
    }

    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto exclude_obj = [&filter] (TGraphAsymmErrors* tmp_obj) -> bool
        {
            bool exclude = false;
            if (filter.size())
            {
                exclude = true;
                for (auto& elm : filter)
                    if (!strcmp(tmp_obj->GetName(), elm.c_str()))
                        exclude = false;
            }
            return exclude;
        };

        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors") || !strcmp(key->GetClassName(), "TGraphErrors"))
        {
            auto tmp_obj = static_cast<TGraphAsymmErrors *>(key->ReadObj());
            if (exclude_obj(tmp_obj))
                continue;
            
            if (!std::strcmp(provider, "crdb") && std::strcmp(tmp_obj->GetName(), "gr_exp1_errtot"))
                continue;

            tmp_obj->SetLineColor(mycolor->GetNumber());
            tmp_obj->SetMarkerColor(mycolor->GetNumber());
            if (std::string(tmp_obj->GetName()).find("upper_limit") != std::string::npos) marker_style = 23;
            tmp_obj->SetMarkerStyle(marker_style);
            tmp_obj->SetMarkerSize(marker_size);
            tmp_obj->SetLineWidth(2);
            tmp_obj->SetTitle(gr_title);

            if (draw_syst_band) {
                tmp_obj->SetFillColor(kGray);
                tmp_obj->SetFillStyle(3001);
            }

            mgE3->Add(tmp_obj);
        }
    }
}

void write(
    const char* outfile, 
    std::shared_ptr<TMultiGraph> mgE3,
    std::shared_ptr<TMultiGraph> mgE3_hess)
{
    TFile* myoutfile = TFile::Open(outfile, "RECREATE");
    if (myoutfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outfile << "]\n";
        exit(100);
    }

    mgE3->Write();
    mgE3_hess->Write();
    myoutfile->Close();
}

/*
void drawflux(const char* infile)
{
    TApplication theApp("App", 0, 0);

    TFile* inFile = TFile::Open(infile, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << infile << "]\n";
        exit(100);
    }

    std::unique_ptr<TMultiGraph> mgr = std::unique_ptr<TMultiGraph>(static_cast<TMultiGraph*>(inFile->Get(multigraph_name)));
    TCanvas cflux("Flux", "Flux");
    cflux.cd();
    //mgr->Draw("AP");
    mgr->Draw("P");
    gPad->Modified(); 
    gPad->Update();
    mgr->GetXaxis()->SetLimits(emin, emax);
    mgr->GetYaxis()->SetLimits(phi_min, phi_max);
    mgr->SetMinimum(phi_min);
    mgr->SetMaximum(phi_max);
    gPad->Modified(); 
    gPad->Update();

    theApp.Run(); 
    gSystem->ProcessEvents();
}
*/

void drawflux(std::shared_ptr<TMultiGraph> mgr, std::shared_ptr<TMultiGraph> mgr_hess)
{
    TApplication theApp("App", 0, 0);
    
    gStyle->SetLineWidth(3);

    TCanvas cflux("Flux", "Flux", 900, 700);
    cflux.cd();

    mgr_hess->Draw("a3");

    gPad->Modified(); 
    gPad->Update();
    mgr_hess->GetXaxis()->SetLimits(emin, emax);
    mgr_hess->GetYaxis()->SetLimits(phi_min, phi_max);
    mgr_hess->SetMinimum(phi_min);
    mgr_hess->SetMaximum(phi_max);
    gPad->Modified(); 
    gPad->Update();
    cflux.SetLogx();
    //cflux.SetLogy();
    cflux.SetTicks();
    gPad->Modified(); 
    gPad->Update();
    mgr_hess->GetXaxis()->SetLabelSize(0.045);
    mgr_hess->GetXaxis()->SetTitleSize(0.045);
    mgr_hess->GetXaxis()->SetTitleOffset(1.3);
    mgr_hess->GetYaxis()->SetLabelSize(0.045);
    mgr_hess->GetYaxis()->SetTitleSize(0.045);
    mgr_hess->GetYaxis()->SetTitleOffset(1.6);
    gPad->Modified(); 
    gPad->Update();

    mgr->Draw("P");
    //mgr_hess->Draw("3P");

    auto legend = cflux.BuildLegend();
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    auto primitives = legend->GetListOfPrimitives();
    for (auto primitiveObj :  *primitives)
    {
        auto primitive = (TLegendEntry*)primitiveObj;
        primitive->SetOption("p");
    }
    gPad->Modified(); 
    gPad->Update();

    theApp.Run(); 
    gSystem->ProcessEvents();
}