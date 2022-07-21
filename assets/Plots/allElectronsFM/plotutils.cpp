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

const bool DAMPE_electron = true;
const bool AMS02_electron = true;
const bool CALET_electron = true;
const bool FERMILAT_electron = true;
const bool ATIC2_electron = false;
const bool HESS_electron = true;

TColor *kAMS02_electron = gROOT->GetColor(kMagenta+1);
TColor *kCALET_electron = gROOT->GetColor(kGreen+2);
TColor *kDAMPE_electron = gROOT->GetColor(kBlue);
TColor *kDAMPE_updated_electron = gROOT->GetColor(kRed);
TColor *kFERMILAT_electron = gROOT->GetColor(kOrange-3);
TColor *kHESS_electron = gROOT->GetColor(kGray+1);
TColor *kATIC2_electron = gROOT->GetColor(kBlack);

int ams02_marker_style = 20;
int dampe_marker_style = 20;
int calet_marker_style = 20;
int fermi_marker_style = 20;
int atic2_marker_style = 20;
int hess_marker_style = 20;

const char* ams02_allele_legendentry = "e^{-}+e^{+} AMS02 (2021)";
const char* calet_allele_legendentry = "e^{-}+e^{+} CALET (2018)";
const char* fermilat_allele_legendentry = "e^{-}+e^{+} FERMI-LAT (2017)";
const char* dampe_allele_legendentry = "e^{-}+e^{+} DAMPE (2017)";
const char* dampe_allele_updated_legendentry = "e^{-}+e^{+} DAMPE (This Work)";
const char* atic2_allele_legendentry = "e^{-}+e^{+} ATIC-2";
const char* hess_allele_legendentry = "e^{-}+e^{+} HESS (2008-2009)";

const char* ams02_electron_path = "data/ssdc_AMS02_e_E3.root";
const char* dampe_electron_path = "data/ssdc_DAMPE_e_E3.root";
const char* dampe_updated_electron_path = "../../../fluxBuilder/flux_310522.root";

//const char* calet_electron_path = "data/crdb_CALET_alle_2015_2017_E3.root";
const char* calet_electron_path = "data/ssdc_CALET_e_E3.root";

const char* fermilat_electron_path = "data/ssdc_FERMILAT_e_E3.root";
const char* atic2_electron_path = "data/ssdc_ATIC2_e_E3.root";
const char* hess_allelectron_path = "data/hess_he_syst_allele.root";
const char* hess_allelectron_path_le = "data/hess_le_syst_allele.root";

const double emin = 20;
const double emax = 1e+4;
const double phi_min = 0;
const double phi_max = 300;

const char* multigraph_name = "flux_mg_E3";

std::vector<std::string> dample_filter = {"gr_flux_E3"};

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
    auto mgE3_hess = buildflux("flux_mg_E3_hess", true);
    write(outfile, mgE3, mgE3_hess);
    if (plot) drawflux(mgE3, mgE3_hess);
}

std::shared_ptr<TMultiGraph> buildflux(const char* name, const bool hess)
{
    std::shared_ptr<TMultiGraph> mgE3 = std::make_shared<TMultiGraph>();
    mgE3->SetName(name);
    
    bool draw_band {false};

    if (!hess) {
        if (HESS_electron) {
            add_to_mg(hess_allelectron_path, kHESS_electron, hess_marker_style, 1, hess_allele_legendentry, mgE3, draw_band, "other");
            add_to_mg(hess_allelectron_path_le, kHESS_electron, hess_marker_style, 1, hess_allele_legendentry, mgE3, draw_band, "other");
        }

        if (AMS02_electron) add_to_mg(ams02_electron_path, kAMS02_electron, ams02_marker_style, 1, ams02_allele_legendentry, mgE3, draw_band, "ssdc");
        if (CALET_electron) add_to_mg(calet_electron_path, kCALET_electron, calet_marker_style, 1, calet_allele_legendentry, mgE3, draw_band, "ssdc");
        if (FERMILAT_electron) add_to_mg(fermilat_electron_path, kFERMILAT_electron, fermi_marker_style, 1, fermilat_allele_legendentry, mgE3, draw_band, "ssdc");
        
        if (ATIC2_electron) add_to_mg(atic2_electron_path, kATIC2_electron, atic2_marker_style, 1, atic2_allele_legendentry, mgE3, draw_band, "ssdc");

        if (DAMPE_electron) {
            add_to_mg(dampe_electron_path, kDAMPE_electron, dampe_marker_style, 1, dampe_allele_legendentry, mgE3, draw_band, "ssdc");
            add_to_mg(dampe_updated_electron_path, kDAMPE_updated_electron, dampe_marker_style, 1, dampe_allele_updated_legendentry, mgE3, draw_band, "ssdc", dample_filter);
        }
    }
    else {
        if (HESS_electron) {
            draw_band = true;
            add_to_mg(hess_allelectron_path, kHESS_electron, hess_marker_style, 1, hess_allele_legendentry, mgE3, draw_band, "other");
            add_to_mg(hess_allelectron_path_le, kHESS_electron, hess_marker_style, 1, hess_allele_legendentry, mgE3, draw_band, "other");
        }
    }

    mgE3->GetXaxis()->SetTitle("Energy [GeV]");
    mgE3->GetYaxis()->SetTitle("E^{3} #times #Phi (GeV^{2} [m^{2} sr s]^{-1})");
    
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
    mgr->Draw("AP");
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