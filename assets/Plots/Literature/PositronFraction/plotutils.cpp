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

const bool AMS01_ps = true;
const bool AMS02_ps = true;
const bool FERMILAT_ps = true;
const bool PAMELA_ps = true;
const bool CAPRICE_ps = true;
const bool HEAT_ps = true;
const bool TS93_ps = true;

TColor *kAMS01_ps = gROOT->GetColor(kOrange+1);
TColor *kAMS02_ps = gROOT->GetColor(kRed);
TColor *kFERMILAT_ps = gROOT->GetColor(kGreen+2);
TColor *kCAPRICE_ps = gROOT->GetColor(kAzure-3);
TColor *kCAPRICE_ps_new = gROOT->GetColor(kViolet-7);
TColor *kPAMELA_ps = gROOT->GetColor(kMagenta-4);
TColor *kHEAT_ps = gROOT->GetColor(kCyan-6);
TColor *kTS93_ps = gROOT->GetColor(kYellow-6);

int ams01_marker_style = 20;
int ams02_marker_style = 21;
int fermi_marker_style = 22;
int caprice_marker_style = 33;
int pamela_marker_style = 41;
int heat_marker_style = 34;
int ts93_marker_style = 29;

const char* ams01_ps_legendentry = "AMS01 (2007)";
const char* ams02_ps_legendentry = "AMS02 (2021)";
const char* fermilat_ps_legendentry = "FERMI-LAT (2012)";
const char* pamela_ps_legendentry = "PAMELA (2013)";
const char* caprice_ps_1994_legendentry = "CAPRICE (2000)";
const char* caprice_ps_1998_legendentry = "CAPRICE (2001)";
const char* heat_ps_legendentry = "HEAT (1997)";
const char* ts93_ps_legendentry = "TS93 (1996)";

const char* ams01_ps_path = "data/crdb_AMS01_ps_1998.root";
const char* ams02_ps_path = "data/ssdc_AMS02_ps_2011_2018.root";
const char* fermilat_ps_path = "data/ssdc_FERMI_LAT_ps_2008_2011.root";
const char* pamela_ps_path = "data/ssdc_PAMELA_ps_2006_2010.root";
const char* caprice_1994_ps_path = "data/crdb_CAPRICE_ps_1994.root";
const char* caprice_1998_ps_path = "data/crdb_CAPRICE_ps_1998.root";
const char* heat_ps_path = "data/crdb_HEAT_ps_1994_1995.root";
const char* ts93_ps_path = "data/ssdc_TS93_ps_1993.root";

const double emin = 0.5;
const double emax = 1000;
const double phi_min = 0;
const double phi_max = 0.4;

const char* multigraph_name = "flux_mg";

void add_to_mg(
    const char* path, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3,
    const char* provider,
    const std::vector<std::string> &filter=std::vector<std::string>());
std::shared_ptr<TMultiGraph> buildflux();
void write(const char* outfile, std::shared_ptr<TMultiGraph> mgE3);
void drawflux(std::shared_ptr<TMultiGraph> mgr);
void drawflux(const char* infile);

void draw(
    const char* outfile = "ps.root",
    const bool plot = true)
{
    auto mgE3 = buildflux();
    write(outfile, mgE3);
    if (plot) drawflux(mgE3);
}

std::shared_ptr<TMultiGraph> buildflux()
{
    std::shared_ptr<TMultiGraph> mgE3 = std::make_shared<TMultiGraph>();
    mgE3->SetName(multigraph_name);
    
    if (kAMS01_ps) add_to_mg(ams01_ps_path, kAMS01_ps, ams01_marker_style, 1, ams01_ps_legendentry, mgE3, "crdb");
    if (kAMS02_ps) add_to_mg(ams02_ps_path, kAMS02_ps, ams02_marker_style, 1, ams02_ps_legendentry, mgE3, "ssdc");
    if (kFERMILAT_ps) add_to_mg(fermilat_ps_path, kFERMILAT_ps, fermi_marker_style, 1, fermilat_ps_legendentry, mgE3, "ssdc");
    if (CAPRICE_ps) {
        add_to_mg(caprice_1994_ps_path, kCAPRICE_ps, caprice_marker_style, 1, caprice_ps_1994_legendentry, mgE3, "crdb");
        add_to_mg(caprice_1998_ps_path, kCAPRICE_ps_new, caprice_marker_style, 1, caprice_ps_1998_legendentry, mgE3, "crdb");
    }
    if (PAMELA_ps) add_to_mg(pamela_ps_path, kPAMELA_ps, pamela_marker_style, 1, pamela_ps_legendentry, mgE3, "ssdc");
    if (HEAT_ps) add_to_mg(heat_ps_path, kHEAT_ps, heat_marker_style, 1, heat_ps_legendentry, mgE3, "crdb");
    if (TS93_ps) add_to_mg(ts93_ps_path, kTS93_ps, ts93_marker_style, 1, ts93_ps_legendentry, mgE3, "ssdc");

    mgE3->GetXaxis()->SetTitle("Energy [GeV]");
    mgE3->GetYaxis()->SetTitle("Positron fraction");
    
    return mgE3;
}

void add_to_mg(
    const char* path, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3,
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
            mgE3->Add(tmp_obj);
        }
    }
}

void write(
    const char* outfile, 
    std::shared_ptr<TMultiGraph> mgE3)
{
    TFile* myoutfile = TFile::Open(outfile, "RECREATE");
    if (myoutfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outfile << "]\n";
        exit(100);
    }

    mgE3->Write();
    myoutfile->Close();
}

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

void drawflux(std::shared_ptr<TMultiGraph> mgr)
{
    TApplication theApp("App", 0, 0);
    
    gStyle->SetLineWidth(3);

    TCanvas cflux("Flux", "Flux", 700, 700);
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
    cflux.SetLogx();
    //cflux.SetLogy();
    cflux.SetTicks();
    gPad->Modified(); 
    gPad->Update();
    mgr->GetXaxis()->SetLabelSize(0.045);
    mgr->GetXaxis()->SetTitleSize(0.045);
    mgr->GetXaxis()->SetTitleOffset(1.3);
    mgr->GetYaxis()->SetLabelSize(0.045);
    mgr->GetYaxis()->SetTitleSize(0.045);
    mgr->GetYaxis()->SetTitleOffset(1.6);
    gPad->Modified(); 
    gPad->Update();

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