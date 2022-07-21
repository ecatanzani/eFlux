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

const bool AMS01_positron = true;
const bool AMS02_positron = true;
const bool FERMILAT_positron = true;
const bool CAPRICE_positron = true;
const bool PAMELA_positron = true;
const bool HEAT_positron = true;

TColor *kAMS01_positron = gROOT->GetColor(kOrange+1);
TColor *kAMS02_positron = gROOT->GetColor(kRed);
TColor *kFERMILAT_positron = gROOT->GetColor(kGreen+2);
TColor *kCAPRICE_positron = gROOT->GetColor(kAzure-3);
TColor *kPAMELA_positron = gROOT->GetColor(kMagenta-4);
TColor *kHEAT_positron = gROOT->GetColor(kCyan-6);

int ams01_marker_style = 20;
int ams02_marker_style = 21;
int fermi_marker_style = 22;
int caprice_marker_style = 33;
int pamela_marker_style = 41;
int heat_marker_style = 34;

const char* ams01_positron_legendentry = "AMS01 (2000)";
const char* ams02_positron_legendentry = "AMS02 (2019)";
const char* fermilat_positron_legendentry = "FERMI-LAT (2012)";
const char* caprice_positron_legendentry = "CAPRICE (2000)";
const char* pamela_positron_legendentry = "PAMELA (2013)";
const char* heat_positron_legendentry = "HEAT (2001)";

const char* ams01_positron_path = "data/crdb_AMS01_positrons_1998_E3.root";
const char* ams02_positron_path = "data/ssdc_AMS02_positrons_2011_2017_E3.root";
const char* fermilat_positron_path = "data/ssdc_FERMI_LAT_positrons_2008_2011_E3.root";
const char* caprice_positron_path = "data/crdb_CAPRICE_positrons_1994_E3.root";
const char* pamela_positron_path = "data/ssdc_PAMELA_positrons_2006_2010_E3.root";
const char* heat_positron_path = "data/crdb_HEAT_positrons_1994_1995_E3.root";

const double emin = 0.5;
const double emax = 1000;
const double phi_min = 0;
const double phi_max = 60;

const char* multigraph_name = "flux_mg_E3";

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
    const char* outfile = "flux.root",
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
    
    if (kAMS01_positron) add_to_mg(ams01_positron_path, kAMS01_positron, ams01_marker_style, 1, ams01_positron_legendentry, mgE3, "crdb");
    if (kAMS02_positron) add_to_mg(ams02_positron_path, kAMS02_positron, ams02_marker_style, 1, ams02_positron_legendentry, mgE3, "ssdc");
    if (kFERMILAT_positron) add_to_mg(fermilat_positron_path, kFERMILAT_positron, fermi_marker_style, 1, fermilat_positron_legendentry, mgE3, "ssdc");
    if (CAPRICE_positron) add_to_mg(caprice_positron_path, kCAPRICE_positron, caprice_marker_style, 1, caprice_positron_legendentry, mgE3, "crdb");
    if (PAMELA_positron) add_to_mg(pamela_positron_path, kPAMELA_positron, pamela_marker_style, 1, pamela_positron_legendentry, mgE3, "ssdc");
    if (HEAT_positron) add_to_mg(heat_positron_path, kHEAT_positron, heat_marker_style, 1, heat_positron_legendentry, mgE3, "crdb");

    mgE3->GetXaxis()->SetTitle("Energy [GeV]");
    mgE3->GetYaxis()->SetTitle("E^{3} #times #Phi_{e^{+}} (GeV^{2} [m^{2} sr s]^{-1})");
    
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