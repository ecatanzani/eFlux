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

const bool ACE = true;
const bool ISOMAX = true;
const bool VOYAGER_1991 = true;
const bool VOYAGER_1998 = true;

TColor *kACE = gROOT->GetColor(kRed);
TColor *kISOMAX = gROOT->GetColor(kBlue);
TColor *kVOYAGER_1991 = gROOT->GetColor(kGreen);
TColor *kVOYAGER_1998 = gROOT->GetColor(kOrange);

int ace_marker_style = 105;
int isomax_marker_style = 20;
int voyager_1991_marker_style = 21;
int voyager_1998_marker_style = 22;

const char* ace_legendentry = "ACE (1999)";
const char* isomax_legendentry = "ISOMAX (1998)";
const char* voyager1991_legendentry = "Voyager 1&2 (1991)";
const char* voyager1998_legendentry = "Voyager 1&2 (1998)";

const double min_value = 0;
const double max_value = 0.5;
const double emin = 1e-2;
const double emax = 10;

const char* multigraph_name = "mg_10Be9Be";
const char* data = "data/in2p3_10Be_9Be.root";

std::vector<const char*> ace_entry_name = {"gr_exp1_errtot", "gr_exp2_errtot"};
const char* isomax_entry_name = "gr_exp3_errtot";
const char* voyager1991_entry_name = "gr_exp4_errtot";
const char* voyager1998_entry_name = "gr_exp5_errtot";

void add_to_mg(
    const char* name, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const double line_width,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3);
void add_to_mg(
    std::vector<const char*> name, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const double line_width,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3);
std::shared_ptr<TMultiGraph> buildflux();
void write(const char* outfile, std::shared_ptr<TMultiGraph> mgE3);
void drawflux(std::shared_ptr<TMultiGraph> mgr);
void drawflux(const char* infile);

void draw(
    const char* outfile = "10Be_9Be.root",
    const bool plot = true)
{
    auto mg = buildflux();
    write(outfile, mg);
    if (plot) drawflux(mg);
}

std::shared_ptr<TMultiGraph> buildflux()
{
    std::shared_ptr<TMultiGraph> mg = std::make_shared<TMultiGraph>();
    mg->SetName(multigraph_name);
    
    if (ACE) add_to_mg(ace_entry_name, kACE, ace_marker_style, 2, 2, ace_legendentry, mg);
    if (ISOMAX) add_to_mg(isomax_entry_name, kISOMAX, isomax_marker_style, 2, 2, isomax_legendentry, mg);
    if (VOYAGER_1991) add_to_mg(voyager1991_entry_name, kVOYAGER_1991, voyager_1991_marker_style, 2, 2, voyager1991_legendentry, mg);
    if (VOYAGER_1998) add_to_mg(voyager1998_entry_name, kVOYAGER_1998, voyager_1998_marker_style, 2, 2, voyager1998_legendentry, mg);

    mg->GetXaxis()->SetTitle("Energy (GeV/n)");
    mg->GetYaxis()->SetTitle("^{10}Be/^{9}Be");
    
    return mg;
}

void add_to_mg(
    const char* name, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const double line_width,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mg)
{
    TFile* inFile = TFile::Open(data, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << data << "]\n";
        exit(100);
    }


    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);

        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
        {
            auto tmp_obj = static_cast<TGraphAsymmErrors *>(key->ReadObj());
            auto iswhatiwant = [&name](TGraphAsymmErrors* obj) -> bool { if(!strcmp(obj->GetName(), name)) return true; else return false; };
            if (!iswhatiwant(tmp_obj))
                continue;
            tmp_obj->SetLineColor(mycolor->GetNumber());
            tmp_obj->SetMarkerColor(mycolor->GetNumber());
            if (std::string(tmp_obj->GetName()).find("upper_limit") != std::string::npos) marker_style = 23;
            tmp_obj->SetMarkerStyle(marker_style);
            tmp_obj->SetMarkerSize(marker_size);
            tmp_obj->SetLineWidth(line_width);
            tmp_obj->SetTitle(gr_title);
            mg->Add(tmp_obj);
        }
    }
}

void add_to_mg(
    std::vector<const char*> name, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const double line_width,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mg)
{
    TFile* inFile = TFile::Open(data, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << data << "]\n";
        exit(100);
    }


    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);

        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
        {
            auto tmp_obj = static_cast<TGraphAsymmErrors *>(key->ReadObj());
            auto iswhatiwant = [&name](TGraphAsymmErrors* obj) -> bool { 
                bool found = false;
                for (auto& elm : name)
                    if(!strcmp(obj->GetName(), elm))
                    {
                        found = true;
                        break;
                    }
                return found; 
            };

            if (!iswhatiwant(tmp_obj))
                continue;
            tmp_obj->SetLineColor(mycolor->GetNumber());
            tmp_obj->SetMarkerColor(mycolor->GetNumber());
            tmp_obj->SetMarkerStyle(marker_style);
            tmp_obj->SetMarkerSize(marker_size);
            tmp_obj->SetLineWidth(line_width);
            tmp_obj->SetTitle(gr_title);
            mg->Add(tmp_obj);
        }
    }
}

void write(
    const char* outfile, 
    std::shared_ptr<TMultiGraph> mg)
{
    TFile* myoutfile = TFile::Open(outfile, "RECREATE");
    if (myoutfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outfile << "]\n";
        exit(100);
    }

    mg->Write();
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
    mgr->GetYaxis()->SetLimits(min_value, max_value);
    mgr->SetMinimum(min_value);
    mgr->SetMaximum(max_value);
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
    mgr->GetYaxis()->SetLimits(min_value, max_value);
    mgr->SetMinimum(min_value);
    mgr->SetMaximum(max_value);
    gPad->Modified(); 
    gPad->Update();
    cflux.SetLogx();
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