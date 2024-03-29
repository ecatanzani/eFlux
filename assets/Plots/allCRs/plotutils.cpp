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

const bool AMS02_proton = true;
const bool CALET_proton = true;
const bool DAMPE_proton = true;
const bool DAMPE_electron = false;
const bool AMS02_electron = true;
const bool AMS02_positron = true;
const bool AMS02_pbar = true;
const bool AMS02_B = true;
const bool AMS02_C = true;
const bool AMS02_O = true;
const bool CALET_electron = false;
const bool FERMILAT_electron = false;
const bool AUGER_allparticle = false;
const bool KASCADEgrande_allparticle = false;
const bool ATIC2_proton = false;
const bool AMS02_helium = true;

TColor *kAMS02_proton = gROOT->GetColor(kRed+3);
TColor *kCALET_proton = gROOT->GetColor(kOrange+1);
TColor *kDAMPE_proton = gROOT->GetColor(kRed);
TColor *kAMS02_electron = gROOT->GetColor(kBlue-2);
TColor *kAMS02_positron = gROOT->GetColor(kBlue-10);
TColor *kAMS02_pbar = gROOT->GetColor(kRed-9);
TColor *kAMS02_B = gROOT->GetColor(kMagenta-2);
TColor *kAMS02_C = gROOT->GetColor(kGray+2);
TColor *kAMS02_O = gROOT->GetColor(kOrange-6);
TColor *kCALET_electron = gROOT->GetColor(kCyan-6);
TColor *kDAMPE_electron = gROOT->GetColor(kBlue);
TColor *kFERMILAT_electron = gROOT->GetColor(kAzure-4);
TColor *kAUGER_allparticle = gROOT->GetColor(kBlack);
TColor *kKASCADEgrande_allparticle = gROOT->GetColor(kGray+2);
TColor *kATIC2_proton = gROOT->GetColor(kMagenta+2);
TColor *kAMS02_helium = gROOT->GetColor(kGreen+2);

int ams02_marker_style = 107;
int dampe_marker_style = 108;
int calet_marker_style = 33;
int fermi_marker_style = 45;
int auger_marker_style = 108;
int kascadegrande_marker_style = 20;
int atic2_marker_style = 110;

const char* ams02_proton_legendentry = "Protons (AMS02)";
const char* calet_proton_legendentry = "Protons (CALET)";
const char* dampe_proton_legendentry = "Protons (DAMPE)";
const char* ams02_electron_legendentry = "e^{-} (AMS02)";
const char* ams02_positron_legendentry = "e^{+} (AMS02)";
const char* ams02_pbar_legendentry = "#bar{p} (AMS02)";
const char* ams02_B_legendentry = "B (AMS02)";
const char* ams02_C_legendentry = "C (AMS02)";
const char* ams02_O_legendentry = "O (AMS02)";
const char* ams02_allele_legendentry = "e^{-}+e^{+} (AMS02)";
const char* calet_allele_legendentry = "e^{-}+e^{+} (CALET)";
const char* fermilat_allele_legendentry = "e^{-}+e^{+} (FERMI-LAT)";
const char* dampe_allele_legendentry = "e^{-}+e^{+} (DAMPE)";
const char* auger_allparticle_legendentry = "All Particles (AUGER)";
const char* kascadegrande_allparticle_legendentry = "All Particles (KASCADE-Grande)";
const char* atic2_proton_legendentry = "Protons (ATIC-2)";
const char* ams02_helium_legendentry = "Helium (AMS02)";

const char* ams02_proton_path = "data/ssdc_AMS02_p_E3.root";
const char* ams02_electron_path = "data/ssdc_AMS02_e_E3.root";
const char* ams02_positron_path = "data/ssdc_AMS02_positron_E3.root";
const char* ams02_pbar_path = "data/in2p3_crdb_AMS02_pbar_E3.root";
const char* ams02_B_path = "data/ssdc_AMS02_B_E3.root";
const char* ams02_C_path = "data/ssdc_AMS02_C_E3.root";
const char* ams02_O_path = "data/ssdc_AMS02_O_E3.root";
const char* calet_proton_path = "data/ssdc_CALET_p_E3.root";
const char* dampe_proton_path = "data/ssdc_DAMPE_p_E3.root";
const char* dampe_electron_path = "data/ssdc_DAMPE_e_E3.root";
const char* calet_electron_path = "data/ssdc_CALET_e_E3.root";
const char* fermilat_electron_path = "data/ssdc_FERMILAT_e_E3.root";
const char* auger_allparticle_path = "data/in2p3_crdb_AUGER_allparticle_E3.root";
const char* kascadegrande_allparticle_path = "data/in2p3_crdb_KASCADEgrande_allparticle_E3.root";
const char* atic2_proton_path = "data/ssdc_ATIC2_p_E3.root";
const char* ams02_helium_path = "data/ssdc_AMS02_He_E3.root";

std::vector<std::string> auger_filter = {"gr_exp1_errtot", "gr_exp1_upper_limit"};
std::vector<std::string> kascadegrande_filter = {"gr_exp1_errtot"};
std::vector<std::string> ams02_pbar_filter = {"gr_exp1_errtot"};

const double emin = 1;
const double emax = 1e+5;
const double phi_min = 1e-1;
const double phi_max = 1e+6;

const char* multigraph_name = "flux_mg_E3";

void add_to_mg(
    const char* path, 
    TColor* mycolor, 
    int marker_style,
    const double marker_size,
    const char* gr_title,
    std::shared_ptr<TMultiGraph> mgE3,
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
    
    if (AMS02_helium) add_to_mg(ams02_helium_path, kAMS02_helium, ams02_marker_style, 1, ams02_helium_legendentry, mgE3);
    if (AMS02_proton) add_to_mg(ams02_proton_path, kAMS02_proton, ams02_marker_style, 1, ams02_proton_legendentry, mgE3);
    if (CALET_proton) add_to_mg(calet_proton_path, kCALET_proton, calet_marker_style, 1.3, calet_proton_legendentry, mgE3);
    if (ATIC2_proton) add_to_mg(atic2_proton_path, kATIC2_proton, atic2_marker_style, 1.3, atic2_proton_legendentry, mgE3);
    if (DAMPE_proton) add_to_mg(dampe_proton_path, kDAMPE_proton, dampe_marker_style, 1, dampe_proton_legendentry, mgE3);
    if (AMS02_electron) add_to_mg(ams02_electron_path, kAMS02_electron, ams02_marker_style, 1, ams02_electron_legendentry, mgE3);
    if (AMS02_positron) add_to_mg(ams02_positron_path, kAMS02_positron, ams02_marker_style, 1, ams02_positron_legendentry, mgE3);
    if (AMS02_pbar) add_to_mg(ams02_pbar_path, kAMS02_pbar, ams02_marker_style, 1, ams02_pbar_legendentry, mgE3, ams02_pbar_filter);
    if (AMS02_B) add_to_mg(ams02_B_path, kAMS02_B, ams02_marker_style, 1, ams02_B_legendentry, mgE3);
    if (AMS02_C) add_to_mg(ams02_C_path, kAMS02_C, ams02_marker_style, 1, ams02_C_legendentry, mgE3);
    if (AMS02_O) add_to_mg(ams02_O_path, kAMS02_O, ams02_marker_style, 1, ams02_O_legendentry, mgE3);
    if (CALET_electron) add_to_mg(calet_electron_path, kCALET_electron, calet_marker_style, 1.3, calet_allele_legendentry, mgE3);
    if (FERMILAT_electron) add_to_mg(fermilat_electron_path, kFERMILAT_electron, fermi_marker_style, 1, fermilat_allele_legendentry, mgE3);
    if (DAMPE_electron) add_to_mg(dampe_electron_path, kDAMPE_electron, dampe_marker_style, 1, dampe_allele_legendentry, mgE3);
    if (AUGER_allparticle) add_to_mg(auger_allparticle_path, kAUGER_allparticle, auger_marker_style, 1.3, auger_allparticle_legendentry, mgE3, auger_filter);
    if (KASCADEgrande_allparticle) add_to_mg(kascadegrande_allparticle_path, kKASCADEgrande_allparticle, kascadegrande_marker_style, 1.3, kascadegrande_allparticle_legendentry, mgE3, kascadegrande_filter);

    mgE3->GetXaxis()->SetTitle("Kinetic Energy/n (GeV/n)");
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
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
        {
            auto tmp_obj = static_cast<TGraphAsymmErrors *>(key->ReadObj());
            if (exclude_obj(tmp_obj))
                continue;
            tmp_obj->SetLineColor(mycolor->GetNumber());
            tmp_obj->SetMarkerColor(mycolor->GetNumber());
            if (std::string(tmp_obj->GetName()).find("upper_limit") != std::string::npos) marker_style = 23;
            tmp_obj->SetMarkerStyle(marker_style);
            tmp_obj->SetMarkerSize(marker_size);
            tmp_obj->SetTitle(gr_title);
            
            // Set error to zero along X axis
            for (int i = 0; i < tmp_obj->GetN(); ++i)
            {
                tmp_obj->SetPointEXlow(i, 0);
                tmp_obj->SetPointEXhigh(i, 0);
            }

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
    cflux.SetLogy();
    cflux.SetTicks();
    gPad->Modified(); 
    gPad->Update();
    mgr->GetXaxis()->SetLabelSize(0.045);
    mgr->GetXaxis()->SetTitleSize(0.045);
    mgr->GetXaxis()->SetTitleOffset(1.3);
    mgr->GetYaxis()->SetLabelSize(0.045);
    mgr->GetYaxis()->SetTitleSize(0.045);
    mgr->GetYaxis()->SetTitleOffset(1.3);
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