#include <memory>
#include <string>
#include <iostream>

#include "TROOT.h"
#include "TRint.h"
#include "TAxis.h"
#include "TKey.h"
#include "TFile.h"
#include "TColor.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TApplication.h"
#include "TGraphAsymmErrors.h"

const bool AMS02_proton = true;
const bool CALET_proton = true;
const bool DAMPE_proton = true;

TColor *kAMS02 = gROOT->GetColor(kBlue);
TColor *kDAMPE = gROOT->GetColor(kRed);
TColor *kCALET = gROOT->GetColor(kOrange);

const char* ams02_proton_path = "data/ssdc_AMS02_p_E3.root";
const char* calet_proton_path = "data/ssdc_CALET_p_E3.root";
const char* dampe_proton_path = "data/ssdc_DAMPE_p_E3.root";

const double emin = 1;
const double emax = 1e+5;
const double phi_min = 1;
const double phi_max = 1e+6;

const char* multigraph_name = "flux_mg_E3";

void add_to_mg(const char* path, TColor* mycolor, std::shared_ptr<TMultiGraph> mgE3);
void buildflux(const char* outfile);
void drawflux(const char* infile);

void buildflux(const char* outfile)
{
    std::shared_ptr<TMultiGraph> mgE3 = std::make_shared<TMultiGraph>();
    mgE3->SetName(multigraph_name);

    if (AMS02_proton)
        add_to_mg(ams02_proton_path, kAMS02, mgE3);
    if (CALET_proton)
        add_to_mg(calet_proton_path, kCALET, mgE3);
    if (DAMPE_proton)
        add_to_mg(dampe_proton_path, kDAMPE, mgE3);


    mgE3->GetXaxis()->SetTitle("Energy (GeV)");
    mgE3->GetYaxis()->SetTitle("E^{3} #times #Phi (GeV [m^{2} sr s]^{-1})");
    
    TFile* myoutfile = TFile::Open(outfile, "RECREATE");
    if (myoutfile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outfile << "]\n";
        exit(100);
    }

    mgE3->Write();

    myoutfile->Close();
}

void add_to_mg(const char* path, TColor* mycolor, std::shared_ptr<TMultiGraph> mgE3)
{
    TFile* inFile = TFile::Open(path, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file [" << path << "]\n";
        exit(100);
    }

    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
        {
            auto tmp_obj = static_cast<TGraphAsymmErrors *>(key->ReadObj());
            tmp_obj->SetLineColor(mycolor->GetNumber());
            tmp_obj->SetMarkerColor(mycolor->GetNumber());
            mgE3->Add(tmp_obj);
        }
    }
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