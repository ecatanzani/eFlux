#include <iostream>

#include "TF1.h"
#include "TFile.h"

void TF1Ratio(
    const char *data_file,
    const char *mc_file,
    const char *tf1_name = "fitfunc",
    const double emin_gev=1000,
    const double emax_gev=2000,
    const char *outputfile="ratiofitresult.root")
{
    TFile *dfile = TFile::Open(data_file, "READ");
    if (dfile->IsZombie())
    {
        std::cerr << "\n\nError opening input data file [" << data_file << "\n\n";
        exit(100);
    }

    auto data_f = static_cast<TF1*>(dfile->Get(tf1_name));
    data_f->SetName("data_f");
    data_f->SetTitle("data_f");
    //data_f->SetDirectory(0);
    //dfile->Close();
    
    TFile *mfile = TFile::Open(mc_file, "READ");
    if (mfile->IsZombie())
    {
        std::cerr << "\n\nError opening input mc file [" << mc_file << "\n\n";
        exit(100);
    }

    auto mc_f = static_cast<TF1*>(mfile->Get(tf1_name));
    mc_f->SetName("mc_f");
    mc_f->SetTitle("mc_f");
    //mc_f->SetDirectory(0);
    //mfile->Close();

    TF1 ratio("ratio", "data_f*TMath::Power(mc_f, -1)", emin_gev, emax_gev);

    TFile *ofile = TFile::Open(outputfile, "RECREATE");
    if (ofile->IsZombie())
    {
        std::cerr << "\n\nError writing output file [" << outputfile << "\n\n";
        exit(100);
    }

    data_f->Write();
    mc_f->Write();
    ratio.Write();

    ofile->Close();
}