#include "histos.h"

#include <string>
#include <iostream>

#include "TFile.h"

std::vector<std::vector<std::shared_ptr<TH1D>>> getrmslayerhistos()
{
    std::vector<std::vector<std::shared_ptr<TH1D>>> rmslayer (nenergybin);
    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        rmslayer[bidx] = std::vector<std::shared_ptr<TH1D>> (bgolayers);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            std::string tmphistoname = std::string("h_rms_energybin_") + std::to_string(bidx+1) + std::string("_layer_") + std::to_string(lidx+1);
            std::string tmphistotitle = std::string("RMS - energybin ") + std::to_string(bidx+1) + std::string(" - BGO layer ") + std::to_string(lidx+1);
            rmslayer[bidx][lidx] = std::make_shared<TH1D>(tmphistoname.c_str(), tmphistotitle.c_str(), 100, -5, 5);
        }
    }
    return rmslayer;
}

std::vector<std::vector<std::shared_ptr<TH1D>>> getrmslayerhistos_ff(const std::string inputfile)
{
    std::vector<std::vector<std::shared_ptr<TH1D>>> rmslayer (nenergybin);
    TFile *infile = TFile::Open(inputfile.c_str(), "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError opening input ROOT file: [" << inputfile << "]\n\n";
        exit(100);
    }
    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        rmslayer[bidx] = std::vector<std::shared_ptr<TH1D>> (bgolayers);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            std::string tmphistoname = std::string("RMS/h_rms_energybin_") + std::to_string(bidx+1) + std::string("_layer_") + std::to_string(lidx+1);
            // Read histo from file
            rmslayer[bidx][lidx] = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get(tmphistoname.c_str())));
            // Fit histo with gaus normalized distribution
            rmslayer[bidx][lidx]->Fit("gaus", "QW");
            // Set histo ownership
            rmslayer[bidx][lidx]->SetDirectory(0);
        }
    }
    infile->Close();
    return rmslayer;
}

std::vector<std::vector<std::shared_ptr<TH1D>>> getenergyfractionlayerhistos()
{
    std::vector<std::vector<std::shared_ptr<TH1D>>> energyfractionlayer (nenergybin);
    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        energyfractionlayer[bidx] = std::vector<std::shared_ptr<TH1D>> (bgolayers);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            std::string tmphistoname = std::string("h_energyfraction_energybin_") + std::to_string(bidx+1) + std::string("_layer_") + std::to_string(lidx+1);
            std::string tmphistotitle = std::string("Energy Fraction - energybin ") + std::to_string(bidx+1) + std::string(" - BGO layer ") + std::to_string(lidx+1);
            energyfractionlayer[bidx][lidx] = std::make_shared<TH1D>(tmphistoname.c_str(), tmphistotitle.c_str(), 100, -5, 5);
        }
    }
    return energyfractionlayer;
}

std::vector<std::vector<std::shared_ptr<TH1D>>> getenergyfractionlayerhistos_ff(const std::string inputfile)
{
    std::vector<std::vector<std::shared_ptr<TH1D>>> energyfractionlayer (nenergybin);
    TFile *infile = TFile::Open(inputfile.c_str(), "READ");
    if (infile->IsZombie())
    {
        std::cerr << "\n\nError opening input ROOT file: [" << inputfile << "]\n\n";
        exit(100);
    }
    for (unsigned int bidx=0; bidx<nenergybin; ++bidx)
    {
        energyfractionlayer[bidx] = std::vector<std::shared_ptr<TH1D>> (bgolayers);
        for (unsigned int lidx=0; lidx<bgolayers; ++lidx)
        {
            std::string tmphistoname = std::string("ELF/h_energyfraction_energybin_") + std::to_string(bidx+1) + std::string("_layer_") + std::to_string(lidx+1);
            // Read histo from file
            energyfractionlayer[bidx][lidx] = std::shared_ptr<TH1D>(static_cast<TH1D*>(infile->Get(tmphistoname.c_str())));
            // Fit histo with gaus normalized distribution
            energyfractionlayer[bidx][lidx]->Fit("gaus", "QW");
            // Set histo ownership
            energyfractionlayer[bidx][lidx]->SetDirectory(0);
        }
    }
    infile->Close();
    return energyfractionlayer;
}