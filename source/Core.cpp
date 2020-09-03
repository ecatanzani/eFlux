#include "myHeader.h"
#include "flux.h"
#include "binning.h"

void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const unsigned int lvTime,
    const std::string accInputPath,
    const std::string p_accInputPath,
    AnyOption &opt,
    const std::string wd)
{   
    if (opt.getValue("ntuple") || opt.getValue('n'))
        produceTuples(
            opt,
            inputPath,
            verbose,
            wd);
    else
    {
        /*
        // Extract electron acceptance
        if (verbose)
            std::cout << "\n\nReading electron acceptance geometry file... [" << accInputPath << "]";
        TFile* electron_mc_file = TFile::Open(accInputPath.c_str(), "READ");
        if (electron_mc_file->IsZombie())
        {
            std::cerr << "\n\nError reading input electron MC reco file: " << accInputPath << std::endl;
            exit(100);
        }
        TH1D* electron_acceptance = static_cast<TH1D*>(electron_mc_file->Get("Acceptance_histos/h_acceptance_all_cut"));
        electron_acceptance->SetDirectory(0);
        electron_mc_file->Close();

        // Extract proton contamination
        if (verbose)
            std::cout << "\nReading proton acceptance geometry file... [" << p_accInputPath << "]" << std::endl;
        TFile* proton_mc_file = TFile::Open(p_accInputPath.c_str(), "READ");
        if (proton_mc_file->IsZombie())
        {
            std::cerr << "\n\nError reading input proton MC reco file: " << p_accInputPath << std::endl;
            exit(100);
        }

        TH1D* proton_background_fraction = static_cast<TH1D*>(proton_mc_file->Get("mc_ancillary/proton_background_ratio"));
        proton_background_fraction->SetDirectory(0);
        proton_mc_file->Close();
        */

        // Create output TFile
        const char* outFilePath = static_cast<const char*>(uniqueOutFile(outputPath, opt).c_str());
        if (verbose)
            std::cout << "\nCreating output ROOT file... [" << outFilePath << "]" << std::endl;
        TFile outFile(outFilePath, "NEW", "Analysis Output File");
        if (!outFile.IsOpen())
        {
            std::cerr << "\n\nError writing output TFile: " << outFilePath << std::endl;
            exit(123);
        }
        /*

        buildFlux(
            inputPath,
            lvTime,
            outFile,
            verbose,
            accInputPath,
            electron_acceptance,
            proton_background_fraction,
            wd);
        */

        buildFlux(
            inputPath,
            lvTime,
            outFile,
            verbose,
            accInputPath,
            wd);

        // Close output file ...
        outFile.Close();
    }
}