#include "myHeader.h"
#include "flux.h"

void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const unsigned int lvTime,
    const std::string accInputPath,
    AnyOption &opt)
{
    // Create output TFile
    TFile outFile(uniqueOutFile(outputPath, opt), "NEW", "Analysis Output File");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << uniqueOutFile(outputPath, opt) << std::endl;
        exit(123);
    }
    
    // Create energy log-binning
    energy_cuts eCuts;
    load_energy_struct(eCuts);
    auto logEBins = createLogBinning(eCuts);
    if (pedantic)
    {
        std::cout << "\nEnergy log binning..." << std::scientific;
        for (auto it = logEBins.begin(); it != logEBins.end(); ++it)
            std::cout << "\n" << *it;
        std::cout << std::defaultfloat;
    }

    // All-Electron flux using xtrl as classifier
    buildXtrlFlux(
        logEBins,
        inputPath,
        lvTime,
        outFile,
        verbose,
        accInputPath);

    // Close output file ...
    outFile.Close();
}