#include "acceptance.h"
#include "myHeader.h"

void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    AnyOption &opt,
    const std::string wd)
{
    // Create output TFile
    const char* outFilePath = static_cast<const char*>(uniqueOutFile(outputPath, opt).c_str());
    TFile outFile(outFilePath, "NEW", "Analysis Output File");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath << std::endl;
        exit(123);
    }

    // Create energy log-binning
    energy_cuts eCuts;
    load_energy_struct(eCuts, wd);
    auto logEBins = createLogBinning(eCuts);
    if (pedantic)
    {
        std::cout << "\nEnergy log binning..." << std::scientific;
        for (auto it = logEBins.begin(); it != logEBins.end(); ++it)
            std::cout << "\n"
                      << *it;
        std::cout << std::defaultfloat;
    }

    // Build acceptance
    buildAcceptance(
        accInputPath,
        verbose,
        logEBins,
        outFile,
        wd);
    
    // Close output file ...
    outFile.Close();
}