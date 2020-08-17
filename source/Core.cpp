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
        // Create output TFile
        const char* outFilePath = static_cast<const char*>(uniqueOutFile(outputPath, opt).c_str());
        TFile outFile(outFilePath, "NEW", "Analysis Output File");
        if (!outFile.IsOpen())
        {
            std::cerr << "\n\nError writing output TFile: " << outFilePath << std::endl;
            exit(123);
        }

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