#include "acceptance.h"
#include "binning.h"
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
    const char *outFilePath = static_cast<const char *>(uniqueOutFile(outputPath, opt).c_str());
    TFile outFile(outFilePath, "NEW", "Analysis Output File");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath << std::endl;
        exit(123);
    }
    
    // Build acceptance
    buildAcceptance(
        accInputPath,
        verbose,
        outFile,
        wd);

    // Close output file ...
    outFile.Close();
}