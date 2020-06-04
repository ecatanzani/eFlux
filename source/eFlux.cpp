#include "myHeader.h"
#include "acceptance.h"
#include "flux.h"
#include "dataLoop.h"

void buildFlux(
    const std::vector<float> &logEBins,
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath)
{
    auto acceptance_tf1 = readAcceptance(
        outFile,
        verbose,
        accInputPath);

    auto eCounts =
        evLoop(
            logEBins,
            inputPath,
            outFile,
            verbose);
    
    
}