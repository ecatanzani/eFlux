#include "myHeader.h"
#include "acceptance.h"
#include "flux.h"
#include "dataLoop.h"

void buildFlux(
    const std::vector<float> &eBins,
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
    
    auto e_dataCounts =
        evLoop(
            eBins,
            inputPath,
            outFile,
            verbose);
    
    
}