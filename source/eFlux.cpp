#include "myHeader.h"
#include "acceptance.h"
#include "flux.h"
#include "data_loop.h"

void buildFlux(
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath,
    const std::string wd)
{
    /*
    auto acceptance_tf1 = readAcceptance(
        outFile,
        verbose,
        accInputPath);
    */
   
    auto eCounts =
        evLoop(
            inputPath,
            outFile,
            verbose,
            wd);
}