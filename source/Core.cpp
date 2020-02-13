#include "myHeader.h"

void eCore(
            const std::string inputPath,
            const std::string outputPath,
            const bool verbose,
            const unsigned int lvTime,
            const bool myAcceptance,
            const std::string accInputPath,
            AnyOption &opt
        )
{
    // Create output TFile
    TFile outFile(uniqueOutFile(outputPath,opt),"NEW","Analysis Output File");
    if(!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << uniqueOutFile(outputPath,opt) << std::endl;
        exit(123);
    }
    
    // Building acceptance
    if(myAcceptance)
        buildAcceptance(accInputPath,verbose);

    // Building eFLux
    buildFlux(inputPath,lvTime,outFile,verbose);

    // Close output file ...
    outFile.Close();

}