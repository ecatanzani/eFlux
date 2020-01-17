#include "myHeader.h"

void eCore(
            const std::string inputPath,
            const std::string outputPath,
            const bool verbose,
            const unsigned int lvTime
            )
{
    //Create output TFile
    TFile outFile(uniqueOutFile(outputPath),"NEW","Analysis Output File");
    if(!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << uniqueOutFile(outputPath) << std::endl;
        exit(123);
    }
    
    //Building eFLux
    buildFlux(inputPath,lvTime);

    //Close output file ...
    outFile.Close();

}