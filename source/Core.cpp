#include "myHeader.h"

void eCore(
            const std::string inputPath,
            const std::string outputPath,
            const bool verbose
            )
{
    //Create output TFile
    TFile outFile(uniqueTree(outputPath),"NEW","Analysis Output TTree");
    if(!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << uniqueTree(outputPath) << std::endl;
        exit(123);
    }
    
    //Building eFLux
    buildFlux(outFile);

    //Close output file ...
    outFile.Close();

}