#include "myHeader.h"

void eCore(
            const std::string inputPath,
            const std::string outputPath,
            const bool verbose,
            const unsigned int lvTime,
            AnyOption &opt
        )
{
    //Create output TFile
    TFile outFile(uniqueOutFile(outputPath,opt),"NEW","Analysis Output File");
    if(!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << uniqueOutFile(outputPath,opt) << std::endl;
        exit(123);
    }
    
    //Building eFLux
    buildFlux(inputPath,lvTime,outFile,verbose);

    //Close output file ...
    outFile.Close();

}