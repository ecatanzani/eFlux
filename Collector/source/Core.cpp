#include "myHeader.h"
#include "flux.h"
#include "produceTuples.h"
#include "data_loop.h"

void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const bool rawdata_flag,
	const bool skimmed_flag,
	const bool ntuple_flag,
    AnyOption &opt,
    const std::string wd)
{   
    if (ntuple_flag)  
        produceTuples(
            opt,
            inputPath,
            verbose,
            wd);
    else
    {
        // Create output TFile
        const char* outFilePath = static_cast<const char*>(uniqueOutFile(outputPath, opt).c_str());
        if (verbose)
            std::cout << "\nCreating output ROOT file... [" << outFilePath << "]" << std::endl;
        TFile outFile(outFilePath, "NEW", "Analysis Output File");
        if (!outFile.IsOpen())
        {
            std::cerr << "\n\nError (100) writing output TFile... [" << outFilePath << "]" << std::endl;
            exit(100);
        }

        if (rawdata_flag)
		    rawDataLoop(
                inputPath,
			    outFile,
			    verbose,
			    skimmed_flag,
			    wd);
        if (skimmed_flag)
            skimmedDataLoop(
                inputPath,
			    outFile,
			    verbose,
			    skimmed_flag,
			    wd);

        // Close output file ...
        outFile.Close();
    }
}