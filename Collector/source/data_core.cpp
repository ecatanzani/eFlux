#include "data.h"
#include "utils.h"
#include "tuple.h"

void dataCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool verbose,
	const bool pedantic,
	const bool rawdata_flag,
	const bool ntuple_flag,
	AnyOption &opt,
	const std::string wd)
{
	// Create output TFile
	const char *outFilePath = static_cast<const char *>(uniqueOutFile(outputPath, opt).c_str());
	if (verbose)
		std::cout << "\nCreating output ROOT file... [" << outFilePath << "]" << std::endl;
	TFile outFile(outFilePath, "NEW", "Analysis Output File");
	if (!outFile.IsOpen())
	{
		std::cerr << "\n\nError (100) writing output TFile... [" << outFilePath << "]" << std::endl;
		exit(100);
	}

	// Ntuple facility
	if (ntuple_flag)
		tupleLoop(
			inputPath,
			outFile,
			verbose,
			wd);

	// Raw data loop
	if (rawdata_flag)
		rawDataLoop(
			inputPath,
			outFile,
			verbose,
			wd);

	// Close output file ...
	outFile.Close();
}