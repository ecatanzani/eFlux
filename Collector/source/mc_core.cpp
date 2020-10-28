#include "mc.h"
#include "utils.h"

void mcCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool verbose,
	const bool pedantic,
	const bool ntuple_flag,
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

	mcLoop(
		inputPath,
		outFile,
		verbose,
		ntuple_flag,
		wd);

	// Close output file ...
	outFile.Close();
}