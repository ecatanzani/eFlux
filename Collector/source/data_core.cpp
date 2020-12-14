#include "data.h"
#include "utils.h"
#include "tuple.h"

void dataCore(
	const std::string inputPath,
	const std::string outputPath,
	const bool verbose,
	const bool pedantic,
	const bool edfilter_flag,
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
	// Load event display filter logger
	std::shared_ptr<ofstream> evlogger (nullptr);
	if (edfilter_flag)
		evlogger = ev_logger(outputPath, opt);

	rawDataLoop(
		inputPath,
		outFile,
		verbose,
		evlogger,
		wd);

	// Close output file ...
	outFile.Close();
}