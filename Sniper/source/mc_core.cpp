#include "mc.h"
#include "utils.h"

void mcCore(
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
	TFile outFile(outFilePath, "NEW", "Analysis Output File");
	if (!outFile.IsOpen())
	{
		std::cerr << "\n\nError writing output TFile: " << outFilePath << std::endl;
		exit(123);
	}
	// Load event display filter logger
	std::shared_ptr<ofstream> evlogger (nullptr);
	if (edfilter_flag)
		evlogger = ev_logger(outputPath, opt);

	mcLoop(
		inputPath,
		outFile,
		verbose,
		evlogger,
		wd);

	// Close output file ...
	outFile.Close();
}