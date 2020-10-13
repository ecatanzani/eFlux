#include "ntuple.h"
#include "unique.h"
#include "read_parallel.h"
#include "read_serial.h"

void read_ntuple(
	AnyOption &opt,
	const std::string inputPath, 
	const std::string outputPath,
	const std::string wd,
	const bool verbose,
	const bool mthreads)
{
	const char* outFilePath = static_cast<const char*>(uniqueOutFile(outputPath, opt).c_str());
	if (verbose)
        std::cout << "\nCreating output ROOT file... [" << outFilePath << "]" << std::endl;
	TFile outFile(outFilePath, "NEW", "NTuple Analysis Output File");
	if (!outFile.IsOpen())
	{
		std::cerr << "\n\nError (100) writing output TFile... [" << outFilePath << "]" << std::endl;
		exit(100);
	}

	if (mthreads)
		read_parallel(
			inputPath,
			outFile,
			wd,
			verbose);
	else
		read_serial(
			inputPath,
			outFile,
			wd,
			verbose);
	
	outFile.Close();
}