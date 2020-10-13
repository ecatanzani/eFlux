#include "main.h"
#include "ntuple.h"
#include "environment.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -i  --input			<path_to_input_list>				(*)	Input nTuple list");
	opt.addUsage(" -w  --workdir		<path_to_software_config_dir>		(*) Config directory");
	opt.addUsage(" -o  --output			<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -d  --outputDir		<path_to_output_TFile_dir>				Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose												Verbose output");
	
	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');

	opt.processCommandArgs(argc, argv);
	
	std::string wd;
	std::string inputPath;
	std::string outputPath;

	bool verbose = false;

	if (!opt.hasOptions())
	{
		opt.printUsage();
		return 0;
	}
	if (opt.getFlag("help") || opt.getFlag('h'))
	{
		opt.printUsage();
		return 0;
	}
	if (opt.getValue("input") || opt.getValue('i'))
		inputPath = opt.getValue('i');
	if (opt.getValue("workdir") || opt.getValue('w'))
		wd = opt.getValue('w');
	if (opt.getValue("output") || opt.getValue('o'))
		outputPath = opt.getValue('o');
	if (opt.getValue("outputDir") || opt.getValue('d'))
		outputPath = opt.getValue('d');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		verbose = opt.getFlag('v');

#if _PARALLEL
	read_ntuple(
		opt, 
		inputPath, 
		outputPath,
		wd,
		verbose);
#else
	read_ntuple(
		opt, 
		inputPath, 
		outputPath,
		wd,
		verbose,
		false);
#endif
	
	return 0;
}