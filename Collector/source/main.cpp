#include "main.h"
#include "mc.h"
#include "data.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -i  --input			<path_to_input_list>				(*)	Input list");
	opt.addUsage(" -w  --workdir		<path_to_software_config_dir>		(*) Config directory");
	opt.addUsage(" -o  --output			<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -d  --outputDir		<path_to_output_TFile_dir>				Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose												Verbose output");
	opt.addUsage(" -p  --pedantic												Pedantic output");

	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc														MC event loop");
	opt.addUsage(" -r  --raw_data												Raw data event loop");
	
	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("pedantic", 'p');
	opt.setFlag("mc", 'm');
	opt.setFlag("raw_data", 'r');
	
	opt.processCommandArgs(argc, argv);

	std::string wd;
	std::string inputPath;
	std::string outputPath;

	bool verbose = false;
	bool pedantic = false;

	// Tasks
	bool mc_flag = false;
	bool rawdata_flag = false;

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
	if (opt.getFlag("pedantic") || opt.getFlag('p'))
		pedantic = opt.getFlag('p');

	if (opt.getFlag("mc") || opt.getFlag('m'))
		mc_flag = true;
	if (opt.getFlag("raw_data") || opt.getFlag('r'))
		rawdata_flag = true;

	if (mc_flag)
		mcCore(
			inputPath,
			outputPath,
			verbose,
			pedantic,
			opt,
			wd);
	else if (rawdata_flag)
		dataCore(
			inputPath,
			outputPath,
			verbose,
			pedantic,
			opt,
			wd);
	else
		std::cerr << "\n\nERROR: Wrong Task selected \n\n";
	
	return 0;
}