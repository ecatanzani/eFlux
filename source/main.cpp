#include "myHeader.h"
#include "acceptance.h"

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
	
	opt.addUsage(" -a  --acceptance												Acceptance calculation");
	opt.addUsage(" -r  --raw_data												Raw data event loop");
	opt.addUsage(" -s  --skimmed												Skimmed data event loop");
	opt.addUsage(" -n  --ntuple													nTuple production facility");
	opt.addUsage(" -c  --collect		<mc/data>								DATA/MC reconstruction facility");
	
	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("pedantic", 'p');
	opt.setFlag("acceptance", 'a');
	opt.setFlag("raw_data", 'r');
	opt.setFlag("skimmed", 's');
	opt.setFlag("ntuple", 'n');
	opt.setOption("collect", 'c');

	opt.processCommandArgs(argc, argv);

	std::string wd;
	std::string inputPath;
	std::string outputPath;

	bool verbose = false;
	bool pedantic = false;

	// Tasks
	bool acceptance_flag = false;
	bool rawdata_flag = false;
	bool skimmed_flag = false;
	bool ntuple_flag = false;
	bool collect_flag = false;
	bool collect_data = false;
	bool collect_mc = false;

	if (!opt.hasOptions())
		opt.printUsage();
	if (opt.getFlag("help") || opt.getFlag('h'))
		opt.printUsage();
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

	if (opt.getFlag("acceptance") || opt.getFlag('a'))
		acceptance_flag = true;
	if (opt.getFlag("raw_data") || opt.getFlag('r'))
		rawdata_flag = true;
	if (opt.getFlag("skimmed") || opt.getFlag('s'))
		skimmed_flag = true;
	if (opt.getFlag("ntuple") || opt.getFlag('n'))
		ntuple_flag = true;
	if (opt.getValue("collect") || opt.getValue('c'))
	{
		if (!strcmp(opt.getValue('r'), "mc"))
			collect_mc = true;
		else if (!strcmp(opt.getValue('r'), "data"))
			collect_data = true;
		else
		{
			std::cerr << "\n\nWrong --resume config...";
			opt.printUsage();
			exit(100);
		}
	}

	if (acceptance_flag)
		computeAcceptance(
			inputPath,
			verbose,
			pedantic,
			outputPath,
			opt,
			wd);
	else if(rawdata_flag || skimmed_flag || ntuple_flag)
		eCore(
			inputPath,
			outputPath,
			verbose,
			pedantic,
			rawdata_flag,
			skimmed_flag,
			ntuple_flag,
			opt,
			wd);
	else if (collect_flag)
	{
		if (collect_mc)
			generateFinalGraph(
				verbose,
				pedantic,
				outputPath,
				inputPath,
				wd);
		if (collect_data)
			generateDataFinalGraph(
				verbose,
				pedantic,
				outputPath,
				inputPath,
				wd);
	}

	return 0;
}