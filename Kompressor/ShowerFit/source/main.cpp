#include "main.h"
#include "utils.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -i  --input			<path_to_input_list>				(*)	Input file list");
	opt.addUsage(" -w  --workdir		<path_to_software_config_dir>		(*) Config directory");
	opt.addUsage(" -o  --output			<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -d  --outputdir		<path_to_output_TFile_dir>				Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose												Verbose output");
	
	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc														MC event loop");

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');

	opt.processCommandArgs(argc, argv);
	
	std::string wd;
	std::string input_list;
	std::string output_path;

	bool _VERBOSE = false;

	// Tasks
	bool mc_flag = false;

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
		input_list = opt.getValue('i');
	if (opt.getValue("workdir") || opt.getValue('w'))
		wd = opt.getValue('w');
	if (opt.getValue("output") || opt.getValue('o'))
		output_path = opt.getValue('o');
	if (opt.getValue("outputdir") || opt.getValue('d'))
		output_path = opt.getValue('d');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		_VERBOSE = opt.getFlag('v');
	
	if (opt.getFlag("mc") || opt.getFlag('m'))
		mc_flag = true;

	reader(
		wd, 
		input_list, 
		expand_output_path(opt, output_path), 
		_VERBOSE, 
		mc_flag);
		
	return 0;
}