#include "MC/mc.h"
#include "Data/data.h"
#include "main.h"
#include "utils.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help         .......... Prints this help");
	opt.addUsage(" -i  --input        .......... <path_to_input_list>              .......... (*) Input list");
	opt.addUsage(" -w  --workdir      .......... <path_to_software_config_dir>     .......... (*) Config directory");
	opt.addUsage(" -o  --output       .......... <path_to_output_TFile>            .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputDir    .......... <path_to_output_TFile_dir>        .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose      .......... Verbose output");
	opt.addUsage(" -p  --pedantic     .......... Pedantic output");

	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc           .......... MC event loop");
	opt.addUsage(" -r  --raw_data     .......... Raw data event loop");
	
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

	// Load args struct
	in_pars input_pars;
	
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
		input_pars.input_path = opt.getValue('i');
	if (opt.getValue("workdir") || opt.getValue('w'))
		input_pars.wd = opt.getValue('w');
	input_pars.output_path = uniqueOutFile(opt);
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_pars.verbose = opt.getFlag('v');
	if (opt.getFlag("pedantic") || opt.getFlag('p'))
		input_pars.pedantic = opt.getFlag('p');
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_pars.mc_flag = true;
	if (opt.getFlag("raw_data") || opt.getFlag('r'))
		input_pars.rawdata_flag = true;
	
	if (input_pars.CheckArgs())
	{
		if (input_pars.verbose)
			input_pars.ExpandArgs();
		if (input_pars.mc_flag)
			mcCore(input_pars);
		else
			dataCore(input_pars);
	}	
	else
		return 100;
	
	return 0;
}