#include "main.h"
#include "preselection.h"

int main(int argc, char **argv) {

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help         .......... Prints this help");
	opt.addUsage(" -i  --input        .......... <path_to_input_list>              .......... (*) Input list");
	opt.addUsage(" -w  --workdir      .......... <path_to_software_config_dir>     .......... (*) Energy Common config directory");
	opt.addUsage(" -d  --outputdir    .......... <path_to_output_TFile_dir>        .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose      .......... Verbose output");

	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc           .......... MC event loop");
	opt.addUsage(" -r  --raw_data     .......... Raw data event loop");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');
	opt.setFlag("raw_data", 'r');
	
	opt.processCommandArgs(argc, argv);

	// Load args struct
	in_pars input_pars;
	
	if (!opt.hasOptions()) {
		opt.printUsage();
		return 0;
	}
	if (opt.getFlag("help") || opt.getFlag('h')) {
		opt.printUsage();
		return 0;
	}
	if (opt.getValue("input") || opt.getValue('i'))
		input_pars.input_path = opt.getValue('i');
	if (opt.getValue("workdir") || opt.getValue('w'))
		input_pars.config_wd = opt.getValue('w');
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_pars.output_wd = opt.getValue('d');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_pars.verbose = opt.getFlag('v');
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_pars.mc_flag = true;
	if (opt.getFlag("raw_data") || opt.getFlag('r'))
		input_pars.rawdata_flag = true;
	
	preselection(input_pars);
	
	return 0;
}