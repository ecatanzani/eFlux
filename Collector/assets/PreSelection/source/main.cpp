#include "main.h"
#include "utils.h"
#include "preselection.h"

int main(int argc, char **argv) {

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help         .......... Prints this help");
	opt.addUsage(" -i  --input        .......... <path_to_input_list>              .......... (*) Input list");
	opt.addUsage(" -w  --workdir      .......... <path_to_software_config_dir>     .......... (*) Local config directory");
	opt.addUsage(" -o  --output       .......... <path_to_output_TFile>            .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputdir    .......... <path_to_output_TFile_dir>        .......... Output ROOT TFile directory");
	opt.addUsage(" -l  --logs         .......... <path_to_output_logs_dir>         .......... Output logs directory");
	opt.addUsage(" -v  --verbose      .......... Verbose output");

	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc           .......... MC event loop");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setOption("logs", 'l');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');
	
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
	if (opt.getValue("output") || opt.getValue('o'))
		input_pars.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_pars.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getValue("logs") || opt.getValue('l'))
		input_pars.logs_dir = opt.getValue('l');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_pars.verbose = opt.getFlag('v');
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_pars.mc = true;
	
	preselection(input_pars);
	
	return 0;
}