#include "main.h"
#include "utils.h"

int main(int argc, char **argv) {

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help         .......... Prints this help");
	opt.addUsage(" -i  --input        .......... <path_to_input_list>              .......... (*) Input list");
	opt.addUsage(" -o  --output       .......... <path_to_output_TFile>            .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputdir    .......... <path_to_output_TFile_dir>        .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose      .......... Verbose output");

	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc           .......... MC event loop");
	
	opt.addUsage("");
    opt.addUsage(" -p  --parallel     .......... <number_of_threads>               .......... Multithreading option");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');
    opt.setOption("parallel", 'p');
	
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
	if (opt.getValue("output") || opt.getValue('o'))
		input_pars.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_pars.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_pars.verbose = opt.getFlag('v');
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_pars.mc = true;
    if (opt.getValue("parallel") || opt.getValue('p'))
		input_pars.threads =  std::stoul(opt.getValue('p'), nullptr, 0);
	
	
	
	
	return 0;
}