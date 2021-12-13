#include "main.h"
#include "utils.h"
#include "bdtinfo.h"

int main(int argc, char **argv)
{
	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                        .......... Prints this help");
	opt.addUsage(" -i  --input                       .......... <path_to_input_list>                             .......... (*) Input file list");
	opt.addUsage(" -c  --config                      .......... <path_to_config_dir>                             .......... Collector working directory");
	opt.addUsage(" -o  --output                      .......... <path_to_output_TFile>                           .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputDir                   .......... <path_to_output_TFile_dir>                       .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose                     .......... Verbose output");

	opt.addUsage("");
	opt.addUsage("TMVA learning method: ");
	opt.addUsage("");

	opt.addUsage(" -m  --method                      .......... <TMVA_learning_method>                           .......... TMVA learning/classifying method");

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("config", 'c');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setOption("method", 'm');
	
	opt.processCommandArgs(argc, argv);
	
	// Initialize args struct
	in_args input_args;
	
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
		input_args.input_list = opt.getValue('i');
	if (opt.getValue("config") || opt.getValue('c'))
		input_args.config_dir = opt.getValue('c');
	if (opt.getValue("workdir") || opt.getValue('w'))
		input_args.config_dir = opt.getValue('w');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_args.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getValue("method") || opt.getValue('m'))
		input_args.learning_method = opt.getValue('m');
	
	if (!input_args.learning_method.empty()) ExtractBDTInfo(input_args);
	else {
		std::cerr << "\nError ! No learnign method has been specified...\n\n";
		exit(100);
	}
	
	return 0;
}