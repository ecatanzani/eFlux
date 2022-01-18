#include "main.h"
#include "split.h"

int main(int argc, char **argv)
{
	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                        .......... Prints this help");
	opt.addUsage(" -i  --input                       .......... <path_to_input_list>                             .......... (*) Input file list");
	opt.addUsage(" -c  --config                      .......... <path_to_energy_config_file>                     .......... Path to energy config file");
	opt.addUsage(" -o  --output                      .......... <path_to_output_TFile_dir>                       .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose                     .......... Verbose output");

	opt.addUsage("");
	opt.addUsage(" -p  --parallel                    .......... <number_of_threads>                              .......... Multithreading option");

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("config", 'c');
	opt.setOption("output", 'o');
	opt.setFlag("verbose", 'v');
	opt.setOption("parallel", 'p');
	
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
		input_args.energy_config_file = opt.getValue('c');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_directory = opt.getValue('o');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getValue("parallel") || opt.getValue('p'))
		input_args.threads =  std::stoul(opt.getValue('p'), nullptr, 0);

	Split(input_args);
	
	return 0;
}