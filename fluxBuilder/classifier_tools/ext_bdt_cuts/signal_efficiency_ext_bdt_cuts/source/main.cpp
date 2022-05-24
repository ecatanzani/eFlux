#include "main.h"
#include "utils.h"
#include "signal.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                  .......... Prints this help");
	opt.addUsage(" -i  --input                 .......... <path_to_input_list>                                 .......... (*) Input file list");
	opt.addUsage(" -c  --config                .......... <path_to_software_config_dir>                        .......... (*) eFlux energy config file");
	opt.addUsage(" -b  --bdt_cut_values        .......... <path_to_BDT_cut_values_TTree>                       .......... (*) Input BDT cut values TTree");
	opt.addUsage(" -f  --function              .......... <path_to_ROOT_correction_file>                       .......... (*) ROOT efficiency correction file");
	opt.addUsage(" -o  --output                .......... <path_to_output_TFile>                               .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputdir             .......... <path_to_output_TFile_dir>                           .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose               .......... Verbose output");
	
	opt.addUsage("");
	opt.addUsage(" -p  --parallel              .......... <number_of_threads>                                  .......... Multithreading option");
	

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("config", 'c');
	opt.setOption("bdt_cut_values", 'b');
	opt.setOption("function", 'f');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setOption("parallel", 'p');

	opt.processCommandArgs(argc, argv);
	
	// Load args struct
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
	if (opt.getValue("bdt_cut_values") || opt.getValue('b'))
		input_args.bdt_cut_values = opt.getValue('b');
	if (opt.getValue("function") || opt.getValue('f'))
		input_args.eff_corr_function = opt.getValue('f');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_args.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getValue("parallel") || opt.getValue('p'))
		input_args.threads =  std::stoul(opt.getValue('p'), nullptr, 0);

	if (!input_args.output_path.empty())
		signal_efficiency(input_args);
	else
	{
		std::cerr << "\n\nError ! No output path has been specified...\n\n";
		exit(100);
	}
		
	return 0;
}