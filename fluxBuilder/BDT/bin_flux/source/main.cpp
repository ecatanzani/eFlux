#include "main.h"
#include "utils.h"
#include "bin_flux.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                  .......... Prints this help");
	opt.addUsage(" -i  --input                 .......... <path_to_DATA_input_list>                            .......... (*) Input DATA file list");
	opt.addUsage(" -s  --simu                  .......... <path_to_electron_MC_input_list>                     .......... (*) Input electron MC file list");
	opt.addUsage(" -c  --config                .......... <path_to_software_config_dir>                        .......... (*) eFlux energy config file");
	opt.addUsage(" -b  --config-bdt            .......... <path_to_bdt_config_file>                            .......... BDT config file");
	opt.addUsage(" -a  --acceptance            .......... <path_to_acceptance_ROOT_file>                       .......... (*) DAMPE acceptance ROOT file");
	opt.addUsage(" -e  --exposure              .......... <exposure_time>                                      .......... (*) DAMPE exposure time");
	opt.addUsage(" -r  --energy-range          .......... <energy_bin>                                         .......... (*) Energy Bin Number");
	opt.addUsage(" -f  --function              .......... <path_to_ROOT_correction_file>                       .......... (*) ROOT efficiency correction file");
	opt.addUsage(" -o  --output                .......... <path_to_output_TFile>                               .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputdir             .......... <path_to_output_TFile_dir>                           .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose               .......... Verbose output");
	
	opt.addUsage("");
	opt.addUsage(" -p  --parallel              .......... <number_of_threads>                                  .......... Multithreading option");

	opt.addUsage("");
	opt.addUsage("TMVA learning method: ");
	opt.addUsage("");

	opt.addUsage(" -m  --method                      .......... <TMVA_learning_method>                           .......... TMVA learning/classifying method");
	
	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("simu", 's');
	opt.setOption("config", 'c');
	opt.setOption("config-bdt", 'b');
	opt.setOption("acceptance", 'a');
	opt.setOption("exposure", 'e');
	opt.setOption("energy-range", 'r');
	opt.setOption("function", 'f');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setOption("parallel", 'p');
	opt.setOption("method", 'm');
	
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
	if (opt.getValue("input") || opt.getValue('i')) 			input_args.input_list = opt.getValue('i');
	if (opt.getValue("simu") || opt.getValue('s')) 				input_args.simu_input_list = opt.getValue('s');
	if (opt.getValue("config") || opt.getValue('c')) 			input_args.energy_config_file = opt.getValue('c');
	if (opt.getValue("config-bdt") || opt.getValue('b')) 		input_args.bdt_config_file = opt.getValue('b');
	if (opt.getValue("acceptance") || opt.getValue('a')) 		input_args.acceptance_file = opt.getValue('a');
	if (opt.getValue("exposure") || opt.getValue('e')) 			input_args.exposure = std::stod(opt.getValue('e'));
	if (opt.getValue("energy-range") || opt.getValue('r')) 		input_args.energy_bin = std::stoul(opt.getValue('r'));
	if (opt.getValue("function") || opt.getValue('f'))			input_args.eff_corr_function = opt.getValue('f');
	if (opt.getValue("output") || opt.getValue('o')) 			input_args.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d')) 		input_args.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v')) 			input_args.verbose = opt.getFlag('v');
	if (opt.getValue("parallel") || opt.getValue('p'))			input_args.threads =  std::stoul(opt.getValue('p'), nullptr, 0);
	if (opt.getValue("method") || opt.getValue('m'))			input_args.learning_method = opt.getValue('m');

	if (!input_args.output_path.empty())
		bin_flux(input_args);
	else
	{
		std::cerr << "\n\nError ! No output path has been specified...\n\n";
		exit(100);
	}
		
	return 0;
}