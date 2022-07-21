#include "main.h"
#include "preselection.h"

#include "anyoption.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                    .......... Prints this help");
	opt.addUsage(" -i  --input                   .......... <path_to_input_list>                   .......... (*) Input list");
	opt.addUsage(" -d  --outputDir               .......... <path_to_output_TFile_dir>             .......... Output ROOT TFile directory");
	opt.addUsage(" -e  --energy                  .......... <path_to_energy_config_file>           .......... (*) Energy config file");
	opt.addUsage(" -s  --spectral                .......... <input_spectral_index>                 .......... (*) Input spectral index if MC");
	opt.addUsage(" -v  --verbose                 .......... Verbose output");

	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc                      .......... MC event loop");

	opt.addUsage("");
	opt.addUsage(" -p  --parallel                .......... <number_of_threads>                    .......... Multithreading option");

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("outputDir", 'd');
	opt.setOption("energy", 'e');
	opt.setOption("spectral", 's');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');
	opt.setFlag("raw_data", 'r');
	opt.setOption("parallel", 'p');
	
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
	if (opt.getValue("outputDir") || opt.getValue('d')) {
		input_pars.output_directory = opt.getValue('d');
		input_pars.GetOutputFile();
	}
	if (opt.getValue("energy") || opt.getValue('e'))
		input_pars.energy_config_file = opt.getValue('e');
	if (opt.getValue("spectral") || opt.getValue('s'))
		input_pars.input_spectral_index = std::stod(opt.getValue('s'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_pars.verbose = opt.getFlag('v');
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_pars.mc_flag = true;
	if (opt.getValue("parallel") || opt.getValue('p'))
		input_pars.threads =  std::stoul(opt.getValue('p'), nullptr, 0);

	
	if (input_pars.CheckArgs())
	{
		if (input_pars.verbose) {
			input_pars.ExpandArgs();
        }
        preselection(
			input_pars.input_path.c_str(), 
			input_pars.output_path.c_str(), 
			input_pars.verbose, 
			input_pars.mc_flag, 
			input_pars.threads,
			input_pars.energy_config_file.c_str(),
			input_pars.input_spectral_index);
	}	
	else
		return 100;
	
	return 0;
}