#include "main.h"
#include "utils.h"
#include "efficiency.h"

int main(int argc, char **argv)
{
	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                  .......... Prints this help");
	opt.addUsage(" -i  --input                 .......... <path_to_input_list>                     .......... Input file list");
	opt.addUsage(" -o  --output                .......... <path_to_output_TFile>                   .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputdir             .......... <path_to_output_TFile_dir>               .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose               .......... Verbose output");
	
	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");
	
	opt.addUsage(" -p  --parallel              .......... <number_of_threads>                      .......... Multithreading option");
	
	opt.addUsage("");
	opt.addUsage("TMVA learning method: ");
	opt.addUsage("");

	opt.addUsage(" -b  --config-bdt              .......... <path_to_bdt_config_file>              .......... BDT config file");
	opt.addUsage(" -e  --config-energy           .......... <path_to_energy_config_file>           .......... eFlux energy config file");
	opt.addUsage(" -l  --learning-method         .......... <TMVA_learning_method>                 .......... TMVA learning/classifying method");
	opt.addUsage(" -c  --cosine-regularize       .......... <path_to_proton_summary_fit_TTree>     .......... Regularize angular variables behaviour");
	opt.addUsage(" -t  --box-cox-regularize      .......... <path_to_best_box-cox_lambda_Tree>     .......... Regularize lambda variables behaviour");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setOption("parallel", 'p');
	opt.setOption("config-bdt", 'b');
	opt.setOption("config-energy", 'e');
	opt.setOption("learning-method", 'l');
	opt.setOption("cosine-regularize", 'c');
	opt.setOption("box-cox-regularize", 't');

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
	if (opt.getValue("config-bdt") || opt.getValue('b'))
		input_args.bdt_config_file = opt.getValue('b');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_args.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getValue("parallel") || opt.getValue('p'))
		input_args.threads =  std::stoul(opt.getValue('p'), nullptr, 0);
	if (opt.getValue("config-energy") || opt.getValue('e'))
		input_args.energy_config_file = opt.getValue('e');
	if (opt.getValue("learning-method") || opt.getValue('l'))
		input_args.bdt_learning_method = opt.getValue('l');
	if (opt.getValue("cosine-regularize") || opt.getValue('c'))
		input_args.cosine_regularize_path = opt.getValue('c');
	if (opt.getValue("box-cox-regularize") || opt.getValue('t'))
		input_args.box_cox_regularize_path = opt.getValue('t');

	// Set the energy bin
	if (!input_args.output_path.empty())
		buildEfficiency(input_args);
	else
	{
		std::cerr << "\n\nError ! Please check input parameters...\n\n";
		exit(100);
	}
		
	return 0;
}