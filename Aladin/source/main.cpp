#include "main.h"
#include "utils.h"
#include "reader.h"

int main(int argc, char **argv)
{
	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help                  .......... Prints this help");
	opt.addUsage(" -i  --input                 .......... <path_to_input_list>                                 .......... (*) Input file list");
	opt.addUsage(" -w  --workdir               .......... <path_to_software_config_dir>                        .......... (*) Collector config directory");
	opt.addUsage(" -o  --output                .......... <path_to_output_TFile>                               .......... Output ROOT TFile");
	opt.addUsage(" -d  --outputdir             .......... <path_to_output_TFile_dir>                           .......... Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose               .......... Verbose output");
	
	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc                    .......... MC event loop");
	opt.addUsage(" -r  --regularize            .......... <path_to_summary_fit_TTree>                          .......... Regularize variables behaviour");
	opt.addUsage("");
	opt.addUsage(" -g  --gaussianize           .......... Gaussianize input TMVA variables");
	opt.addUsage(" -l  --likelihood            .......... Build LogLikelihood profile");
	opt.addUsage(" -f  --fit                   .......... Fit gaussianized input TMVA variables");
	
	opt.addUsage("");
	opt.addUsage(" -b  --bin                   .......... Energy bin");

	opt.addUsage("");
	opt.addUsage(" -t  --tmva                  .......... <path_to_best_lambda_Tree>                           .......... Regularize variables behaviour");

	opt.addUsage("");
	opt.addUsage(" -p  --parallel              .......... <number_of_threads>                                  .......... Multithreading option");
	

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');
	opt.setOption("regularize", 'r');
	opt.setFlag("gaussianize", 'g');
	opt.setFlag("fit", 'f');
	opt.setOption("bin", 'b');
	opt.setOption("tmva", 't');
	opt.setFlag("likelihood", 'l');
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
	if (opt.getValue("workdir") || opt.getValue('w'))
		input_args.collector_wd = opt.getValue('w');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_args.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_args.mc_flag = true;
	if (opt.getValue("regularize") || opt.getValue('r'))
		input_args.regularize_tree_path = opt.getValue('r');
	if (opt.getFlag("gaussianize") || opt.getFlag('g'))
		input_args.gaussianize = opt.getFlag('g');
	if (opt.getFlag("fit") || opt.getFlag('f'))
		input_args.fit = opt.getFlag('f');
	if (opt.getValue("bin") || opt.getValue('b'))
		input_args.energybin = std::stoul(opt.getValue('b'), nullptr, 0);
	if (opt.getValue("tmva") || opt.getValue('t'))
		input_args.best_lambda_tree_path = opt.getValue('t');
	if (opt.getFlag("likelihood") || opt.getFlag('l'))
		input_args.loglikelihood = opt.getFlag('l');
	if (opt.getValue("parallel") || opt.getValue('p'))
		input_args.threads =  std::stoul(opt.getValue('p'), nullptr, 0);

	// Set the energy bin
	if (!input_args.output_path.empty() && input_args.check_input())
		reader(input_args);
	else
	{
		std::cerr << "\n\nError ! Please check input parameters...\n\n";
		exit(100);
	}
		
	return 0;
}