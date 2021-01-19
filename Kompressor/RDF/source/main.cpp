#include "main.h"
#include "utils.h"
#include "reader.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -i  --input			<path_to_input_list>				(*)	Input file list");
	opt.addUsage(" -w  --workdir		<path_to_software_config_dir>		(*) Config directory");
	opt.addUsage(" -o  --output			<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -d  --outputdir		<path_to_output_TFile_dir>				Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose												Verbose output");
	
	opt.addUsage("");
	opt.addUsage("Tasks: ");
	opt.addUsage("");

	opt.addUsage(" -m  --mc														MC event loop");
	opt.addUsage(" -r  --regularize		<path_to_summary_fit_TTree>				Regularize variables behaviour");
	opt.addUsage(" -t  --tmva-set		<s(signal)/b(background)>				Create TMVA Test/Training sets");
	

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputdir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("mc", 'm');
	opt.setOption("regularize", 'r');
	opt.setOption("tmva-set", 't');

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
		input_args.wd = opt.getValue('w');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = expand_output_path(opt, opt.getValue('o'));
	if (opt.getValue("outputdir") || opt.getValue('d'))
		input_args.output_path = expand_output_path(opt, opt.getValue('d'));
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args._VERBOSE = opt.getFlag('v');
	
	if (opt.getFlag("mc") || opt.getFlag('m'))
		input_args.mc_flag = true;
	if (opt.getValue("regularize") || opt.getValue('r'))
		input_args.fit_tree_path = opt.getValue('r');
	if (opt.getValue("tmva-set") || opt.getValue('t'))
	{
		input_args.tmva_set = true;
		input_args.SetSeType(opt.getValue('t'));
	}

	if (!input_args.output_path.empty())
		reader(input_args);
	else
	{
		std::cerr << "\n\nError ! No output path has been specified...\n\n";
		exit(100);
	}
		
	return 0;
}