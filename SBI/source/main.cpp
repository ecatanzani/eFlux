#include "main.h"
#include "utils.h"
#include "buildsbi.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -i  --input			<path_to_input_list>				(*)	Input list");
	opt.addUsage(" -o  --output			<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -d  --outputDir		<path_to_output_TFile_dir>				Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose												Verbose output");
	opt.addUsage(" -p  --pedantic												Pedantic output");

	opt.addUsage("");
	
	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setFlag("pedantic", 'p');
	
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
	input_pars.output_path = uniqueOutFile(opt);
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_pars.verbose = opt.getFlag('v');
	if (opt.getFlag("pedantic") || opt.getFlag('p'))
		input_pars.pedantic = opt.getFlag('p');
	
	if (input_pars.CheckArgs())
		buildSBI(input_pars);
	else
		return 100;
	
	return 0;
}