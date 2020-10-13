#include "main.h"
#include "extractor.h"
#include "RooFitter.h"

int main(int argc, char **argv)
{
	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -d  --data			<path_to_input_DATA_root_file>	    (*)	Input DATA root file");
	opt.addUsage(" -e  --electron		<path_to_input_electron_MC_file>    (*) Input electron MC file");
	opt.addUsage(" -p  --proton		    <path_to_input_proton_MC_file>		(*) Input proton MC file");
	opt.addUsage(" -o  --output			<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -v  --verbose												Verbose output");

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("data", 'd');
	opt.setOption("electron", 'e');
	opt.setOption("proton", 'p');
	opt.setOption("output", 'o');
	opt.setFlag("verbose", 'v');

	opt.processCommandArgs(argc, argv);
	deps_paths paths;
	std::string output_path;
	bool verbose = false;

	if (!opt.hasOptions() || opt.getFlag("help") || opt.getFlag('h'))
	{
		opt.printUsage();
		return 0;
	}
	if (opt.getValue("data") || opt.getValue('d'))
		paths.data_path = opt.getValue('d');
	if (opt.getValue("electron") || opt.getValue('e'))
		paths.electron_path = opt.getValue('e');
	if (opt.getValue("proton") || opt.getValue('p'))
		paths.proton_path = opt.getValue('p');
	if (opt.getValue("output") || opt.getValue('o'))
		output_path = opt.getValue('o');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		verbose = opt.getFlag('v');

	// Reading deps files
	extractor deps(paths, verbose);

	// Call RooFitter class
	RooFitter fitter(
		deps.GetData(),
		deps.GetElectronTemplates(),
		deps.GetProtonTemplates(),
		deps.GetBins(),
		_RooFitVerbosity,
		_RooFitClamping);
	// Perform fit
	fitter.Fit();
	// Write results
	fitter.SaveResults(output_path);

	return 0;
}