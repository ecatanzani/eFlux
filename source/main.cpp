#include "myHeader.h"
#include "acceptance.h"

#include <sstream>

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");
	opt.addUsage(" -h  --help                                                   		Prints this help");
	opt.addUsage(" -i  --input          <path_to_input_DATA_list>           		(*) Input DATA list");
	opt.addUsage(" -w  --workdir        <path_to_software_config_dir>       		(*) Config directory");
	opt.addUsage(" -o  --output         <path_to_output_TFile>                          Output ROOT TFile");
	opt.addUsage(" -d  --outputDir      <path_to_output_TFile_dir>                      Output ROOT TFile directory");
	opt.addUsage(" -t  --lvtime         <live_time-value>                   		(*) DAMPE live-time ");
	opt.addUsage(" -a  --acceptance     <path_to_MC_list>                   		(*) Acceptance calculation");
	opt.addUsage(" -n  --ntuple                                             		(*) nTuple production facility");
	opt.addUsage(" -c  --collect        <path_to_complete_histo>            		(*) Generate TGraph from final histo");
	opt.addUsage(" -r  --resume         <mc/data>                           		(*) MC/DATA facility");
	opt.addUsage(" -g  --geometry       <path_to_electron_acceptance_ROOT_file>     (*) Electron Acceptance file - flux calculation only");
	opt.addUsage(" -b  --background     <path_to_proton_acceptance_ROOT_file>       (*) Proton Acceptance file - flux calculation only");
	opt.addUsage(" -v  --verbose                                                        Verbose output");
	opt.addUsage(" -p  --pedantic                                                       Pedantic output");
	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("input", 'i');
	opt.setOption("workdir", 'w');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setOption("lvtime", 't');
	opt.setOption("acceptance", 'a');
	opt.setOption("geometry", 'g');
	opt.setOption("background", 'b');
	opt.setOption("collect", 'c');
	opt.setOption("resume", 'r');
	opt.setFlag("verbose", 'v');
	opt.setFlag("pedantic", 'p');
	opt.setFlag("ntuple", 'n');

	opt.processCommandArgs(argc, argv);

	/*
        Input variables

    */

	std::string wd;
	std::string inputPath;
	std::string outputPath;
	std::string accInputPath;
	std::string p_accInputPath;
	std::string inputCompleteHisto;
	stringstream str_lvTime;

	bool verbose = false;
	bool pedantic = false;
	bool myAcceptance = false;
	bool myFlux = false;
	bool myNtuple = false;
	unsigned int lvTime = 0;
	bool genGraph = false;
	bool pasteMC = false;
	bool pasteDATA = false;

	if (!opt.hasOptions())
		opt.printUsage();

	if (opt.getFlag("help") || opt.getFlag('h'))
		opt.printUsage();
	if (opt.getValue("input") || opt.getValue('i'))
	{
		inputPath = opt.getValue('i');
		myFlux = true;
	}
	if (opt.getValue("workdir") || opt.getValue('w'))
		wd = opt.getValue('w');
	if (opt.getValue("output") || opt.getValue('o'))
		outputPath = opt.getValue('o');
	if (opt.getValue("outputDir") || opt.getValue('d'))
		outputPath = opt.getValue('d');
	if (opt.getValue("lvtime") || opt.getValue('t'))
	{
		str_lvTime << opt.getValue('t');
		str_lvTime >> lvTime;
	}
	if (opt.getValue("acceptance") || opt.getValue('a'))
	{
		myAcceptance = true;
		accInputPath = opt.getValue('a');
	}
	if (opt.getValue("ntuple") || opt.getValue('n'))
		myNtuple = true;
	if (opt.getValue("geometry") || opt.getValue('g'))
		accInputPath = opt.getValue('g');
	if (opt.getValue("background") || opt.getValue('b'))
		p_accInputPath = opt.getValue('b');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		verbose = opt.getFlag('v');
	if (opt.getValue("collect") || opt.getValue('c'))
	{
		genGraph = true;
		inputCompleteHisto = opt.getValue('c');
	}
	if (opt.getValue("resume") || opt.getValue('r'))
	{
		if (!strcmp(opt.getValue('r'), "mc"))
			pasteMC = true;
		else if (!strcmp(opt.getValue('r'), "data"))
			pasteDATA = true;
		else
		{
			std::cerr << "\n\nWrong --resume config...";
			opt.printUsage();
			exit(100);
		}
	}
	if (opt.getFlag("pedantic") || opt.getFlag('p'))
		pedantic = opt.getFlag('p');

	if (myAcceptance)
		computeAcceptance(
			accInputPath,
			verbose,
			pedantic,
			outputPath,
			opt,
			wd);

	if (genGraph)
	{
		if (pasteMC)
			generateFinalGraph(
				verbose,
				pedantic,
				outputPath,
				inputCompleteHisto,
				wd);
				
		if (pasteDATA)
			generateDataFinalGraph(
				verbose,
				pedantic,
				outputPath,
				inputCompleteHisto,
				wd);
	}

	if (myFlux || myNtuple)
		eCore(
			inputPath,
			outputPath,
			verbose,
			pedantic,
			lvTime,
			accInputPath,
			p_accInputPath,
			opt,
			wd);

	return 0;
}