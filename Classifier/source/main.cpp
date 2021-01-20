#include "main.h"
#include "train.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help													Prints this help");
	opt.addUsage(" -S  --train-signal			<path_to_input_signal_TTree>		(*)	Input signal train TTree");
	opt.addUsage(" -B  --train-background		<path_to_input_background_TTree>	(*) Input background train TTree");
	opt.addUsage(" -s  --test-signal			<path_to_input_signal_TTree>		(*)	Input signal test TTree");
	opt.addUsage(" -b  --test-background		<path_to_input_background_TTree>	(*) Input background test TTree");
	opt.addUsage(" -o  --output					<path_to_output_TFile>					Output ROOT TFile");
	opt.addUsage(" -d  --outputDir				<path_to_output_TFile_dir>				Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose														Verbose output");

	opt.addUsage("");
	opt.addUsage("TMVA learning method: ");
	opt.addUsage("");

	opt.addUsage(" -m  --method					<TMVA_learning_method>					TMVA learning/classifying method");
	opt.addUsage(" -t  --data-test														Test with data");

	opt.addUsage("");

	opt.setFlag("help", 'h');
	opt.setOption("train-signal", 'S');
	opt.setOption("train-background", 'B');
	opt.setOption("test-signal", 's');
	opt.setOption("test-background", 'b');
	opt.setOption("output", 'o');
	opt.setOption("outputDir", 'd');
	opt.setFlag("verbose", 'v');
	opt.setOption("method", 'm');
	opt.setFlag("data-test", 't');
	
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
	if (opt.getValue("train-signal") || opt.getValue('S'))
		input_args.train_signal_path = opt.getValue('S');
	if (opt.getValue("train-background") || opt.getValue('B'))
		input_args.train_background_path = opt.getValue('B');
	if (opt.getValue("test-signal") || opt.getValue('s'))
		input_args.test_signal_path = opt.getValue('s');
	if (opt.getValue("test-background") || opt.getValue('b'))
		input_args.test_background_path = opt.getValue('b');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = opt.getValue('o');
	if (opt.getValue("outputDir") || opt.getValue('d'))
		input_args.output_path = opt.getValue('d');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getValue("method") || opt.getValue('m'))
		input_args.learning = opt.getValue('m');
	if (opt.getFlag("data-test") || opt.getFlag('t'))
		input_args.test_with_data = opt.getFlag('t');

	if (!input_args.output_path.empty())
	{
		if (!input_args.learning.empty())
		{
			if (input_args.verbose)
				PrintArgs(input_args);
			Train(input_args);
		}
		else
		{
			std::cerr << "\nError ! No learnign method has been specified...\n\n";
			exit(100);
		}
	}
	else
	{
		std::cerr << "\nError ! No output path has been specified...\n\n";
		exit(100);
	}
	
	return 0;
}

void PrintArgs(in_args input_args)
{
	std::cout << "\n***** TMVA input arguments *****\n";
	std::cout << "\nSignal input training DSet: [" << input_args.train_signal_path << "]";
	std::cout << "\nBackground input training DSet: [" << input_args.train_background_path << "]";
	std::cout << "\nSignal input test DSet: [" << input_args.test_signal_path << "]";
	std::cout << "\nBackground input test DSet: [" << input_args.test_background_path << "]";
	std::cout << "\nOutput path: [" << input_args.output_path << "]";
	std::cout << "\n\nLearning method: [" << input_args.learning << "]";
	input_args.test_with_data ? std::cout << "\nTest with data: [True]" : std::cout << "\nTest with data: [False]";
	std::cout << "\n\n******************************\n\n";
}