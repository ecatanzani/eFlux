#include "main.h"
#include "train.h"

int main(int argc, char **argv)
{

	AnyOption opt;

	opt.addUsage("Usage: ");
	opt.addUsage("");

	opt.addUsage(" -h  --help																Prints this help");
	opt.addUsage(" -S  --train-signal			<path_to_input_signal_TTree_list>		(*)	Input signal train list");
	opt.addUsage(" -B  --train-background		<path_to_input_background_TTree_list>	(*) Input background train list");
	opt.addUsage(" -s  --test-signal			<path_to_input_signal_TTree_list>		(*)	Input signal test list");
	opt.addUsage(" -b  --test-background		<path_to_input_background_TTree_list>	(*) Input background test list");
	opt.addUsage(" -o  --output					<path_to_output_TFile>						Output ROOT TFile");
	opt.addUsage(" -d  --outputDir				<path_to_output_TFile_dir>					Output ROOT TFile directory");
	opt.addUsage(" -v  --verbose															Verbose output");

	opt.addUsage("");
	opt.addUsage("TMVA learning method: ");
	opt.addUsage("");

	opt.addUsage(" -m  --method					<TMVA_learning_method_1:method_2...>					TMVA learning/classifying method");
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
		input_args.train_signal_input_list = opt.getValue('S');
	if (opt.getValue("train-background") || opt.getValue('B'))
		input_args.train_background_input_list = opt.getValue('B');
	if (opt.getValue("test-signal") || opt.getValue('s'))
		input_args.test_signal_input_list = opt.getValue('s');
	if (opt.getValue("test-background") || opt.getValue('b'))
		input_args.test_background_input_list = opt.getValue('b');
	if (opt.getValue("output") || opt.getValue('o'))
		input_args.output_path = opt.getValue('o');
	if (opt.getValue("outputDir") || opt.getValue('d'))
		input_args.output_path = opt.getValue('d');
	if (opt.getFlag("verbose") || opt.getFlag('v'))
		input_args.verbose = opt.getFlag('v');
	if (opt.getValue("method") || opt.getValue('m'))
		input_args.ParseLearningMethods(opt.getValue('m'));
	if (opt.getFlag("data-test") || opt.getFlag('t'))
		input_args.test_with_data = opt.getFlag('t');
	
	if (input_args.test_input_lists())
	{
		if (!input_args.learning_method.empty())
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
		std::cerr << "\nError ! Not all the requested input lists have been provided...\n\n";
		exit(100);
	}
	
	return 0;
}

void PrintArgs(in_args input_args)
{	
	auto print_methods = [](const std::vector<std::string> methods) -> std::string 
	{ 
		std::string ret;
		if (methods.size())
			for (unsigned int elm=0; elm< methods.size(); ++elm) 
				if (!elm)
					ret = methods[elm];
				else
					ret += std::string(":") + methods[elm];
		else
			ret = methods[0];
		return ret;
	};

	std::cout << "\n***** TMVA input arguments *****\n";
	std::cout << "\nSignal input training DSet list: [" << input_args.train_signal_input_list << "]";
	std::cout << "\nBackground input training DSet list: [" << input_args.train_background_input_list << "]";
	std::cout << "\nSignal input test DSet list: [" << input_args.test_signal_input_list << "]";
	std::cout << "\nBackground input test DSet list: [" << input_args.test_background_input_list << "]";
	std::cout << "\nOutput path: [" << input_args.output_path << "]";
	
	std::cout << "\n\nLearning method: [" << print_methods(input_args.learning_method) << "]";
	

	input_args.test_with_data ? std::cout << "\nTest with data: [True]" : std::cout << "\nTest with data: [False]";
	std::cout << "\n\n******************************\n\n";
}