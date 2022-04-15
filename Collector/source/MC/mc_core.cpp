#include "MC/mc.h"

#include <memory>

void mcCore(in_pars input_pars)
{
	// Create output TFile
	if (input_pars.verbose)
		std::cout << "\nCreating output ROOT file... [" << input_pars.output_path << "]" << std::endl;
	TFile out_file(input_pars.output_path.c_str(), "NEW", "Analysis Output File");
	if (!out_file.IsOpen())
	{
		std::cerr << "\n\nError writing output TFile: " << input_pars.output_path << std::endl;
		exit(100);
	}
	
	mcLoop(
		input_pars.input_path,
		out_file,
		input_pars.verbose,
		input_pars.wd,
		input_pars.stk_correction);

	// Close output file ...
	out_file.Close();
}