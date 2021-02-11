#include "sbicore.h"
#include "buildsbi.h"

void SBIcore(in_pars input_pars)
{
    // Load output file
    TFile outfile(input_pars.output_path.c_str(), "NEW", "SBI Output File");
    if (!outfile.IsOpen())
    {
        std::cerr << "\n\nERROR: output file not created [" << input_pars.output_path << "]\n\n";
        exit(100);
    }

    buildSBI(input_pars, outfile);
    outfile.Close();
}