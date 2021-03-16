#include <iostream>

#include "TMVA/TMVAGui.h"

void TMVAGui(std::string tmva_outfile)
{
    if (!tmva_outfile.empty())
        TMVA::TMVAGui(tmva_outfile.c_str());
    else
        std::cerr << "\n\nEmpty TMVA outfile, please specify one...\n\n";
} 