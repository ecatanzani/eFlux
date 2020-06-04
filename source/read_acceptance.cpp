#include "myHeader.h"
#include "acceptance_fit.h"

TF1 readAcceptance(
    TFile &outFile, 
    const bool verbose,
    const std::string accInputPath)
{
    // Reading acceptance
    TFile accInputFile(accInputPath.c_str(), "READ");
    if (!accInputFile.IsOpen())
    {
        std::cerr << "\n\nError opening input acceptance TFile: " << accInputPath << std::endl;
        exit(123);
    }
    auto acceptanceGr = static_cast<TGraphAsymmErrors*> (accInputFile.Get("Acceptance/gr_acceptance_all_cut"));
    auto fitter = fitAcceptance<TGraphAsymmErrors>(acceptanceGr,outFile,verbose);
    return fitter;
}