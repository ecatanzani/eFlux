#include "acceptance_fit.h"

inline TGraphAsymmErrors removeZeros(TGraphAsymmErrors* acceptanceGr)
{
    // Read TGraph number of points
    auto nPoints = acceptanceGr->GetN();
    // Create struct object
    gr_point tmp_point;
    // Create vectors to store graph values
    std::vector<double> energyValues;
    std::vector<double> acceptanceValues;
    std::vector<double> acceptanceErrorLow;
    std::vector<double> acceptanceErrorHigh;
    std::vector<double> energyErrorLow;
    std::vector<double> energyErrorHigh;

    // Read input graph
    for(auto idx=0; idx<nPoints; ++idx)
    {
        acceptanceGr->GetPoint(idx, tmp_point.X, tmp_point.Y);
        if (tmp_point.Y)
        {
            energyValues.push_back(tmp_point.X);
            acceptanceValues.push_back(tmp_point.Y);
            acceptanceErrorLow.push_back(acceptanceGr->GetErrorYlow(idx));
            acceptanceErrorHigh.push_back(acceptanceGr->GetErrorYhigh(idx));
            energyErrorLow.push_back(acceptanceGr->GetErrorXlow(idx));
            energyErrorHigh.push_back(acceptanceGr->GetErrorXhigh(idx));
        }   
    }

    // Create cleaned graph
    TGraphAsymmErrors cleanedAcceptance(energyValues.size(), &energyValues[0], &acceptanceValues[0], &energyErrorLow[0], &energyErrorHigh[0], &acceptanceErrorLow[0], &acceptanceErrorHigh[0]);

    cleanedAcceptance.SetName("gr_acceptance_cleaned");
    cleanedAcceptance.SetTitle("Acceptance - cleaned");

    return cleanedAcceptance;
    
}

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
    auto acceptanceGr_cleaned = removeZeros(acceptanceGr);
    auto fitter = fitAcceptance(&acceptanceGr_cleaned,outFile,verbose);
    return fitter;
}