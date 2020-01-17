#include "myHeader.h"

void buildFlux(TFile &outFile,const std::string inputPath)
{
    
    double minValue = 3e-10;
    double maxValue = 5e+10;
    const int nBins = 10;

    std::vector<double> logEBins = createLogBinning(minValue,maxValue,nBins);
    
    minValue = 10;
    maxValue = 1e+5;
    
    std::vector<double> eCounts = createLinBinning(minValue,maxValue,nBins);

    buildXtrlFlux(outFile,logEBins,eCounts,inputPath);

}

void buildXtrlFlux(TFile &outFile,std::vector<double> &eBins,std::vector<double> &cBins,const std::string inputPath)
{
    /*
        All-Electron flux using xtrl as classifier
    */

    TH1D eCounts("eCounts","All Electron counts - xtrl classifier",eBins.size(),&(eBins[0]));
    evLoop(outFile,eCounts,inputPath,true);


}