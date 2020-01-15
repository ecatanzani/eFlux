#include "myHeader.h"

void buildFlux(TFile &outFile)
{
    
    const double minValue = 3e-10;
    const double maxValue = 5e+10;
    const int nBins = 10;

    std::vector<double> logEBins = createLogBinning(minValue,maxValue,nBins);

    buildXtrlFlux(outFile,logEBins);

}

void buildXtrlFlux(TFile &outFile,std::vector<double> &logBins)
{

}