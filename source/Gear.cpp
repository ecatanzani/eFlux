#include "myHeader.h"

std::vector<double> createLogBinning(const double minValue,const double maxValue,const int nBins)
{
    std::vector<double> binning;

    const double minLog = std::log(minValue);
    const double maxLog = std::log(maxValue);

    binning.resize(nBins);

    const double logIncrement = (maxLog - minLog)/nBins;

    double value = minValue ;
    double logValue = minLog ;
    
    for(int bIdx = 0 ; bIdx<nBins ; ++bIdx )
    {
        binning[bIdx] = value;
        logValue += logIncrement;
        value = std::exp(logValue);
    }

    return binning;
}

std::vector<double> createLinBinning(const double minValue,const double maxValue,const int nBins)
{
    std::vector<double> binning;

    binning.resize(nBins);
    
    const double linIncrement = (maxValue - minValue)/nBins;
    double value = minValue;
    
    for(int bIdx = 0 ; bIdx<nBins ; ++bIdx )
    {
        binning[bIdx] = value;
        value += linIncrement;
    }
    
    return binning;
}