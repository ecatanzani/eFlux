#include "myHeader.h"

std::vector<float> createLogBinning(const double minValue, const double maxValue, const int nBins)
{
    std::vector<float> binning;

    const double minLog = std::log(minValue);
    const double maxLog = std::log(maxValue);

    binning.resize(nBins + 1);

    const double logIncrement = (maxLog - minLog) / nBins;

    double value = minValue;
    double logValue = minLog;

    for (int bIdx = 0; bIdx < nBins; ++bIdx)
    {
        binning[bIdx] = value;
        logValue += logIncrement;
        value = std::exp(logValue);
    }
    binning[nBins] = value;

    return binning;
}

std::vector<float> createLinBinning(const double minValue, const double maxValue, const int nBins)
{
    std::vector<float> binning;

    binning.resize(nBins + 1);

    const double linIncrement = (maxValue - minValue) / nBins;
    double value = minValue;

    for (int bIdx = 0; bIdx < nBins; ++bIdx)
    {
        binning[bIdx] = value;
        value += linIncrement;
    }
    binning[nBins] = value;

    return binning;
}

bool chechFlags(
    AnyOption &opt,
    const std::string inputPath,
    const std::string outputPath,
    const unsigned int lvTime)
{
    bool status = true;

    if (opt.getValue("input") || opt.getValue('i'))

        if (inputPath.empty())
        {
            status *= false;
            std::cerr << "\n\t !!! Empty input value !\n";
            opt.printUsage();
        }
    if (outputPath.empty())
    {
        status *= false;
        if (opt.getValue("outputDir") || opt.getValue('d'))
            std::cerr << "\n\t Empty output directory value -- this value is needed if `-d` flag is used \n";
        if (opt.getValue("output") || opt.getValue('o'))
            std::cerr << "\n\t !!! Empty output value -- Set to the default value `myAnalysisOut` !\n";
        opt.printUsage();
    }
    if (lvTime == 0)
    {
        status *= false;
        std::cerr << "\n\t !!! Empty live-time value !\n";
        opt.printUsage();
    }

    return status;
}