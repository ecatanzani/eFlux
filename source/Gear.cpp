#include "myHeader.h"

#include <fstream>
#include <sstream>

std::string getWorkingDir(const char* exePath)
{
    std::string tmpPath(exePath);
    std::size_t index = tmpPath.find("eFlux");
    auto wd = tmpPath.substr(0, index + 5);
    return wd;
}

std::string GetCurrentWorkingDir(void)
{
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
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

std::vector<float> createLogBinning(const energy_cuts e_cuts)
{
    std::vector<float> binning(e_cuts.e_bins + 1, 0);

    const double minLog = std::log(e_cuts.cut_min_energy);
    const double maxLog = std::log(e_cuts.cut_max_energy);
    const double logIncrement = (maxLog - minLog) / e_cuts.e_bins;
    double value = e_cuts.cut_min_energy;
    double logValue = minLog;

    for (auto it = binning.begin(); it != (binning.end() - 1); ++it)
    {
        auto bIdx = std::distance(binning.begin(), it);
        binning[bIdx] = value;
        logValue += logIncrement;
        value = std::exp(logValue);
    }
    binning[e_cuts.e_bins] = value;

    return binning;
}

std::vector<float> readLogBinning(const std::string wd)
{
    std::string cwd = wd;
    std::size_t index = cwd.find("eFlux");
    std::string configPath = cwd.substr(0, index + 5);
    configPath += "/config/binningConfig.txt";
    std::ifstream input_file(configPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << configPath << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    std::vector<float> binning;
    while (input_stream >> tmp_str)
        binning.push_back(stod(tmp_str, &sz));
    return binning;
}