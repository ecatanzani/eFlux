#include "myHeader.h"

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