#include "myHeader.h"

#include <ctime>
#include <sstream>

const std::string uniqueOutFile(
    const std::string outputPath, 
    AnyOption &opt)
{
    std::time_t ctime = std::time(0);
    std::stringstream fPath;
    if (opt.getValue("outputDir") || opt.getValue('d'))
        fPath << outputPath << "/analysisOutFile_" << ctime << ".root";
    else if (opt.getValue("output") || opt.getValue('o'))
        fPath << outputPath;
    else
        fPath << "analysisOutFile_" << ctime << ".root";

    return fPath.str();
}

#if 0
const std::string uniqueTupleOutFile(
    AnyOption &opt, 
    const int year, 
    const int month,
    const int emin,
    const int emax)
{
    std::stringstream fPath;
    if (opt.getValue("outputDir") || opt.getValue('d'))
        fPath << opt.getValue("outputDir") << "/DmpNtup_" << year << month << "_" << emin << "_" << emax << ".root";
    else if (opt.getValue("output") || opt.getValue('o'))
        fPath << opt.getValue("outputDir");
    else
        fPath << "DmpNtup_" << year << month << "_" << emin << "_" << emax << ".root";
    return fPath.str();
}
#else
const std::string uniqueTupleOutFile(
    AnyOption &opt, 
    const std::string year, 
    const std::string month)
{
    std::time_t ctime = std::time(0);
    std::stringstream fPath;
    if (opt.getValue("outputDir") || opt.getValue('d'))
        fPath << opt.getValue("outputDir") << "/DmpNtup_" << year << month << "_" << ctime << ".root";
    else if (opt.getValue("output") || opt.getValue('o'))
        fPath << opt.getValue("outputDir");
    else
        fPath << "DmpNtup_" << year << month << "_" << ctime << ".root";
    return fPath.str();
}
#endif