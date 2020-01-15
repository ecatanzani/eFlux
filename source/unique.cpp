#include "myHeader.h"

#include <ctime>
#include <sstream>

const char* uniqueOutFile(const std::string outputPath)
{
    std::time_t ctime = std::time(0);
    std::stringstream fPath;
    fPath << outputPath << "/analysisOutFile_" << ctime;
    return (fPath.str()).c_str();
}