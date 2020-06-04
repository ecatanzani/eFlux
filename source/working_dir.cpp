#include "working_dir.h"

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