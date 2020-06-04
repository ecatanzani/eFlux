#ifndef WORKINGDIR_H
#define WORKINGDIR_H

#include <string>
#include <unistd.h>

#define GetCurrentDir getcwd

extern std::string getWorkingDir(const char *exePath);
extern std::string GetCurrentWorkingDir(void);

#endif