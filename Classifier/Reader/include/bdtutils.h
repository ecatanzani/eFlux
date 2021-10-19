#ifndef BDTUTILS_H
#define BDTUTILS_H

#include <map>
#include <memory>
#include <string>

#include "TMVA/Reader.h"

extern std::map<std::string, int> GetTMVAMethods(const std::string mymethod);

#endif