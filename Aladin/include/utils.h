#ifndef UTILS_H
#define UTILS_H

#include <string>

#include "anyoption.h"

extern std::string expand_output_path(
    AnyOption& opt, 
    std::string output);

extern std::string st1_path(std::string output);
extern std::string st2_path(std::string output);

#endif
