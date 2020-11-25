#ifndef UTILS_H
#define UTILS_H

#include <string>

#include "anyoption.h"

extern std::string expand_output_path(
    AnyOption& opt, 
    std::string output);

#endif
