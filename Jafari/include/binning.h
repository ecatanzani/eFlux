#ifndef BINNING_H
#define BINNING_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>

extern std::vector<float> read_binning_from_config(
    const std::string wd);

#endif