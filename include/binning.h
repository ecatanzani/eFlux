#ifndef BINNING_H
#define BINNING_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "energyCuts.h"

extern std::vector<float> createLogBinning(const energy_cuts e_cuts);
extern std::vector<float> readLogBinning(const std::string wd);

#endif