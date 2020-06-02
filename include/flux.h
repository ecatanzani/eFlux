#ifndef FLUX_H
#define FLUX_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"

extern void buildFlux(
    const std::vector<float> &eBins,
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath);

#endif