#ifndef FLUX_H
#define FLUX_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"

struct flux_conf
{
    double min_energy;
    double max_energy;
    unsigned int n_energy_bins;
};

extern std::vector<float> createLogBinning(const flux_conf flux_params);

extern void buildFlux(
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const bool pedantic,
    const std::string accInputPath,
    const bool myAcceptance);

extern void buildXtrlFlux(
    const std::vector<float> &eBins,
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const bool myAcceptance);

extern void load_flux_struct(flux_conf &flux_params);

#endif