#ifndef FLUX_H
#define FLUX_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "data_cuts.h"

extern void buildFlux(
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath,
    const std::string wd);

extern std::vector<float> load_flux_struct(
    cuts_conf &flux_cuts,
    data_active_cuts &active_cuts,
    const std::string wd);

#endif