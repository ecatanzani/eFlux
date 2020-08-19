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
#include "anyoption.h"

extern void buildFlux(
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath,
    const TH1D* electron_acceptance,
    const TH1D* proton_background_fraction,
    const std::string wd);

extern void produceTuples(
    AnyOption &opt,
    const std::string inputPath,
    const bool verbose,
    const std::string wd);

extern std::vector<float> load_flux_struct(
    cuts_conf &flux_cuts,
    data_active_cuts &active_cuts,
    const std::string wd);

#endif