#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"

#include "anyoption.h"
#include "data_cuts.h"

struct mc_ancillary_cuts
{
    bool compute_proton_background = false;
};

extern void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    AnyOption &opt,
    const std::string wd);

extern void buildAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const std::vector<float> &logEBins,
    TFile &outFile,
    const std::string wd);

extern void buildAcceptance_vector(
    const std::string accInputPath,
    const bool verbose,
    const std::vector<float> &logEBins,
    TFile &outFile,
    const std::string wd);

extern TF1 readAcceptance(
    TFile &outFile,
    const bool verbose,
    const std::string accInputPath);

extern void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic);

extern void load_acceptance_struct(
    cuts_conf &acceptance_cuts,
    data_active_cuts &active_cuts,
    mc_ancillary_cuts &ancillary_cuts,
    const std::string wd);

#endif