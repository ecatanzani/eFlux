#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>

#include "anyoption.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TChain.h"

// DAMPESW includes

#include "DmpChain.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"

#pragma once

#define nGFitLoops 1000000
#define chiSQLLimit 0.000001

#define GetCurrentDir getcwd
struct energy_cuts
{
    double cut_min_energy;
    double cut_max_energy;
    unsigned int e_bins;
};

extern std::string getWorkingDir(const char *exePath);
extern void load_energy_struct(energy_cuts &e_cuts, const std::string wd);
extern std::vector<float> createLogBinning(const energy_cuts e_cuts);

extern void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const unsigned int lvTime,
    const std::string accInputPath,
    AnyOption &opt,
    const std::string wd);

extern const std::string uniqueOutFile(
    const std::string outputPath,
    AnyOption &opt);

extern std::string GetCurrentWorkingDir(void);

extern std::string getListPath(
    const std::string accInputPath,
    const bool MC = false);

extern bool chechFlags(
    AnyOption &opt,
    const std::string inputPath,
    const std::string outputPath,
    const unsigned int lvTime);

extern void generateFinalGraph(
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    const std::string complete_histo_path,
    const std::string wd);

#endif