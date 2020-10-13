#ifndef NTUPLE_H
#define NTUPLE_H

#include <string>

#include "anyoption.h"

#include "TFile.h"

extern void read_ntuple(
    AnyOption &opt,
    const std::string inputPath, 
    const std::string outputPath, 
    const std::string wd,
    const bool verbose,
    const bool mthreads=true);

#endif