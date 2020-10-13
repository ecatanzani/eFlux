#ifndef READ_PARALLEL_H
#define READ_PARALLEL_H

#include <string>
#include <memory>

#include "TFile.h"
#include "TChain.h"
#include <ROOT/TSeq.hxx>
#include "TROOT.h"

extern void read_parallel(
    const std::string inputPath,
    TFile& outFile,
    const std::string wd,
    const bool verbose);

#endif