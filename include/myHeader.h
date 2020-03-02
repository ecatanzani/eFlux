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

// DAMPE struct

#define DAMPE_bgo_nLayers 14

#define GetCurrentDir getcwd

extern void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const unsigned int lvTime,
    const bool myAcceptance,
    const std::string accInputPath,
    AnyOption &opt);

extern const char *uniqueOutFile(
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

#endif