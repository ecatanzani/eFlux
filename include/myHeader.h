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

extern void eCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    const bool pedantic,
    const unsigned int lvTime,
    const bool myAcceptance,
    const std::string accInputPath,
    AnyOption &opt);

extern void readInputTree(
    const std::string inputPath,
    std::vector<float> &dataValues,
    TTree &dTree);

extern void branchTree(TTree &myDataTree, std::vector<float> &dataValues);
extern const char *uniqueOutFile(const std::string outputPath, AnyOption &opt);
extern std::vector<float> createLogBinning(const double minValue, const double maxValue, const int nBins);
extern std::vector<float> createLinBinning(const double minValue, const double maxValue, const int nBins);

extern std::string getListPath(const std::string accInputPath, const bool MC = false);

extern void buildFlux(
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose,
    const bool pedantic,
    const std::string accInputPath,
    const bool myAcceptance);

extern bool chechFlags(
    AnyOption &opt,
    const std::string inputPath,
    const std::string outputPath,
    const unsigned int lvTime);

extern void buildXtrlFlux(
    std::vector<float> &eBins,
    std::vector<float> &cBins,
    const std::string inputPath,
    const unsigned int lvTime,
    TFile &outFile,
    const bool verbose);

extern void evLoop(
    TH1D &inputHisto,
    const std::string inputPath,
    TFile &outFile,
    const bool verbose,
    const bool eClassifier = false,
    const double xtrlCut = 8.5);

extern void readAcceptance(TH1D &acceptance, TFile &outFile, const bool verbose);

#endif