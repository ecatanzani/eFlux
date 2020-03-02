#ifndef DATALOOP_H
#define DATALOOP_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#define kStep 10000

extern std::vector<unsigned int> evLoop(
    const std::vector<float> &eBins,
    const std::string inputPath,
    TFile &outFile,
    const bool verbose,
    const bool eClassifier = false,
    const double xtrlCut = 8.5);

extern void readInputTree(
    const std::string inputPath,
    std::vector<float> &dataValues,
    TTree &dTree);

extern std::shared_ptr<TTree> getDataTree(
    const std::string inputPath,
    std::vector<float> &dataValues);

extern void branchTree(
    TTree &myDataTree,
    std::vector<float> &dataValues);

#endif