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
#include "TClonesArray.h"

// DAMPESW includes
#include "DmpChain.h"
#include "DmpBgoContainer.h"
#include "DmpEvtHeader.h"
#include "DmpEvtSimuPrimaries.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdHits.h"

extern TH1D evLoop(
    const std::vector<float> &logEBins,
    const std::string inputPath,
    TFile &outFile,
    const bool verbose);

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