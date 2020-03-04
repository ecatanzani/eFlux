#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <functional>
#include <algorithm>
#include <unistd.h>
#include <numeric>

#include "DmpChain.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtSimuPrimaries.h"

#include "TDirectory.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TH1D.h"

#include "anyoption.h"

struct acceptance_conf
{
    double min_event_energy;
    double max_event_energy;
    double energy_lRatio;
    int shower_axis_delta;
    double vertex_radius;
};

#define _memType "graph"
//#define _memType "histo"

extern void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic,
    const std::string outputPath,
    AnyOption &opt);

extern void buildAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const std::vector<float> &logEBins,
    TFile &outFile);

extern void load_acceptance_struct(acceptance_conf &acceptance_cuts);

extern std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string accInputPath,
    const bool verbose);

extern std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string accInputPath,
    const bool verbose);

extern bool maxElater_cut(
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf &acceptance_cuts,
    const double bgoTotalE);

extern bool maxBarLayer_cut(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    const int nBgoHits);

extern bool BGOTrackContainment_cut(
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf &acceptance_cuts,
    bool passEvent);

extern TF1 readAcceptance( 
    TFile &outFile, 
    const bool verbose,
    const std::string accInputPath);

extern void computeAcceptance(
    const std::string accInputPath,
    const bool verbose,
    const bool pedantic);

#endif