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

#include "DmpChain.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtSimuPrimaries.h"

#include "TDirectory.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TGraph.h"

#define GetCurrentDir getcwd
struct acceptance_conf
{
    double event_energy;
    double energy_lRatio;
    int shower_axis_delta;
    double vertex_radius;
};

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

#endif