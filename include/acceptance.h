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
    int max_rms_shower_width;
    int track_X_clusters;
    int track_Y_clusters;
    int track_missingHit_X;
    int track_missingHit_Y;
    int STK_BGO_delta_track;
    int STK_BGO_delta_position;
};

//#define _memType "graph"
#define _memType "histo"
#define _kStep 10000

extern double wtsydp(
    const float minene,
    const float maxene,
    const float index);
    
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
    acceptance_conf &acceptance_cuts,
    const std::string wd);

extern std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string accInputPath,
    const bool verbose);

extern std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string accInputPath,
    const bool verbose);

#endif