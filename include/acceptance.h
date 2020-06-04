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

#include "TClonesArray.h"
#include "TString.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"

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

#include "anyoption.h"

// DAMPE struct
const int DAMPE_bgo_nLayers = 14;
const int nSTKladders = 192;
const double BGO_TopZ = 46;
const double BGO_BottomZ = 448;
const double BGO_SideXY = 301.25;

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
    double xtrl;
    int STK_PSD_delta_position;
    double PSD_bar_min_energy_release;
};

struct acceptance_active_cuts
{
    bool geometry = false;
    bool BGO_fiducial = false;
    bool nBarLayer13 = false;
    bool maxRms = false;
    bool track_selection = false;
    bool xtrl = false;
    bool psd_charge = false;
    unsigned int nActiveCuts = 0;
};

//#define _memType "graph"
//#define _memType "histo"

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
    acceptance_active_cuts &active_cuts,
    const std::string wd);

#endif