#ifndef ACCEPTANCE_CUTS_H
#define ACCEPTANCE_CUTS_H

#include <memory>
#include <vector>

#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"

#include "acceptance.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtSimuPrimaries.h"
#include "DmpStkTrackHelper.h"

struct best_track
{
    unsigned int n_points = 0;
    std::vector<unsigned int> n_holes = {0, 0};
    std::vector<double> track_slope = {-999, -999};
    std::vector<double> track_intercept = {-999, -999};
    TVector3 track_direction;
    double extr_BGO_topX = -999;
    double extr_BGO_topY = -999;
    double STK_BGO_topX_distance = -999;
    double STK_BGO_topY_distance = -999;
    double angular_distance_STK_BGO = -999;
};

extern bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);

extern bool maxElayer_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const double bgoTotalE);

extern void evaluateEnergyRatio(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const double bgoTotalE,
    TH1D &h_layer_energy_ratio);

extern bool maxBarLayer_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<int> iMaxLayer,
    const std::vector<int> idxBarMaxLayer);

extern bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts);

extern bool nBarLayer13_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::vector<short> layerBarNumber,
    const double bgoTotalE);

extern bool maxRms_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<double> rmsLayer,
    const double bgoTotalE,
    const acceptance_conf acceptance_cuts);

extern bool track_selection_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const acceptance_conf acceptance_cuts,
    best_track &event_best_track);

extern bool xtrl_cut(
    const double sumRms,
    const std::vector<double> fracLayer,
    const acceptance_conf acceptance_cuts);

extern bool psd_charge_cut(
    const std::shared_ptr<DmpEvtPsdHits> psdhits,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const best_track event_best_track);

// Analysis functions
extern void evaluateTopPosition(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    TH1D &h_BGOrec_topX_vs_realX,
    TH1D &h_BGOrec_topY_vs_realY,
    TH1D &h_real_slopeX,
    TH1D &h_real_slopeY,
    TH1D &h_BGOrec_slopeX,
    TH1D &h_BGOrec_slopeY,
    TH1D &h_real_interceptX,
    TH1D &h_real_interceptY,
    TH1D &h_BGOrec_interceptX,
    TH1D &h_BGOrec_interceptY,
    TH2D &h_real_topMap,
    TH2D &h_BGOreco_topMap);

extern bool checkBGOreco(const std::shared_ptr<DmpEvtBgoRec> bgorec);

extern void fillExternalMap(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries, 
    TH2D &h_noBGOenergy_real_topMap);

#endif