#ifndef ACCEPTANCE_CUTS_H
#define ACCEPTANCE_CUTS_H

#include <memory>
#include <vector>

#include "TVector3.h"

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

// DAMPE struct
#define DAMPE_bgo_nLayers 14
#define nSTKladders 192
#define BGO_TopZ 46
#define BGO_BottomZ 448
#define BGO_SideXY 301.25

extern bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);

extern bool maxElayer_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts,
    const double bgoTotalE);

extern bool maxBarLayer_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const int nBgoHits);

extern bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const acceptance_conf acceptance_cuts);

extern bool nBarLayer13_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::vector<short> &layerBarNumber,
    const double bgoTotalE);

extern bool maxRms_cut(
    const std::vector<std::vector<short>> &layerBarNumber,
    const std::vector<double> &rmsLayer,
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

#endif