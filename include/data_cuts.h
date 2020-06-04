#ifndef DATA_CUTS_H
#define DATA_CUTS_H

#include <memory>
#include <vector>

#include "DAMPE_geo_structure.h"

#include "TClonesArray.h"

// DAMPESW includes
#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdHits.h"
#include "DmpStkTrackHelper.h"


struct cuts_conf
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

struct data_active_cuts
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

extern bool maxElayer_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf acceptance_cuts,
    const double bgoTotalE);

extern bool maxBarLayer_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<int> iMaxLayer,
    const std::vector<int> idxBarMaxLayer);

extern bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts);

extern bool BGOTrackContainment_top_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts);

extern bool nBarLayer13_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::vector<short> layerBarNumber,
    const double bgoTotalE);

extern bool maxRms_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<double> rmsLayer,
    const double bgoTotalE,
    const cuts_conf data_cuts);

extern bool track_selection_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const cuts_conf data_cuts,
    best_track &event_best_track);


extern bool xtrl_cut(
    const double sumRms,
    const std::vector<double> fracLayer,
    const cuts_conf data_cuts);

extern bool psd_charge_cut(
    const std::shared_ptr<DmpEvtPsdHits> psdhits,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf acceptance_cuts,
    const best_track data_cuts);

#endif