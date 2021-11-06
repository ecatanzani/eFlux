#ifndef CUTS_H
#define CUTS_H

#include <memory>
#include <vector>

#include "histos.h"

#include "TVector3.h"
#include "TClonesArray.h"

#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpStkTrack.h"

#include "Dmp/DmpStruct.h"

struct best_track {
	int n_points = -999;
	std::vector<int> n_holes{-999, -999};
	std::vector<double> track_slope{-999, -999};
	std::vector<double> track_intercept{-999, -999};
	TVector3 track_direction{-999, -999, -999};
	double extr_BGO_topX = -999;
	double extr_BGO_topY = -999;
	double STK_BGO_topX_distance = -999;
	double STK_BGO_topY_distance = -999;
	double angular_distance_STK_BGO = -999;
	DmpStkTrack myBestTrack;
};

struct psd_cluster_match {
	std::vector<int> icloPsdClu;
	std::vector<double> dxCloPsdClu;
	std::vector<int> icloPsdCluMaxHit;
	std::vector<double> dxCloPsdCluMaxHit;
	std::vector<int> icloPsdClu_bgoRec;
	std::vector<double> dxCloPsdClu_bgoRec;
	std::vector<int> icloPsdClu_track;
	std::vector<double> dxCloPsdClu_track;
	std::vector<int> icloPsdClu2_track;
	std::vector<double> dxCloPsdClu2_track;

	psd_cluster_match() : icloPsdClu(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu(DAMPE_psd_nLayers, -999),
						  icloPsdCluMaxHit(DAMPE_psd_nLayers, -999),
						  dxCloPsdCluMaxHit(DAMPE_psd_nLayers, -999),
						  icloPsdClu_bgoRec(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu_bgoRec(DAMPE_psd_nLayers, -999),
						  icloPsdClu_track(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu_track(DAMPE_psd_nLayers, -999),
						  icloPsdClu2_track(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu2_track(DAMPE_psd_nLayers, -999)
	{}
};

extern const bool maxElayer_cut(
    const std::vector<double> layer_energies, 
    const double max_energy_lRatio, 
    const double bgoTotalE);

extern const bool maxBarLayer_cut(
	const std::vector<std::vector<short>> layerBarNumber,
	const std::vector<int> iMaxLayer,
	const std::vector<int> idxBarMaxLayer);

extern const bool BGOTrackContainment_cut(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const double shower_axis_delta);

extern const bool nBarLayer13_cut(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::vector<short> layerBarNumber, 
    const double bgoTotalE);

extern const bool maxRms_cut(
	const std::vector<double> layer_energy,
	const std::vector<double> rmsLayer,
	const double bgoTotalE,
	const double max_rms_shower_width);

extern const bool track_selection_cut(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtBgoHits> bgohits,
	const std::shared_ptr<TClonesArray> stkclusters,
	const std::shared_ptr<TClonesArray> stktracks,
    best_track &event_best_track,
	const double STK_BGO_delta_position,
    const double STK_BGO_delta_track,
    const int track_X_clusters,
    const int track_Y_clusters,
	const int track_X_holes,
    const int track_Y_holes);

extern const bool psd_stk_match_cut(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::vector<std::vector<short>> psdCluster_idxBeg,
	const std::vector<std::vector<double>> psdCluster_Z,
	const std::vector<std::vector<double>> psdCluster_maxEcoordinate,
    const best_track &event_best_track,
    psd_cluster_match &clu_matching,
    const double STK_PSD_delta_position);

extern const bool psd_charge_cut(
	const std::vector<std::vector<double>> psdCluster_maxE,
	best_track &event_best_track,
	psd_cluster_match &clu_matching,
	const double PSD_sharge_sum,
	const double PSD_single_charge);

extern const bool stk_charge_cut(
	const std::shared_ptr<TClonesArray> stkclusters,
	best_track &event_best_track,
	const double STK_single_charge);


#endif