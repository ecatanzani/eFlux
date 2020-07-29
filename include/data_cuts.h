#ifndef DATA_CUTS_H
#define DATA_CUTS_H

#include <memory>
#include <vector>
#include <numeric>

#include "DAMPE_geo_structure.h"

#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"

// DAMPESW includes
#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtPsdHits.h"
#include "DmpStkTrack.h"
#include "DmpSvcPsdEposCor.h"
#include "DmpVSvc.h"

struct cuts_conf
{
	double min_event_energy;
	double max_event_energy;
	double energy_lRatio;
	int shower_axis_delta;
	double vertex_radius;
	int max_rms_shower_width;
	double layer_min_energy;
	int track_X_clusters;
	int track_Y_clusters;
	int track_missingHit_X;
	int track_missingHit_Y;
	int STK_BGO_delta_track;
	int STK_BGO_delta_position;
	double xtrl;
	int STK_PSD_delta_position;
	double PSD_bar_min_energy_release;
	double PSD_charge;
	double PSD_charge_sum;
	double STK_charge;
};

struct data_active_cuts
{
	bool nBarLayer13 = false;
	bool maxRms = false;
	bool track_selection = false;
	bool psd_stk_match = false;
	bool psd_charge = false;
	bool stk_charge = false;
	bool xtrl = false;
	unsigned int nActiveCuts = 0;
};

struct event_filter
{
	bool geometric = false;
	bool BGO_fiducial = false;
	bool all_cut = false;
	bool all_cut_no_xtrl = false;
	bool BGO_fiducial_maxElayer_cut = false;
	bool BGO_fiducial_maxBarLayer_cut = false;
	bool BGO_fiducial_BGOTrackContainment_cut = false;
	bool nBarLayer13_cut = false;
	bool maxRms_cut = false;
	bool track_selection_cut = false;
	bool psd_stk_match_cut = false;
	bool psd_charge_cut = false;
	bool stk_charge_cut = false;
	bool xtrl_cut = false;
	bool psd_charge_measurement = false;
	bool stk_charge_measurement = false;
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
	DmpStkTrack myBestTrack;
};

struct psd_cluster_match
{
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

	psd_cluster_match()	:	icloPsdClu (DAMPE_psd_nLayers, -1), 
    						dxCloPsdClu (DAMPE_psd_nLayers, 9999),
    						icloPsdCluMaxHit (DAMPE_psd_nLayers, -1),
    						dxCloPsdCluMaxHit (DAMPE_psd_nLayers, 9999),
    						icloPsdClu_bgoRec (DAMPE_psd_nLayers, -1),
    						dxCloPsdClu_bgoRec (DAMPE_psd_nLayers, 9999),
    						icloPsdClu_track (DAMPE_psd_nLayers, -1),
    						dxCloPsdClu_track (DAMPE_psd_nLayers, 9999),
    						icloPsdClu2_track (DAMPE_psd_nLayers, -1),
    						dxCloPsdClu2_track (DAMPE_psd_nLayers, 9999)
							{

							}
};

struct stk_charge
{
	double chargeX = -999;
	double chargeY = -999;
};
struct psd_charge
{
	double chargeX = -999;
	double chargeY = -999;
};

extern void print_filter_status(data_active_cuts active_cuts);

extern bool checkBGOreco_data(const std::shared_ptr<DmpEvtBgoRec> bgorec);
extern bool geometric_cut_data(const std::shared_ptr<DmpEvtBgoRec> bgorec);

extern void evaluateTopBottomPosition_data(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	TH1D &h_BGOrec_slopeX,
	TH1D &h_BGOrec_slopeY,
	TH1D &h_BGOrec_interceptX,
	TH1D &h_BGOrec_interceptY,
	TH2D &h_BGOreco_topMap,
	TH2D &h_BGOreco_bottomMap);

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

extern bool psd_stk_match_cut(
	const std::shared_ptr<DmpEvtBgoRec> bgorec,
	const cuts_conf acceptance_cuts,
	const best_track data_cuts,
	psd_cluster_match &clu_matching,
	const std::vector<std::vector<short>> psdCluster_idxBeg,
	const std::vector<std::vector<double>> psdCluster_Z,
	const std::vector<std::vector<double>> psdCluster_maxEcoordinate);

extern bool psd_charge_cut(
	const best_track track,
	const psd_cluster_match clu_matching,
	const std::vector<std::vector<double>> psdCluster_maxE,
	const std::vector<std::vector<short>> psdCluster_idxMaxE,
	const std::vector<double> hitZ,
	const std::vector<short> globalBarID,
	const cuts_conf data_cuts,
	psd_charge &extracted_psd_charge);

extern bool stk_charge_cut(
	const best_track track,
	const std::shared_ptr<TClonesArray> stkclusters,
	const cuts_conf data_cuts,
	stk_charge &extracted_stk_charge);

extern bool xtrl_cut(
	const double sumRms,
	const double lastFracLayer,
	const cuts_conf data_cuts);

#endif