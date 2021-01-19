#ifndef DMPFILTERCONTAINER_H
#define DMPFILTERCONTAINER_H

#include <vector>
#include <memory>

#include "DAMPE_geo_structure.h"

#include "DmpVSvc.h"
#include "DmpStkTrack.h"
#include "DmpStkTrack.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtPsdHits.h"
#include "DmpFilterOrbit.h"
#include "DmpSvcPsdEposCor.h"
#include "DmpStkTrackHelper.h"
#include "DmpEvtSimuPrimaries.h"

#include "TChain.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include "config.h"
#include "DmpBgoContainer.h"
#include "DmpPsdContainer.h"

struct best_track
{
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

struct bgo_classifiers
{
	double xtr = -999;
	double xtrl = -999;
};

struct trigger_info
{
	bool mip1 = false;
	bool mip2 = false;
	bool MIP = false;
	bool HET = false;
	bool LET = false;
	bool general = false;
};

struct filter_output
{
	bool out_energy_range = false;
	bool geometric_before_trigger = false;
	bool evt_in_saa = false;
	bool trigger_check = false;
	bool evt_triggered = false;
	bool correct_bgo_reco = false;
	bool good_event = false;
	bool geometric = false;
	bool BGO_fiducial = false;
	bool BGO_fiducial_maxElayer_cut = false;
	bool BGO_fiducial_maxBarLayer_cut = false;
	bool BGO_fiducial_BGOTrackContainment_cut = false;
	bool nBarLayer13_cut = false;
	bool maxRms_cut = false;
	bool track_selection_cut = false;
	bool psd_stk_match_cut = false;
	bool psd_charge_cut = false;
	bool stk_charge_cut = false;
	bool psd_charge_measurement = false;
	bool stk_charge_measurement = false;
	bool all_cut = false;
	bool xtrl_tight_cut = false;
	bool xtrl_loose_cut = false;
};

struct statistics
{
	unsigned int event_counter = 0;
	unsigned int events_in_range = 0;
	unsigned int events_out_range = 0;
	unsigned int triggered_events = 0;
	unsigned int selected_events = 0;
	unsigned int events_in_saa = 0;
};

struct data_evt_time
{
	int start_second = 0;
	int end_second = 0;
	short start_msecond = 0;
	short end_msecond = 0;
	int second = 0;
	short msecond = 0;
};

class DmpFilterContainer
{
public:
	DmpFilterContainer(){};
	~DmpFilterContainer(){};
	void Pipeline(
		const std::shared_ptr<DmpEvtBgoRec> &bgorec,
		const std::shared_ptr<DmpEvtBgoHits> &bgohits,
		const cuts_conf &cuts,
		const double bgoTotalE,
		const double bgoTotalE_corr,
		DmpBgoContainer &bgoVault,
		DmpPsdContainer &psdVault,
		const std::shared_ptr<TClonesArray> &stkclusters,
		const std::shared_ptr<TClonesArray> &stktracks,
		const active_cuts &acuts);
	void Reset();
	void CheckGeometry(
		const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::shared_ptr<DmpEvtSimuPrimaries>(nullptr),
		const std::vector<double> bgoRec_slope = std::vector<double> (2, -999),
		const std::vector<double> bgoRec_intercept = std::vector<double> (2, -999));
	const bool CheckIncomingEvent(
		const std::shared_ptr<DmpEvtHeader> evt_header,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::shared_ptr<DmpEvtSimuPrimaries>(nullptr));
	void EnergyCheck(
		const cuts_conf &cuts,
		const double bgoTotalE_corr,
		const double min_energy,
		const double max_energy);
	void UpdateEvtCounter();
	const filter_output GetFilterOutput();
	const psd_charge GetPSDCharge();
	const stk_charge GetSTKCharge();
	const best_track GetBestTrack();
	DmpStkTrack GetBestTrackObj();
	const bgo_classifiers GetClassifiers();
	const data_evt_time GetDataTime();
	const trigger_info GetTrigger();
	const unsigned int GetStatEvtCounter();
	const unsigned int GetStatEvtInRange();
	const unsigned int GetStatEvtOutRange();
	const unsigned int GetStatiEvtTrigger();
	const unsigned int GetStatEvtSelection();
	const unsigned int GetStatEvtSAA();
	void UpdateEvtTime(
		const std::shared_ptr<TChain> dmpch,
		const std::shared_ptr<DmpEvtHeader> evt_header);
	void PrintDataInfo(
		const std::shared_ptr<TChain> dmpch,
		const bool mc = false);
	void SAACheck(
		const std::shared_ptr<DmpEvtHeader> evt_header,
		const std::shared_ptr<DmpFilterOrbit> pFilter);

private:
	const bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);
	const bool geometric_cut_data(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);
	void check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header);
	void checkBGOreco(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::shared_ptr<DmpEvtSimuPrimaries>(nullptr));
	const bool maxElayer_cut(
		const std::vector<double> layer_energies,
		const cuts_conf acceptance_cuts,
		const double bgoTotalE);
	const bool maxBarLayer_cut(
		const std::vector<std::vector<short>> layerBarNumber,
		const std::vector<int> iMaxLayer,
		const std::vector<int> idxBarMaxLayer);
	const bool BGOTrackContainment_cut(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const cuts_conf data_cuts);
	const bool nBarLayer13_cut(
		const std::shared_ptr<DmpEvtBgoHits> bgohits,
		const std::vector<short> layerBarNumber,
		const double bgoTotalE);
	const bool maxRms_cut(
		const std::vector<std::vector<short>> layerBarNumber,
		const std::vector<double> rmsLayer,
		const double bgoTotalE,
		const cuts_conf data_cuts);
	const bool track_selection_cut(
		const std::shared_ptr<DmpEvtBgoRec> bgorec,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const std::shared_ptr<DmpEvtBgoHits> bgohits,
		const std::shared_ptr<TClonesArray> stkclusters,
		const std::shared_ptr<TClonesArray> stktracks,
		const cuts_conf data_cuts);
	const bool psd_stk_match_cut(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const cuts_conf acceptance_cuts,
		const std::vector<std::vector<short>> psdCluster_idxBeg,
		const std::vector<std::vector<double>> psdCluster_Z,
		const std::vector<std::vector<double>> psdCluster_maxEcoordinate);
	const bool psd_charge_cut(
		const std::vector<std::vector<double>> psdCluster_maxE,
		const std::vector<std::vector<short>> psdCluster_idxMaxE,
		const std::vector<double> hitZ,
		const std::vector<short> globalBarID,
		const cuts_conf data_cuts);
	const bool stk_charge_cut(
		const std::shared_ptr<TClonesArray> stkclusters,
		const cuts_conf data_cuts);
	const bool xtrl_tight_cut(
		const double input_xtrl,
		const double cut_value = 8.5);
	const bool xtrl_loose_cut(
		const double input_xtrl,
		const double energy);

	void reset_stk_best_track();
	void reset_psd_clusters();
	void reset_psd_charges();
	void reset_stk_charges();
	void reset_classifiers();
	void reset_filter_output();
	void reset_time();
	void reset_trigger();

	best_track event_best_track;
	psd_cluster_match clu_matching;
	psd_charge extracted_psd_charge;
	stk_charge extracted_stk_charge;
	bgo_classifiers classifier;
	filter_output output;
	statistics particle_counter;
	data_evt_time time;
	trigger_info evt_trigger;
};

#endif