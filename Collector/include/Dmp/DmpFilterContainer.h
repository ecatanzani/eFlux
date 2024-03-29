#ifndef DMPFILTERCONTAINER_H
#define DMPFILTERCONTAINER_H

#include <tuple>
#include <vector>
#include <memory>

#include "Dmp/DmpGeoStruct.h"

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

#include "TF1.h"
#include "TChain.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include "config.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpPsdContainer.h"

struct best_track
{
	int n_points 								{-999};
	std::vector<int> n_holes					{-999, -999};
	std::vector<double> track_slope				{-999, -999};
	std::vector<double> track_intercept			{-999, -999};
	TVector3 track_direction					{-999, -999, -999};
	double extr_BGO_topX 						{-999};
	double extr_BGO_topY 						{-999};
	double STK_BGO_topX_distance 				{-999};
	double STK_BGO_topY_distance 				{-999};
	double angular_distance_STK_BGO 			{-999};
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
	std::vector<int> icloPsdClu_track_fiducial;
	std::vector<double> dxCloPsdClu_track;
	std::vector<double> dxCloPsdClu_track_fiducial;

	bool X_match {false};
	bool Y_match {false};

	psd_cluster_match() : icloPsdClu					(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu					(DAMPE_psd_nLayers, -999),
						  icloPsdCluMaxHit				(DAMPE_psd_nLayers, -999),
						  dxCloPsdCluMaxHit				(DAMPE_psd_nLayers, -999),
						  icloPsdClu_bgoRec				(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu_bgoRec			(DAMPE_psd_nLayers, -999),
						  icloPsdClu_track				(DAMPE_psd_nLayers, -999),
						  icloPsdClu_track_fiducial		(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu_track				(DAMPE_psd_nLayers, -999),
						  dxCloPsdClu_track_fiducial 	(DAMPE_psd_nLayers, -999)
						{
						}
};

struct stk_charge
{
	double chargeX {-999};
	double chargeY {-999};
};
struct psd_charge
{
	double chargeX {-999};
	double chargeY {-999};
};

struct bgo_classifiers
{
	double xtr 	{-999};
	double xtrl {-999};
};

struct trigger_info
{
	bool mip1 			{false};
	bool mip2 			{false};
	bool MIP 			{false};
	bool HET 			{false};
	bool LET 			{false};
	bool general 		{false};
	bool unbiased 		{false};
};

struct filter_output
{
	bool out_energy_range 							{false};
	bool geometric_before_trigger 					{false};
	bool evt_in_saa 								{false};
	bool trigger_check 								{false};
	bool evt_triggered 								{false};
	bool correct_bgo_reco 							{false};
	bool good_event 								{false};
	bool geometric 									{false};
	bool BGO_fiducial 								{false};
	bool BGO_fiducial_HET 							{false};
	bool BGO_fiducial_maxElayer_cut 				{false};
	bool BGO_fiducial_maxBarLayer_cut 				{false};
	bool BGO_fiducial_BGOTrackContainment_cut 		{false};
	bool nBarLayer13_cut 							{false};
	bool maxRms_cut 								{false};
	bool sumrms_low_energy_cut						{false};
	bool stk_fiducial_volume						{false};
	bool stk_fiducial_volume_X						{false};
	bool stk_fiducial_volume_Y						{false};
	bool track_selection_cut 						{false};
	bool track_selection_cut_no_3hit_recover 		{false};
	bool three_cluster_only_track					{false};
	bool stk_1rm_cut								{false};
	bool psd_stk_match_cut 							{false};
	bool psd_fiducial_volume 						{false};
	bool psd_fiducial_volume_X 						{false};
	bool psd_fiducial_volume_Y 						{false};
	bool psd_stk_match_cut_x 						{false};
	bool psd_stk_match_cut_y 						{false};
	bool psd_charge_cut 							{false};
	bool psd_charge_cut_no_single_view_recover		{false};
	bool stk_charge_cut 							{false};
	bool psd_charge_measurement 					{false};
	bool stk_charge_measurement 					{false};
	bool all_cut 									{false};
	bool xtrl_tight_cut 							{false};
	bool xtrl_loose_cut 							{false};
};

struct statistics
{
	unsigned int event_counter 			{0};
	unsigned int events_in_range 		{0};
	unsigned int events_out_range 		{0};
	unsigned int triggered_events 		{0};
	unsigned int selected_events 		{0};
	unsigned int events_in_saa 			{0};
};

struct data_evt_time
{
	int start_second 			{0};
	int end_second 				{0};
	short start_msecond 		{0};
	short end_msecond 			{0};
	int second 					{0};
	short msecond 				{0};
};

struct stk_correction_functions
{
	std::shared_ptr<TF1> f_stk_correction_20_100 = std::shared_ptr<TF1>(nullptr);
	std::shared_ptr<TF1> f_stk_correction_100_250 = std::shared_ptr<TF1>(nullptr);
	std::shared_ptr<TF1> f_stk_correction_250_500 = std::shared_ptr<TF1>(nullptr);
	std::shared_ptr<TF1> f_stk_correction_500_1000 = std::shared_ptr<TF1>(nullptr);
	std::shared_ptr<TF1> f_stk_correction_1000_3000 = std::shared_ptr<TF1>(nullptr);
	std::shared_ptr<TF1> f_stk_correction_3000 = std::shared_ptr<TF1>(nullptr);
};

class DmpFilterContainer
{
public:
	DmpFilterContainer(){};
	~DmpFilterContainer(){};
	void LoadStkCorrectionFunctions(std::string path);
	void Pipeline(
		const std::shared_ptr<DmpEvtBgoRec> &bgorec,
		const std::shared_ptr<DmpEvtBgoHits> &bgohits,
		const cuts_conf &cuts,
		const double bgoTotalE,
		const double bgoTotalE_corr,
		DmpBgoContainer &bgoVault,
		DmpStkContainer &stkVault,
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
	//const bool CheckIncomingEventNoTrigger();
	void EnergyCheck(
		const cuts_conf &cuts,
		const double bgoTotalE_corr,
		const double min_energy,
		const double max_energy);
	void UpdateEvtCounter();
	const filter_output GetFilterOutput();
	const psd_charge GetPSDCharge();
	const std::tuple<double, double, double, double> GetPSDSTKMatchDistances();
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

protected:
	const bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries);
	const bool geometric_cut_data(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);
	void check_trigger(
		const std::shared_ptr<DmpEvtHeader> evt_header,
		const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::shared_ptr<DmpEvtSimuPrimaries>(nullptr));
	void checkBGOreco(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);

	// Analysis cuts
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
		const std::vector<short> layerBarIndex,
		const double bgoTotalE);
	const bool maxRms_cut(
		const std::vector<double> layer_energy,
		const std::vector<double> rmsLayer,
		const double bgoTotalE,
		const cuts_conf data_cuts);

	
	const bool sumrms_low_energy_cut(
		const double bgoTotalE,
		const double sumrms,
		const double bgo_direction_cosine);
	const bool stk_1rm_cut(
		const std::vector<double> stkEcore1Rm,
		const std::vector<unsigned int> nStkClu1Rm,
		const double bgoTotalE_corr);
	const bool rvalue_cut(const double rvalue);
	const bool lvalue_cut(const double lvalue, const double bgoTotalE, const double bgoTotalE_corr);

	void link_ladders(std::vector<int> &LadderToLayer);
	void fill_BGO_vectors(
		TVector3 &bgoRecEntrance,
		TVector3 &bgoRecDirection,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept);
	void get_track_points(
		DmpStkTrack *track,
		const std::shared_ptr<TClonesArray> stkclusters,
		const std::vector<int> LadderToLayer,
		std::vector<int> &track_nHoles,
		best_track &event_best_track,
		const bool best_track = false);
	const bool track_selection_cut(
		const std::shared_ptr<DmpEvtBgoRec> bgorec,
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const std::shared_ptr<DmpEvtBgoHits> bgohits,
		const std::shared_ptr<TClonesArray> stkclusters,
		const std::shared_ptr<TClonesArray> stktracks,
		const cuts_conf data_cuts,
		const bool recover_3_hits_tracks = true,
		const bool update_struct = true);
	const bool psd_fiducial_volume_cut();
	const bool psd_stk_match_cut(
		const std::vector<double> bgoRec_slope,
		const std::vector<double> bgoRec_intercept,
		const cuts_conf acceptance_cuts,
		const std::vector<std::vector<short>> psdCluster_idxBeg,
		const std::vector<std::vector<double>> psdCluster_Z,
		const std::vector<std::vector<double>> psdCluster_maxEcoordinate);
	void psd_charge_measurement(
		const std::vector<std::vector<double>> psdCluster_maxE,
		const std::vector<std::vector<short>> psdCluster_idxMaxE,
		const std::vector<double> hitZ,
		const std::vector<short> globalBarID,
		const bool update_struct = true);
	void stk_charge_measurement(
		const std::shared_ptr<TClonesArray> stkclusters,
		const bool update_struct = true);
	const bool psd_charge_cut(
		const cuts_conf data_cuts,
		const bool recover_one_view = true);
	const bool stk_charge_cut(const double hi_cut, const double med_cut, const double low_cut, const bool reject_if_no_clusters = false);
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

	stk_correction_functions stk_cleaning_functions;
};

#endif