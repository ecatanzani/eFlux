#ifndef TUPLE_H
#define TUPLE_H

#include "config.h"
#include "Dmp/DmpGeoStruct.h"
#include "Dmp/DmpFilterContainer.h"
#include "Efficiency/efficiency.h"
#include "Preselection/preselection.h"

#include "TTree.h"
#include "TVector3.h"

#include <memory>
#include <vector>

class ntuple
{
public:
	ntuple() : eLayer				(DAMPE_bgo_nLayers, -999),
			   layerBarEnergy		(DAMPE_bgo_nLayers, std::vector<double>(DAMPE_bgo_bars_layer, -999)),
			   rmsLayer				(DAMPE_bgo_nLayers, -999),
			   fracLayer			(DAMPE_bgo_nLayers, -999),
			   nud_adc				(DAMPE_NUD_channels, -999),
			   STK_plane_clusters	(DAMPE_stk_planes, -999)
			   {
			   };
	~ntuple(){};
	void Write(TFile &outfile);

protected:
	void set_active_cuts(const active_cuts &acuts);
	void fill_trigger_info(const trigger_info &evt_trigger);
	void fill_psd_info(
		const double psdstkmatch_x,
		const double psdstkmatch_y,
		const double psdstkmatch_x_fiducial,
		const double psdstkmatch_y_fiducial);
	void fill_stk_info(
		const best_track &event_best_track,
		const std::vector<int> _stk_clusters_on_plane);
	void fill_bgo_info(
		const double raw_energy,
		const double corr_energy,
		const std::vector<double> &energy_release_layer,
		const std::vector<std::vector<double>> &energy_release_layer_bar,
		const std::vector<double> &bgoRec_slope,
		const std::vector<double> &bgoRec_intercept,
		const TVector3 &bgoRec_trajectory2D,
		const double sumRMS,
		const std::vector<double> &rms_layer,
		const std::vector<double> &bgo_fracLayer,
		const double lastFracLayer,
		const double frac_layer_13,
		const int last_bgo_layer,
		const int bgo_entries);
	void fill_psdcharge_info(const psd_charge &extracted_psd_charge);
	void fill_stkcharge_info(const stk_charge &extracted_stk_charge);
	void fill_classifier_info(const bgo_classifiers &classifier);
	void fill_nud_info(
		const std::vector<int> adc,
		const int total_adc,
		const int max_adc,
		const int max_channel_id);
	void fill_preselectionfilters_info(const presel_output &preselect);
	void fill_preselectionfiltersefficiency_info(const eff_output &eff_preselect);
	void core_reset();

	// Tree
	std::unique_ptr<TTree> DmpNtupTree;
	// Trigger
	bool mip1_trigger 														{false};
	bool mip2_trigger 														{false};
	bool HET_trigger 														{false};
	bool LET_trigger 														{false};
	bool MIP_trigger 														{false};
	bool general_trigger 													{false};
	bool unbiased_trigger 													{false};
	// STK
	int STK_bestTrack_npoints 												{-999};
	int STK_bestTrack_nholesX 												{-999};
	int STK_bestTrack_nholesY 												{-999};
	double STK_bestTrack_slopeX 											{-999};
	double STK_bestTrack_slopeY 											{-999};
	double STK_bestTrack_interceptX 										{-999};
	double STK_bestTrack_interceptY 										{-999};
	double STK_bestTrack_costheta 											{-999};
	double STK_bestTrack_phi 												{-999};
	double STK_bestTrack_extr_BGO_topX 										{-999};
	double STK_bestTrack_extr_BGO_topY 										{-999};
	double STK_bestTrack_STK_BGO_topX_distance 								{-999};
	double STK_bestTrack_STK_BGO_topY_distance 								{-999};
	double STK_bestTrack_angular_distance_STK_BGO 							{-999};
	double STK_chargeX 														{-999};
	double STK_chargeY 														{-999};
	double STK_charge 														{-999};
	std::vector<int> STK_plane_clusters;
	// BGO
	double energy 															{-999};
	double energy_corr 														{-999};
	std::vector<double> eLayer;
	std::vector<std::vector<double>> layerBarEnergy;
	double BGOrec_slopeX 													{-999};
	double BGOrec_slopeY 													{-999};
	double BGOrec_interceptX 												{-999};
	double BGOrec_interceptY 												{-999};
	TVector3 trajectoryDirection2D											{-999, -999, -999};
	double sumRms 															{-999};
	std::vector<double> rmsLayer;
	std::vector<double> fracLayer;
	double fracLast 														{-999};
	double fracLast_13 														{-999};
	int lastBGOLayer 														{-999};
	int nBGOentries 														{-999};
	// PSD
	double PSD_chargeX 														{-999};
	double PSD_chargeY 														{-999};
	double PSD_charge 														{-999};
	double SPD_STK_match_X_distance											{-999};
	double SPD_STK_match_Y_distance											{-999};
	double SPD_STK_match_X_distance_fiducial_volume							{-999};
	double SPD_STK_match_Y_distance_fiducial_volume							{-999};
	// NUD
	std::vector<int> nud_adc;
	int nud_total_adc														{-999};
	int nud_max_adc 														{-999};
	int nud_max_channel_id 													{-999};
	// Classifiers
	double xtr 																{-999};
	double xtrl 															{-999};
	// Filters
	bool evtfilter_out_energy_range 										{false};
	bool evtfilter_evt_triggered 											{false};
	bool evtfilter_correct_bgo_reco 										{false};
	bool evtfilter_good_event 												{false};
	bool evtfilter_geometric 												{false};
	bool evtfilter_BGO_fiducial 											{false};
	bool evtfilter_BGO_fiducial_HET 										{false};
	bool evtfilter_BGO_fiducial_maxElayer_cut 								{false};
	bool evtfilter_BGO_fiducial_maxBarLayer_cut 							{false};
	bool evtfilter_BGO_fiducial_BGOTrackContainment_cut 					{false};
	bool evtfilter_nBarLayer13_cut 											{false};
	bool evtfilter_maxRms_cut 												{false};
	bool evtfilter_track_selection_cut 										{false};
	bool evtfilter_track_selection_cut_no_3hit_recover 						{false};
	bool evtfilter_psd_stk_match_cut 										{false};
	bool evtfilter_psd_stk_match_cut_X_view 								{false};
	bool evtfilter_psd_stk_match_cut_Y_view 								{false};
	bool evtfilter_psd_stk_match_cut_psd_fiducial_volume					{false};
	bool evtfilter_psd_stk_match_cut_psd_fiducial_volume_X					{false};
	bool evtfilter_psd_stk_match_cut_psd_fiducial_volume_Y					{false};
	bool evtfilter_psd_charge_cut 											{false};
	bool evtfilter_psd_charge_cut_no_single_view_recover					{false};
	bool evtfilter_stk_charge_cut 											{false};
	bool evtfilter_psd_charge_measurement 									{false};
	bool evtfilter_stk_charge_measurement 									{false};
	bool evtfilter_xtrl_tight_cut 											{false};
	bool evtfilter_xtrl_loose_cut 											{false};
	bool evtfilter_all_cut 													{false};
	// Preselection Cuts
	bool cut_nBarLayer13 													{false};
	bool cut_maxRms 														{false};
	bool cut_track_selection 												{false};
	bool cut_psd_stk_match 													{false};
	bool cut_psd_charge 													{false};
	bool cut_stk_charge 													{false};
	int nActiveCuts 														{-999};
	// Preselection Filters
	bool preselection_maxelayer_cut 										{false};
    bool preselection_maxbarlayer_cut 										{false};
    bool preselection_bgotrack_cut 											{false};
    bool preselection_bgofiducial_cut 										{false};
    bool preselection_nbarlayer13_cut 										{false};
    bool preselection_maxrms_cut 											{false};
    bool preselection_trackselection_cut 									{false};
    bool preselection_psdstkmatch_cut 										{false};
    bool preselection_psdcharge_cut 										{false};
    bool preselection_stkcharge_cut 										{false};
	bool preselection_maxelayer_lastcut										{false};
	bool preselection_maxbarlayer_lastcut									{false};
	bool preselection_bgotrack_lastcut										{false};
	bool preselection_bgofiducial_lastcut									{false};
	bool preselection_nbarlayer13_lastcut									{false};
	bool preselection_maxrms_lastcut										{false};
	bool preselection_trackselection_lastcut								{false};
	bool preselection_psdstkmatch_lastcut									{false};
	// Efficiency Preselection Filters
	bool trigger_efficiency_preselection                        			{false};
	bool trigger_efficiency_preselection_is_het                 			{false};
	bool trigger_efficiency_preselection_is_let                 			{false};
	bool trigger_efficiency_preselection_is_unb                 			{false};
	bool maxrms_efficiency_preselection                         			{false};
	bool maxrms_efficiency_preselection_accepted                			{false};
	bool nbarlayer13_efficiency_preselection                    			{false};
	bool nbarlayer13_efficiency_preselection_accepted           			{false};
	bool maxrms_and_nbarlayer13_efficiency_preselection						{false};
	bool maxrms_and_nbarlayer13_efficiency_preselection_accepted 			{false};
	bool track_efficiency_preselection                          			{false};
	bool track_efficiency_preselection_accepted                 			{false};
	bool psdstkmatch_efficiency_preselection                    			{false};
	bool psdstkmatch_efficiency_preselection_accepted           			{false};
	bool psdcharge_efficiency_preselection                      			{false};
	bool psdcharge_efficiency_preselection_accepted             			{false};
};

#endif