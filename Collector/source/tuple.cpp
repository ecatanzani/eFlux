#include "tuple.h"

void ntuple::set_active_cuts(const active_cuts &acuts)
{
	cut_nBarLayer13 = acuts.nBarLayer13;
	cut_maxRms = acuts.maxRms;
	cut_track_selection = acuts.track_selection;
	cut_psd_stk_match = acuts.psd_stk_match;
	cut_psd_charge = acuts.psd_charge;
	cut_stk_charge = acuts.stk_charge;
	nActiveCuts = (int)acuts.nActiveCuts;
}

void ntuple::fill_trigger_info(const trigger_info &evt_trigger)
{
	mip1_trigger = evt_trigger.mip1;
	mip2_trigger = evt_trigger.mip2;
	HET_trigger = evt_trigger.HET;
	LET_trigger = evt_trigger.LET;
	MIP_trigger = evt_trigger.MIP;
	general_trigger = evt_trigger.general;
}

void ntuple::fill_stk_info(const best_track &event_best_track)
{
	if (evtfilter_track_selection_cut)
	{
		auto track = (event_best_track.myBestTrack);
		STK_bestTrack_npoints = event_best_track.n_points;
		STK_bestTrack_nholesX = event_best_track.n_holes[0];
		STK_bestTrack_nholesY = event_best_track.n_holes[1];
		STK_bestTrack_slopeX = event_best_track.track_slope[0];
		STK_bestTrack_slopeY = event_best_track.track_slope[1];
		STK_bestTrack_interceptX = event_best_track.track_intercept[0];
		STK_bestTrack_interceptY = event_best_track.track_intercept[1];
		STK_bestTrack_costheta = track.getDirection().CosTheta();
		STK_bestTrack_phi = track.getDirection().Phi();
		STK_bestTrack_extr_BGO_topX = event_best_track.extr_BGO_topX;
		STK_bestTrack_extr_BGO_topY = event_best_track.extr_BGO_topY;
		STK_bestTrack_STK_BGO_topX_distance = event_best_track.STK_BGO_topX_distance;
		STK_bestTrack_STK_BGO_topY_distance = event_best_track.STK_BGO_topY_distance;
		STK_bestTrack_angular_distance_STK_BGO = event_best_track.angular_distance_STK_BGO;
	}
}

void ntuple::fill_bgo_info(
	const double raw_energy,
	const double corr_energy,
	const std::vector<double> &energy_release_layer,
	const std::vector<double> &bgoRec_slope,
	const std::vector<double> &bgoRec_intercept,
	const double sumRMS,
	const std::vector<double> rms_layer,
	const std::vector<double> bgo_fracLayer,
	const double lastFracLayer,
	const double frac_layer_13,
	const int last_bgo_layer,
	const int bgo_entries,
	const std::vector<double> energy_1_moliere_radius,
	const std::vector<double> energy_2_moliere_radius,
	const std::vector<double> energy_3_moliere_radius,
	const std::vector<double> energy_5_moliere_radius)
{
	if (evtfilter_correct_bgo_reco)
	{
		energy = raw_energy;
		energy_corr = corr_energy;
		eLayer = energy_release_layer;
		BGOrec_slopeX = bgoRec_slope[0];
		BGOrec_slopeY = bgoRec_slope[1];
		BGOrec_interceptX = bgoRec_intercept[0];
		BGOrec_interceptY = bgoRec_intercept[1];
		sumRms = sumRMS;
		rmsLayer = rms_layer;
		fracLayer = bgo_fracLayer;
		fracLast = lastFracLayer;
		fracLast_13 = frac_layer_13;
		lastBGOLayer = last_bgo_layer;
		nBGOentries = bgo_entries;
		energy_1R_radius = energy_1_moliere_radius;
		energy_2R_radius = energy_2_moliere_radius;
		energy_3R_radius = energy_3_moliere_radius;
		energy_5R_radius = energy_5_moliere_radius;
	}
}

void ntuple::fill_psdcharge_info(const psd_charge &extracted_psd_charge)
{
	if (evtfilter_psd_charge_measurement)
		if (extracted_psd_charge.chargeX != -999 && extracted_psd_charge.chargeY != -999)
		{
			PSD_chargeX = extracted_psd_charge.chargeX;
			PSD_chargeY = extracted_psd_charge.chargeY;
			PSD_charge = 0.5 * (extracted_psd_charge.chargeX + extracted_psd_charge.chargeY);
		}
}

void ntuple::fill_stkcharge_info(const stk_charge &extracted_stk_charge)
{
	if (evtfilter_stk_charge_measurement)
		if (extracted_stk_charge.chargeX != -999 && extracted_stk_charge.chargeY != -999)
		{
			STK_chargeX = extracted_stk_charge.chargeX;
			STK_chargeY = extracted_stk_charge.chargeY;
			STK_charge = 0.5 * (extracted_stk_charge.chargeX + extracted_stk_charge.chargeY);
		}
}

void ntuple::fill_classifier_info(const bgo_classifiers &classifier)
{
	xtr = classifier.xtr;
	xtrl = classifier.xtrl;
}

void ntuple::Write(TFile &outfile)
{
	outfile.cd();
	DmpNtupTree->Write();
}

void ntuple::fill_nud_info(
	const std::vector<double> adc,
	const double total_adc,
	const double max_adc,
	const int max_channel_id)
{
	nud_adc = adc;
	nud_total_adc = total_adc;
	nud_max_adc = max_adc;
	nud_max_channel_id = max_channel_id;
}

void ntuple::core_reset()
{
	// Trigger
	mip1_trigger = false;
	mip2_trigger = false;
	HET_trigger = false;
	LET_trigger = false;
	MIP_trigger = false;
	general_trigger = false;
	// STK
	STK_bestTrack_npoints = -999;
	STK_bestTrack_nholesX = -999;
	STK_bestTrack_nholesY = -999;
	STK_bestTrack_slopeX = -999;
	STK_bestTrack_slopeY = -999;
	STK_bestTrack_interceptX = -999;
	STK_bestTrack_interceptY = -999;
	STK_bestTrack_costheta = -999;
	STK_bestTrack_phi = -999;
	STK_bestTrack_extr_BGO_topX = -999;
	STK_bestTrack_extr_BGO_topY = -999;
	STK_bestTrack_STK_BGO_topX_distance = -999;
	STK_bestTrack_STK_BGO_topY_distance = -999;
	STK_bestTrack_angular_distance_STK_BGO = -999;
	STK_chargeX = -999;
	STK_chargeY = -999;
	STK_charge = -999;
	// BGO
	energy = -999;
	energy_corr = -999;
	eLayer = std::vector<double> (DAMPE_bgo_nLayers, -999);
	BGOrec_slopeX = -999;
	BGOrec_slopeY = -999;
	BGOrec_interceptX = -999;
	BGOrec_interceptY = -999;
	sumRms = -999;
	rmsLayer = std::vector<double> (DAMPE_bgo_nLayers, -999);
	fracLayer = std::vector<double> (DAMPE_bgo_nLayers, -999);
	fracLast = -999;
	fracLast_13 = -999;
	lastBGOLayer = -999;
	nBGOentries = -999;
	energy_1R_radius = std::vector<double> (DAMPE_bgo_nLayers, -999);
	energy_2R_radius = std::vector<double> (DAMPE_bgo_nLayers, -999);
	energy_3R_radius = std::vector<double> (DAMPE_bgo_nLayers, -999);
	energy_5R_radius = std::vector<double> (DAMPE_bgo_nLayers, -999);
	// PSD
	PSD_chargeX = -999;
	PSD_chargeY = -999;
	PSD_charge = -999;
	// NUD
	nud_adc = std::vector<double> (DAMPE_NUD_channels, -999);
	nud_total_adc = -999;
	nud_max_adc = -999;
	nud_max_channel_id = -999;
	// Classifiers
	xtr = -999;
	xtrl = -999;
	// Filters
	evtfilter_out_energy_range = false;
	evtfilter_evt_triggered = false;
	evtfilter_correct_bgo_reco = false;
	evtfilter_good_event = false;
	evtfilter_geometric = false;
	evtfilter_BGO_fiducial = false;
	evtfilter_BGO_fiducial_maxElayer_cut = false;
	evtfilter_BGO_fiducial_maxBarLayer_cut = false;
	evtfilter_BGO_fiducial_BGOTrackContainment_cut = false;
	evtfilter_nBarLayer13_cut = false;
	evtfilter_maxRms_cut = false;
	evtfilter_track_selection_cut = false;
	evtfilter_psd_stk_match_cut = false;
	evtfilter_psd_charge_cut = false;
	evtfilter_stk_charge_cut = false;
	evtfilter_psd_charge_measurement = false;
	evtfilter_stk_charge_measurement = false;
	evtfilter_xtrl_tight_cut = false;
	evtfilter_xtrl_loose_cut = false;
	evtfilter_all_cut = false;
	// Preselection Cuts
	cut_nBarLayer13 = false;
	cut_maxRms = false;
	cut_track_selection = false;
	cut_psd_stk_match = false;
	cut_psd_charge = false;
	cut_stk_charge = false;
	nActiveCuts = -999;
}