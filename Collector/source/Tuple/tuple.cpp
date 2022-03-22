#include "Tuple/tuple.h"

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
	unbiased_trigger = evt_trigger.unbiased;
}

void ntuple::fill_psd_info(
	const double psdstkmatch_x,
	const double psdstkmatch_y,
	const double psdstkmatch_x_fiducial,
	const double psdstkmatch_y_fiducial)
	{
		if (evtfilter_track_selection_cut)
		{
			SPD_STK_match_X_distance = psdstkmatch_x;
			SPD_STK_match_Y_distance = psdstkmatch_y;
			SPD_STK_match_X_distance_fiducial_volume = psdstkmatch_x_fiducial;
			SPD_STK_match_Y_distance_fiducial_volume = psdstkmatch_y_fiducial;
		}
	}

void ntuple::fill_stk_info(
	const best_track &event_best_track,
	const std::vector<int> _stk_clusters_on_plane)
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
	STK_plane_clusters = _stk_clusters_on_plane;
}

void ntuple::fill_bgo_info(
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
	const int bgo_entries)
{
	if (evtfilter_correct_bgo_reco)
	{
		energy = raw_energy;
		energy_corr = corr_energy;
		eLayer = energy_release_layer;
		layerBarEnergy = energy_release_layer_bar;
		BGOrec_slopeX = bgoRec_slope[0];
		BGOrec_slopeY = bgoRec_slope[1];
		BGOrec_interceptX = bgoRec_intercept[0];
		BGOrec_interceptY = bgoRec_intercept[1];
		trajectoryDirection2D = bgoRec_trajectory2D;
		sumRms = sumRMS;
		rmsLayer = rms_layer;
		fracLayer = bgo_fracLayer;
		fracLast = lastFracLayer;
		fracLast_13 = frac_layer_13;
		lastBGOLayer = last_bgo_layer;
		nBGOentries = bgo_entries;
	}
}

void ntuple::fill_psdcharge_info(const psd_charge &extracted_psd_charge)
{
	PSD_chargeX = extracted_psd_charge.chargeX;
	PSD_chargeY = extracted_psd_charge.chargeY;
	PSD_charge = 0.5 * (extracted_psd_charge.chargeX + extracted_psd_charge.chargeY);		
}

void ntuple::fill_stkcharge_info(const stk_charge &extracted_stk_charge)
{
	STK_chargeX = extracted_stk_charge.chargeX;
	STK_chargeY = extracted_stk_charge.chargeY;
	STK_charge = 0.5 * (extracted_stk_charge.chargeX + extracted_stk_charge.chargeY);		
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
	const std::vector<int> adc,
	const int total_adc,
	const int max_adc,
	const int max_channel_id)
{
	if (evtfilter_good_event)
	{
		nud_adc = adc;
		nud_total_adc = total_adc;
		nud_max_adc = max_adc;
		nud_max_channel_id = max_channel_id;
	}
}

void ntuple::fill_preselectionfilters_info(const presel_output &preselect) 
{
	if (evtfilter_good_event) 
	{
		preselection_maxelayer_cut = preselect.maxelayer_cut;
		preselection_maxbarlayer_cut = preselect.maxbarlayer_cut;
		preselection_bgotrack_cut = preselect.bgotrack_cut;
		preselection_bgofiducial_cut = preselect.bgofiducial_cut;
		preselection_nbarlayer13_cut = preselect.nbarlayer13_cut;
		preselection_maxrms_cut = preselect.maxrms_cut;
		preselection_trackselection_cut = preselect.trackselection_cut;
		preselection_psdstkmatch_cut = preselect.psdstkmatch_cut;
		preselection_psdcharge_cut = preselect.psdcharge_cut;
		preselection_stkcharge_cut = preselect.stkcharge_cut;
		preselection_maxelayer_lastcut = preselect.maxelayer_lastcut;
		preselection_maxbarlayer_lastcut = preselect.maxbarlayer_lastcut;
		preselection_bgotrack_lastcut = preselect.bgotrack_lastcut;
		preselection_bgofiducial_lastcut = preselect.bgofiducial_lastcut;
		preselection_nbarlayer13_lastcut = preselect.nbarlayer13_lastcut;
		preselection_maxrms_lastcut = preselect.maxrms_lastcut;
		preselection_trackselection_lastcut = preselect.trackselection_lastcut;
		preselection_psdstkmatch_lastcut = preselect.psdstkmatch_lastcut;
	}
}

void ntuple::fill_preselectionfiltersefficiency_info(const eff_output &eff_preselect)
{
	if (evtfilter_correct_bgo_reco)
	{
		trigger_efficiency_preselection 								= eff_preselect.trigger_efficiency_preselection;
		trigger_efficiency_preselection_is_het 							= eff_preselect.trigger_efficiency_preselection_is_het;
		trigger_efficiency_preselection_is_let 							= eff_preselect.trigger_efficiency_preselection_is_let;
		trigger_efficiency_preselection_is_unb 							= eff_preselect.trigger_efficiency_preselection_is_unb;
		maxrms_efficiency_preselection 									= eff_preselect.maxrms_efficiency_preselection;
		maxrms_efficiency_preselection_accepted 						= eff_preselect.maxrms_efficiency_preselection_accepted;
		nbarlayer13_efficiency_preselection 							= eff_preselect.nbarlayer13_efficiency_preselection;
		nbarlayer13_efficiency_preselection_accepted 					= eff_preselect.nbarlayer13_efficiency_preselection_accepted;
		maxrms_and_nbarlayer13_efficiency_preselection 					= eff_preselect.maxrms_and_nbarlayer13_efficiency_preselection;
		maxrms_and_nbarlayer13_efficiency_preselection_accepted 		= eff_preselect.maxrms_and_nbarlayer13_efficiency_preselection_accepted;
		track_efficiency_preselection 									= eff_preselect.track_efficiency_preselection;
		track_efficiency_preselection_accepted 							= eff_preselect.track_efficiency_preselection_accepted;
		psdstkmatch_efficiency_preselection 							= eff_preselect.psdstkmatch_efficiency_preselection;
		psdstkmatch_efficiency_preselection_accepted 					= eff_preselect.psdstkmatch_efficiency_preselection_accepted;
		psdcharge_efficiency_preselection 								= eff_preselect.psdcharge_efficiency_preselection;
		psdcharge_efficiency_preselection_accepted 						= eff_preselect.psdcharge_efficiency_preselection_accepted;
		stkcharge_efficiency_preselection 								= eff_preselect.stkcharge_efficiency_preselection;
		stkcharge_efficiency_preselection_accepted 						= eff_preselect.stkcharge_efficiency_preselection_accepted;
	}
}

void ntuple::core_reset()
{
	// Trigger
	mip1_trigger 										= false;
	mip2_trigger 										= false;
	HET_trigger 										= false;
	LET_trigger	 										= false;
	MIP_trigger 										= false;
	general_trigger 									= false;
	unbiased_trigger 									= false;
	// STK
	STK_bestTrack_npoints 								= -999;
	STK_bestTrack_nholesX 								= -999;
	STK_bestTrack_nholesY 								= -999;
	STK_bestTrack_slopeX 								= -999;
	STK_bestTrack_slopeY								= -999;
	STK_bestTrack_interceptX 							= -999;
	STK_bestTrack_interceptY 							= -999;
	STK_bestTrack_costheta 								= -999;
	STK_bestTrack_phi 									= -999;
	STK_bestTrack_extr_BGO_topX 						= -999;
	STK_bestTrack_extr_BGO_topY 						= -999;
	STK_bestTrack_STK_BGO_topX_distance 				= -999;
	STK_bestTrack_STK_BGO_topY_distance 				= -999;
	STK_bestTrack_angular_distance_STK_BGO 				= -999;
	STK_chargeX 										= -999;
	STK_chargeY 										= -999;
	STK_charge 											= -999;
	STK_plane_clusters 									= std::vector<int> (DAMPE_stk_planes, -999);
	// BGO
	energy 												= -999;
	energy_corr	 										= -999;
	eLayer 												= std::vector<double>(DAMPE_bgo_nLayers, -999);
	layerBarEnergy 										= std::vector<std::vector<double>> (DAMPE_bgo_nLayers, std::vector<double> (DAMPE_bgo_bars_layer, -999));
	BGOrec_slopeX 										= -999;
	BGOrec_slopeY										= -999;
	BGOrec_interceptX 									= -999;
	BGOrec_interceptY 									= -999;
	trajectoryDirection2D								.SetXYZ(-999, -999, -999);
	sumRms = -999;
	rmsLayer 											= std::vector<double>(DAMPE_bgo_nLayers, -999);
	fracLayer 											= std::vector<double>(DAMPE_bgo_nLayers, -999);
	fracLast 											= -999;
	fracLast_13 										= -999;
	lastBGOLayer 										= -999;
	nBGOentries 										= -999;
	// PSD
	PSD_chargeX 										= -999;
	PSD_chargeY 										= -999;
	PSD_charge 											= -999;
	SPD_STK_match_X_distance							= -999;
	SPD_STK_match_Y_distance							= -999;
	SPD_STK_match_X_distance_fiducial_volume			= -999;
	SPD_STK_match_Y_distance_fiducial_volume			= -999;
	// Classifiers
	xtr 												= -999;
	xtrl 												= -999;
	// NUD
	nud_adc 											= std::vector<int>(DAMPE_NUD_channels, -999);
	nud_total_adc 										= -999;
	nud_max_adc 										= -999;
	nud_max_channel_id 									= -999;
	// Filters
	evtfilter_out_energy_range 							= false;
	evtfilter_evt_triggered 							= false;
	evtfilter_correct_bgo_reco 							= false;
	evtfilter_good_event 								= false;
	evtfilter_geometric 								= false;
	evtfilter_BGO_fiducial 								= false;
	evtfilter_BGO_fiducial_HET 							= false;
	evtfilter_BGO_fiducial_maxElayer_cut 				= false;
	evtfilter_BGO_fiducial_maxBarLayer_cut 				= false;
	evtfilter_BGO_fiducial_BGOTrackContainment_cut 		= false;
	evtfilter_nBarLayer13_cut 							= false;
	evtfilter_maxRms_cut 								= false;
	evtfilter_stk_fiducial_volume 						= false;
	evtfilter_stk_fiducial_volume_X						= false;
	evtfilter_stk_fiducial_volume_Y						= false;
	evtfilter_track_selection_cut 						= false;
	evtfilter_track_selection_cut_no_3hit_recover 		= false;
	evtfilter_psd_fiducial_volume 						= false;
	evtfilter_psd_fiducial_volume_X						= false;
	evtfilter_psd_fiducial_volume_Y						= false;
	evtfilter_psd_stk_match_cut 						= false;
	evtfilter_psd_stk_match_cut_X_view 					= false;
	evtfilter_psd_stk_match_cut_Y_view 					= false;
	evtfilter_psd_charge_cut 							= false;
	evtfilter_psd_charge_cut_no_single_view_recover 	= false;
	evtfilter_stk_charge_cut					 		= false;
	evtfilter_psd_charge_measurement 					= false;
	evtfilter_stk_charge_measurement 					= false;
	evtfilter_xtrl_tight_cut 							= false;
	evtfilter_xtrl_loose_cut 							= false;
	evtfilter_all_cut 									= false;
	// Preselection Cuts
	cut_nBarLayer13 									= false;
	cut_maxRms 											= false;
	cut_track_selection 								= false;
	cut_psd_stk_match 									= false;
	cut_psd_charge 										= false;
	cut_stk_charge 										= false;
	nActiveCuts 										= -999;
	// Preselection Filters
	preselection_maxelayer_cut 							= false;
    preselection_maxbarlayer_cut 						= false;
    preselection_bgotrack_cut 							= false;
    preselection_bgofiducial_cut 						= false;
    preselection_nbarlayer13_cut 						= false;
    preselection_maxrms_cut 							= false;
    preselection_trackselection_cut 					= false;
    preselection_psdstkmatch_cut 						= false;
    preselection_psdcharge_cut 							= false;
    preselection_stkcharge_cut 							= false;
	preselection_maxelayer_lastcut						= false;
	preselection_maxbarlayer_lastcut					= false;
	preselection_bgotrack_lastcut						= false;
	preselection_bgofiducial_lastcut					= false;
	preselection_nbarlayer13_lastcut					= false;
	preselection_maxrms_lastcut							= false;
	preselection_trackselection_lastcut					= false;
	preselection_psdstkmatch_lastcut					= false;
	// Efficiency Preselection Filters
	trigger_efficiency_preselection                     		= false;
	trigger_efficiency_preselection_is_het              		= false;
	trigger_efficiency_preselection_is_let              		= false;
	trigger_efficiency_preselection_is_unb              		= false;
	maxrms_efficiency_preselection                      		= false;
	maxrms_efficiency_preselection_accepted             		= false;
	nbarlayer13_efficiency_preselection                 		= false;
	nbarlayer13_efficiency_preselection_accepted        		= false;
	maxrms_and_nbarlayer13_efficiency_preselection	  			= false;
	maxrms_and_nbarlayer13_efficiency_preselection_accepted		= false;
	track_efficiency_preselection                       		= false;
	track_efficiency_preselection_accepted              		= false;
	psdstkmatch_efficiency_preselection                 		= false;
	psdstkmatch_efficiency_preselection_accepted        		= false;
	psdcharge_efficiency_preselection                   		= false;
	psdcharge_efficiency_preselection_accepted          		= false;
	stkcharge_efficiency_preselection                   		= false;
	stkcharge_efficiency_preselection_accepted          		= false;
}