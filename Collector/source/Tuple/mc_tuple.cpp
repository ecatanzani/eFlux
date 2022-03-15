#include "Tuple/mc_tuple.h"

void mc_tuple::init(const active_cuts &acuts)
{
	// Init Tree
	DmpNtupTree = std::make_unique<TTree>("DmpMCEvtNtup", "DAMPE MC Event nTuple Tree");
	// Set active cuts
	set_active_cuts(acuts);
	// Branch tree
	branch_tree();
}

void mc_tuple::branch_tree()
{
	// Trigger
	DmpNtupTree->Branch("mip1_trigger", &mip1_trigger, "mip1_trigger/O");
	DmpNtupTree->Branch("mip2_trigger", &mip2_trigger, "mip2_trigger/O");
	DmpNtupTree->Branch("HET_trigger", &HET_trigger, "HET_trigger/O");
	DmpNtupTree->Branch("LET_trigger", &LET_trigger, "LET_trigger/O");
	DmpNtupTree->Branch("MIP_trigger", &MIP_trigger, "MIP_trigger/O");
	DmpNtupTree->Branch("general_trigger", &general_trigger, "general_trigger/O");
	DmpNtupTree->Branch("unbiased_trigger", &unbiased_trigger, "unbiased_trigger/O");
	// STK
	DmpNtupTree->Branch("STK_bestTrack_npoints", &STK_bestTrack_npoints, "STK_bestTrack_npoints/I");
	DmpNtupTree->Branch("STK_bestTrack_nholesX", &STK_bestTrack_nholesX, "STK_bestTrack_nholesX/I");
	DmpNtupTree->Branch("STK_bestTrack_nholesY", &STK_bestTrack_nholesY, "STK_bestTrack_nholesY/I");
	DmpNtupTree->Branch("STK_bestTrack_slopeX", &STK_bestTrack_slopeX, "STK_bestTrack_slopeX/D");
	DmpNtupTree->Branch("STK_bestTrack_slopeY", &STK_bestTrack_slopeY, "STK_bestTrack_slopeY/D");
	DmpNtupTree->Branch("STK_bestTrack_interceptX", &STK_bestTrack_interceptX, "STK_bestTrack_interceptX/D");
	DmpNtupTree->Branch("STK_bestTrack_interceptY", &STK_bestTrack_interceptY, "STK_bestTrack_interceptY/D");
	DmpNtupTree->Branch("STK_bestTrack_costheta", &STK_bestTrack_costheta, "STK_bestTrack_costheta/D");
	DmpNtupTree->Branch("STK_bestTrack_phi", &STK_bestTrack_phi, "STK_bestTrack_phi/D");
	DmpNtupTree->Branch("STK_bestTrack_extr_BGO_topX", &STK_bestTrack_extr_BGO_topX, "STK_bestTrack_extr_BGO_topX/D");
	DmpNtupTree->Branch("STK_bestTrack_extr_BGO_topY", &STK_bestTrack_extr_BGO_topY, "STK_bestTrack_extr_BGO_topY/D");
	DmpNtupTree->Branch("STK_bestTrack_STK_BGO_topX_distance", &STK_bestTrack_STK_BGO_topX_distance, "STK_bestTrack_STK_BGO_topX_distance/D");
	DmpNtupTree->Branch("STK_bestTrack_STK_BGO_topY_distance", &STK_bestTrack_STK_BGO_topY_distance, "STK_bestTrack_STK_BGO_topY_distance/D");
	DmpNtupTree->Branch("STK_bestTrack_angular_distance_STK_BGO", &STK_bestTrack_angular_distance_STK_BGO, "STK_bestTrack_angular_distance_STK_BGO/D");
	DmpNtupTree->Branch("STK_chargeX", &STK_chargeX, "STK_chargeX/D");
	DmpNtupTree->Branch("STK_chargeY", &STK_chargeY, "STK_chargeY/D");
	DmpNtupTree->Branch("STK_charge", &STK_charge, "STK_charge/D");
	DmpNtupTree->Branch("STK_plane_clusters", &STK_plane_clusters);	
	// BGO
	DmpNtupTree->Branch("energy", &energy, "energy/D");
	DmpNtupTree->Branch("energy_corr", &energy_corr, "energy_corr/D");
	DmpNtupTree->Branch("energy_corr_w", &energy_corr_w, "energy_corr_w/D");	
	DmpNtupTree->Branch("eLayer", &eLayer);
	DmpNtupTree->Branch("layerBarEnergy", &layerBarEnergy);
	DmpNtupTree->Branch("BGOrec_slopeX", &BGOrec_slopeX, "BGOrec_slopeX/D");
	DmpNtupTree->Branch("BGOrec_slopeY", &BGOrec_slopeY, "BGOrec_slopeY/D");
	DmpNtupTree->Branch("BGOrec_interceptX", &BGOrec_interceptX, "BGOrec_interceptX/D");
	DmpNtupTree->Branch("BGOrec_interceptY", &BGOrec_interceptY, "BGOrec_interceptY/D");
	DmpNtupTree->Branch("BGOrec_trajectoryDirection2D", "TVector3", &trajectoryDirection2D);
	DmpNtupTree->Branch("sumRms", &sumRms, "sumRms/D");
	DmpNtupTree->Branch("rmsLayer", &rmsLayer);
	DmpNtupTree->Branch("fracLayer", &fracLayer);
	DmpNtupTree->Branch("fracLast", &fracLast, "fracLast/D");
	DmpNtupTree->Branch("fracLast_13", &fracLast_13, "fracLast_13/D");
	DmpNtupTree->Branch("lastBGOLayer", &lastBGOLayer, "lastBGOLayer/I");
	DmpNtupTree->Branch("nBGOentries", &nBGOentries, "nBGOentries/I");
	// Simu
	DmpNtupTree->Branch("simu_energy", &simu_energy, "simu_energy/D");
	DmpNtupTree->Branch("simu_energy_w", &simu_energy_w, "simu_energy_w/D");
	DmpNtupTree->Branch("simu_position", "TVector3", &simu_position);
	DmpNtupTree->Branch("simu_momentum",  "TVector3",  &simu_momentum);
	DmpNtupTree->Branch("simu_slope_x", &simu_slope_x, "simu_slope_x/D");
	DmpNtupTree->Branch("simu_slope_y", &simu_slope_y, "simu_slope_y/D");
	DmpNtupTree->Branch("simu_intercept_x", &simu_intercept_x, "simu_intercept_x/D");
	DmpNtupTree->Branch("simu_intercept_y", &simu_intercept_y, "simu_intercept_y/D");
	DmpNtupTree->Branch("simu_radius", &simu_radius, "simu_radius/D");
	DmpNtupTree->Branch("simu_theta", &simu_theta, "simu_theta/D");
	DmpNtupTree->Branch("simu_phi", &simu_phi, "simu_phi/D");
	DmpNtupTree->Branch("simu_flux_w", &simu_flux_w, "simu_flux_w/D");
	DmpNtupTree->Branch("simu_n_particle", &simu_n_particle, "simu_n_particle/I");
	DmpNtupTree->Branch("simu_cos_x", &simu_cos_x, "simu_cos_x/D");
	DmpNtupTree->Branch("simu_cos_y", &simu_cos_y, "simu_cos_y/D");
	DmpNtupTree->Branch("simu_cos_z", &simu_cos_z, "simu_cos_z/D");
	DmpNtupTree->Branch("simu_charge", &simu_charge, "simu_charge/D");
	DmpNtupTree->Branch("simu_zenith", &simu_zenith, "simu_zenith/D");
	DmpNtupTree->Branch("simu_azimuth", &simu_azimuth, "simu_azimuth/D");
	DmpNtupTree->Branch("simu_w", &simu_w, "simu_w/D");
	DmpNtupTree->Branch("simu_PDG", &simu_PDG, "simu_PDG/D");
	DmpNtupTree->Branch("simu_geocut", &simu_geocut, "simu_geocut/D");
	DmpNtupTree->Branch("simu_thruthtrajectory_x", &simu_thruthtrajectory_x);
	DmpNtupTree->Branch("simu_thruthtrajectory_y", &simu_thruthtrajectory_y);
	DmpNtupTree->Branch("simu_thruthtrajectory_z", &simu_thruthtrajectory_z);
	DmpNtupTree->Branch("simu_thruthtrajectory_energy", &simu_thruthtrajectory_energy);
	DmpNtupTree->Branch("simu_thruthtrajectory_start_x", &simu_thruthtrajectory_start_x);
	DmpNtupTree->Branch("simu_thruthtrajectory_start_y", &simu_thruthtrajectory_start_y);
	DmpNtupTree->Branch("simu_thruthtrajectory_start_z", &simu_thruthtrajectory_start_z);
	DmpNtupTree->Branch("simu_thruthtrajectory_stop_x", &simu_thruthtrajectory_stop_x);
	DmpNtupTree->Branch("simu_thruthtrajectory_stop_y", &simu_thruthtrajectory_stop_y);
	DmpNtupTree->Branch("simu_thruthtrajectory_stop_z", &simu_thruthtrajectory_stop_z);
	DmpNtupTree->Branch("simu_thruthtrajectory_trackID", &simu_thruthtrajectory_trackID);
	DmpNtupTree->Branch("simu_thruthtrajectory_parentID", &simu_thruthtrajectory_parentID);
	DmpNtupTree->Branch("simu_thruthtrajectory_charge", &simu_thruthtrajectory_charge);
	DmpNtupTree->Branch("simu_thruthtrajectory_PDG", &simu_thruthtrajectory_PDG);
	DmpNtupTree->Branch("simu_thruthtrajectory_stop_index", &simu_thruthtrajectory_stop_index);
	// PSD
	DmpNtupTree->Branch("PSD_chargeX", &PSD_chargeX, "PSD_chargeX/D");
	DmpNtupTree->Branch("PSD_chargeY", &PSD_chargeY, "PSD_chargeY/D");
	DmpNtupTree->Branch("PSD_charge", &PSD_charge, "PSD_charge/D");
	DmpNtupTree->Branch("SPD_STK_match_X_distance", &SPD_STK_match_X_distance, "SPD_STK_match_X_distance/D");
	DmpNtupTree->Branch("SPD_STK_match_Y_distance", &SPD_STK_match_Y_distance, "SPD_STK_match_Y_distance/D");
	DmpNtupTree->Branch("SPD_STK_match_X_distance_fiducial_volume", &SPD_STK_match_X_distance_fiducial_volume, "SPD_STK_match_X_distance_fiducial_volume/D");
	DmpNtupTree->Branch("SPD_STK_match_Y_distance_fiducial_volume", &SPD_STK_match_Y_distance_fiducial_volume, "SPD_STK_match_Y_distance_fiducial_volume/D");
	// NUD
	DmpNtupTree->Branch("nud_adc", &nud_adc);
	DmpNtupTree->Branch("nud_total_adc", &nud_total_adc, "nud_total_adc/I");
	DmpNtupTree->Branch("nud_max_adc", &nud_max_adc, "nud_max_adc/I");
	DmpNtupTree->Branch("nud_max_channel_id", &nud_max_channel_id, "nud_max_channel_id/I");
	// Classifiers
	DmpNtupTree->Branch("xtr", &xtr, "xtr/D");
	DmpNtupTree->Branch("xtrl", &xtrl, "xtrl/D");
	// Cuts
	DmpNtupTree->Branch("cut_nBarLayer13", &cut_nBarLayer13, "cut_nBarLayer13/O");
	DmpNtupTree->Branch("cut_maxRms", &cut_maxRms, "cut_maxRms/O");
	DmpNtupTree->Branch("cut_track_selection", &cut_track_selection, "cut_track_selection/O");
	DmpNtupTree->Branch("cut_psd_stk_match", &cut_psd_stk_match, "cut_psd_stk_match/O");
	DmpNtupTree->Branch("cut_psd_charge", &cut_psd_charge, "cut_psd_charge/O");
	DmpNtupTree->Branch("cut_stk_charge", &cut_stk_charge, "cut_stk_charge/O");
	DmpNtupTree->Branch("nActiveCuts", &nActiveCuts, "nActiveCuts/I");
	// Filters
	DmpNtupTree->Branch("evtfilter_out_energy_range", &evtfilter_out_energy_range, "evtfilter_out_energy_range/O");
	DmpNtupTree->Branch("evtfilter_geometric_before_trigger", &evtfilter_geometric_before_trigger, "evtfilter_geometric_before_trigger/O");
	DmpNtupTree->Branch("evtfilter_trigger_check", &evtfilter_trigger_check, "evtfilter_trigger_check/O");
	DmpNtupTree->Branch("evtfilter_evt_triggered", &evtfilter_evt_triggered, "evtfilter_evt_triggered/O");
	DmpNtupTree->Branch("evtfilter_correct_bgo_reco", &evtfilter_correct_bgo_reco, "evtfilter_correct_bgo_reco/O");
	DmpNtupTree->Branch("evtfilter_good_event", &evtfilter_good_event, "evtfilter_good_event/O");
	DmpNtupTree->Branch("evtfilter_geometric", &evtfilter_geometric, "evtfilter_geometric/O");
	DmpNtupTree->Branch("evtfilter_BGO_fiducial", &evtfilter_BGO_fiducial, "evtfilter_BGO_fiducial/O");
	DmpNtupTree->Branch("evtfilter_BGO_fiducial_HET", &evtfilter_BGO_fiducial_HET, "evtfilter_BGO_fiducial_HET/O");
	DmpNtupTree->Branch("evtfilter_BGO_fiducial_maxElayer_cut", &evtfilter_BGO_fiducial_maxElayer_cut, "evtfilter_BGO_fiducial_maxElayer_cut/O");
	DmpNtupTree->Branch("evtfilter_BGO_fiducial_maxBarLayer_cut", &evtfilter_BGO_fiducial_maxBarLayer_cut, "evtfilter_BGO_fiducial_maxBarLayer_cut/O");
	DmpNtupTree->Branch("evtfilter_BGO_fiducial_BGOTrackContainment_cut", &evtfilter_BGO_fiducial_BGOTrackContainment_cut, "evtfilter_BGO_fiducial_BGOTrackContainment_cut/O");
	DmpNtupTree->Branch("evtfilter_nBarLayer13_cut", &evtfilter_nBarLayer13_cut, "evtfilter_nBarLayer13_cut/O");
	DmpNtupTree->Branch("evtfilter_maxRms_cut", &evtfilter_maxRms_cut, "evtfilter_maxRms_cut/O");
	DmpNtupTree->Branch("evtfilter_track_selection_cut", &evtfilter_track_selection_cut, "evtfilter_track_selection_cut/O");
	DmpNtupTree->Branch("evtfilter_track_selection_cut_no_3hit_recover", &evtfilter_track_selection_cut_no_3hit_recover, "evtfilter_track_selection_cut_no_3hit_recover/O");
	DmpNtupTree->Branch("evtfilter_psd_stk_match_cut", &evtfilter_psd_stk_match_cut, "evtfilter_psd_stk_match_cut/O");
	DmpNtupTree->Branch("evtfilter_psd_stk_match_cut_X_view", &evtfilter_psd_stk_match_cut_X_view, "evtfilter_psd_stk_match_cut_X_view/O");
	DmpNtupTree->Branch("evtfilter_psd_stk_match_cut_Y_view", &evtfilter_psd_stk_match_cut_Y_view, "evtfilter_psd_stk_match_cut_Y_view/O");
	DmpNtupTree->Branch("evtfilter_psd_stk_match_cut_psd_fiducial_volume", &evtfilter_psd_stk_match_cut_psd_fiducial_volume, "evtfilter_psd_stk_match_cut_psd_fiducial_volume/O");
	DmpNtupTree->Branch("evtfilter_psd_stk_match_cut_psd_fiducial_volume_X", &evtfilter_psd_stk_match_cut_psd_fiducial_volume_X, "evtfilter_psd_stk_match_cut_psd_fiducial_volume_X/O");
	DmpNtupTree->Branch("evtfilter_psd_stk_match_cut_psd_fiducial_volume_Y", &evtfilter_psd_stk_match_cut_psd_fiducial_volume_Y, "evtfilter_psd_stk_match_cut_psd_fiducial_volume_Y/O");
	DmpNtupTree->Branch("evtfilter_psd_charge_cut", &evtfilter_psd_charge_cut, "evtfilter_psd_charge_cut/O");
	DmpNtupTree->Branch("evtfilter_psd_charge_cut_no_single_view_recover", &evtfilter_psd_charge_cut_no_single_view_recover, "evtfilter_psd_charge_cut_no_single_view_recover/O");
	DmpNtupTree->Branch("evtfilter_stk_charge_cut", &evtfilter_stk_charge_cut, "evtfilter_stk_charge_cut/O");
	DmpNtupTree->Branch("evtfilter_psd_charge_measurement", &evtfilter_psd_charge_measurement, "evtfilter_psd_charge_measurement/O");
	DmpNtupTree->Branch("evtfilter_stk_charge_measurement", &evtfilter_stk_charge_measurement, "evtfilter_stk_charge_measurement/O");
	DmpNtupTree->Branch("evtfilter_xtrl_tight_cut", &evtfilter_xtrl_tight_cut, "evtfilter_xtrl_tight_cut/O");
	DmpNtupTree->Branch("evtfilter_xtrl_loose_cut", &evtfilter_xtrl_loose_cut, "evtfilter_xtrl_loose_cut/O");
	DmpNtupTree->Branch("evtfilter_all_cut", &evtfilter_all_cut, "evtfilter_all_cut/O");
	// Preselection Filters
	DmpNtupTree->Branch("preselection_maxelayer_cut", &preselection_maxelayer_cut, "preselection_maxelayer_cut/O");
	DmpNtupTree->Branch("preselection_maxbarlayer_cut", &preselection_maxbarlayer_cut, "preselection_maxbarlayer_cut/O");
	DmpNtupTree->Branch("preselection_bgotrack_cut", &preselection_bgotrack_cut, "preselection_bgotrack_cut/O");
	DmpNtupTree->Branch("preselection_bgofiducial_cut", &preselection_bgofiducial_cut, "preselection_bgofiducial_cut/O");
	DmpNtupTree->Branch("preselection_nbarlayer13_cut", &preselection_nbarlayer13_cut, "preselection_nbarlayer13_cut/O");
	DmpNtupTree->Branch("preselection_maxrms_cut", &preselection_maxrms_cut, "preselection_maxrms_cut/O");
	DmpNtupTree->Branch("preselection_trackselection_cut", &preselection_trackselection_cut, "preselection_trackselection_cut/O");
	DmpNtupTree->Branch("preselection_psdstkmatch_cut", &preselection_psdstkmatch_cut, "preselection_psdstkmatch_cut/O");
	DmpNtupTree->Branch("preselection_psdcharge_cut", &preselection_psdcharge_cut, "preselection_psdcharge_cut/O");
	DmpNtupTree->Branch("preselection_stkcharge_cut", &preselection_stkcharge_cut, "preselection_stkcharge_cut/O");

	DmpNtupTree->Branch("preselection_maxelayer_lastcut", &preselection_maxelayer_lastcut, "preselection_maxelayer_lastcut/O");
	DmpNtupTree->Branch("preselection_maxbarlayer_lastcut", &preselection_maxbarlayer_lastcut, "preselection_maxbarlayer_lastcut/O");
	DmpNtupTree->Branch("preselection_bgotrack_lastcut", &preselection_bgotrack_lastcut, "preselection_bgotrack_lastcut/O");
	DmpNtupTree->Branch("preselection_bgofiducial_lastcut", &preselection_bgofiducial_lastcut, "preselection_bgofiducial_lastcut/O");
	DmpNtupTree->Branch("preselection_nbarlayer13_lastcut", &preselection_nbarlayer13_lastcut, "preselection_nbarlayer13_lastcut/O");
	DmpNtupTree->Branch("preselection_maxrms_lastcut", &preselection_maxrms_lastcut, "preselection_maxrms_lastcut/O");
	DmpNtupTree->Branch("preselection_trackselection_lastcut", &preselection_trackselection_lastcut, "preselection_trackselection_lastcut/O");
	DmpNtupTree->Branch("preselection_psdstkmatch_lastcut", &preselection_psdstkmatch_lastcut, "preselection_psdstkmatch_lastcut/O");
	// Efficiency Preselection Filters
	DmpNtupTree->Branch("trigger_efficiency_preselection", &trigger_efficiency_preselection, "trigger_efficiency_preselection/O");
	DmpNtupTree->Branch("trigger_efficiency_preselection_is_het", &trigger_efficiency_preselection_is_het, "trigger_efficiency_preselection_is_het/O");
	DmpNtupTree->Branch("trigger_efficiency_preselection_is_let", &trigger_efficiency_preselection_is_let, "trigger_efficiency_preselection_is_let/O");
	DmpNtupTree->Branch("trigger_efficiency_preselection_is_unb", &trigger_efficiency_preselection_is_unb, "trigger_efficiency_preselection_is_unb/O");
	DmpNtupTree->Branch("maxrms_efficiency_preselection", &maxrms_efficiency_preselection, "maxrms_efficiency_preselection/O");
	DmpNtupTree->Branch("maxrms_efficiency_preselection_accepted", &maxrms_efficiency_preselection_accepted, "maxrms_efficiency_preselection_accepted/O");
	DmpNtupTree->Branch("nbarlayer13_efficiency_preselection", &nbarlayer13_efficiency_preselection, "nbarlayer13_efficiency_preselection/O");
	DmpNtupTree->Branch("nbarlayer13_efficiency_preselection_accepted", &nbarlayer13_efficiency_preselection_accepted, "nbarlayer13_efficiency_preselection_accepted/O");
	DmpNtupTree->Branch("maxrms_and_nbarlayer13_efficiency_preselection", &maxrms_and_nbarlayer13_efficiency_preselection, "maxrms_and_nbarlayer13_efficiency_preselection/O");
	DmpNtupTree->Branch("maxrms_and_nbarlayer13_efficiency_preselection_accepted", &maxrms_and_nbarlayer13_efficiency_preselection_accepted, "maxrms_and_nbarlayer13_efficiency_preselection_accepted/O");
	DmpNtupTree->Branch("track_efficiency_preselection", &track_efficiency_preselection, "track_efficiency_preselection/O");
	DmpNtupTree->Branch("track_efficiency_preselection_accepted", &track_efficiency_preselection_accepted, "track_efficiency_preselection_accepted/O");
	DmpNtupTree->Branch("psdstkmatch_efficiency_preselection", &psdstkmatch_efficiency_preselection, "psdstkmatch_efficiency_preselection/O");
	DmpNtupTree->Branch("psdstkmatch_efficiency_preselection_accepted", &psdstkmatch_efficiency_preselection_accepted, "psdstkmatch_efficiency_preselection_accepted/O");
	DmpNtupTree->Branch("psdcharge_efficiency_preselection", &psdcharge_efficiency_preselection, "psdcharge_efficiency_preselection/O");
	DmpNtupTree->Branch("psdcharge_efficiency_preselection_accepted", &psdcharge_efficiency_preselection_accepted, "psdcharge_efficiency_preselection_accepted/O");
}

void mc_tuple::Fill(
	const std::shared_ptr<_tmp_filter> _filter_res,
	const std::vector<int> _stk_clusters_on_plane,
	const std::shared_ptr<_tmp_psd> _psd_res,
	const std::shared_ptr<_tmp_bgo> _bgo_res,
	const std::shared_ptr<_tmp_simu> _simu_res,
	const std::shared_ptr<_tmp_energy> _energy_res,
	const std::shared_ptr<_tmp_nud> _nud_res)
{
	fill_trigger_info(_filter_res->evt_trigger_info);
	fill_filter_info(_filter_res->output);
	fill_preselectionfilters_info(_filter_res->preselection_output);
	fill_preselectionfiltersefficiency_info(_filter_res->efficiency_output);
	fill_psd_info(
		_psd_res->SPD_STK_match_X_distance,
		_psd_res->SPD_STK_match_Y_distance,
		_psd_res->SPD_STK_match_X_distance_fiducial_volume,
		_psd_res->SPD_STK_match_Y_distance_fiducial_volume);
	fill_stk_info(
		_filter_res->evt_best_track, 
		_stk_clusters_on_plane);
	fill_bgo_info(
		_energy_res->raw,
		_energy_res->correct,
		_bgo_res->layer_energies,
		_bgo_res->layer_bar_energies,
		_bgo_res->slope,
		_bgo_res->intercept,
		_bgo_res->trajectory2D,
		_bgo_res->sumrms,
		_bgo_res->sumrms_layer,
		_bgo_res->energy_fraction_layer,
		_bgo_res->energy_fraction_last_layer,
		_bgo_res->energy_fraction_13th_layer,
		_bgo_res->last_energy_layer,
		_bgo_res->hits);
	fill_simu_info(_simu_res, _energy_res);
	fill_psdcharge_info(_filter_res->evt_psd_charge);
	fill_stkcharge_info(_filter_res->evt_stk_charge);
	fill_classifier_info(_filter_res->evt_bgo_classifier);
	fill_nud_info(
		_nud_res->adc,
		_nud_res->total_adc,
		_nud_res->max_adc,
		_nud_res->max_channel_ID);
	DmpNtupTree->Fill();
}

void mc_tuple::fill_filter_info(const filter_output &output)
{
	evtfilter_out_energy_range 								= output.out_energy_range;
	evtfilter_geometric_before_trigger 						= output.geometric_before_trigger;
	evtfilter_trigger_check 								= output.trigger_check;
	evtfilter_evt_triggered 								= output.evt_triggered;
	evtfilter_correct_bgo_reco 								= output.correct_bgo_reco;
	evtfilter_good_event 									= output.good_event;
	evtfilter_geometric 									= output.geometric;
	evtfilter_BGO_fiducial 									= output.BGO_fiducial;
	evtfilter_BGO_fiducial_HET 								= output.BGO_fiducial_HET;
	evtfilter_BGO_fiducial_maxElayer_cut 					= output.BGO_fiducial_maxElayer_cut;
	evtfilter_BGO_fiducial_maxBarLayer_cut 					= output.BGO_fiducial_maxBarLayer_cut;
	evtfilter_BGO_fiducial_BGOTrackContainment_cut 			= output.BGO_fiducial_BGOTrackContainment_cut;
	evtfilter_nBarLayer13_cut 								= output.nBarLayer13_cut;
	evtfilter_maxRms_cut 									= output.maxRms_cut;
	evtfilter_track_selection_cut 							= output.track_selection_cut;
	evtfilter_track_selection_cut_no_3hit_recover 			= output.track_selection_cut_no_3hit_recover;
	evtfilter_psd_stk_match_cut 							= output.psd_stk_match_cut;
	evtfilter_psd_stk_match_cut_X_view 						= output.psd_stk_match_cut_x;
	evtfilter_psd_stk_match_cut_Y_view 						= output.psd_stk_match_cut_y;
	evtfilter_psd_stk_match_cut_psd_fiducial_volume	 		= output.psd_stk_match_cut_psd_fiducial_volume;
	evtfilter_psd_stk_match_cut_psd_fiducial_volume_X 		= output.psd_stk_match_cut_psd_fiducial_volume_X;
	evtfilter_psd_stk_match_cut_psd_fiducial_volume_Y 		= output.psd_stk_match_cut_psd_fiducial_volume_Y;
	evtfilter_psd_charge_cut 								= output.psd_charge_cut;
	evtfilter_psd_charge_cut_no_single_view_recover 		= output.psd_charge_cut_no_single_view_recover;
	evtfilter_stk_charge_cut 								= output.stk_charge_cut;
	evtfilter_psd_charge_measurement 						= output.psd_charge_measurement;
	evtfilter_stk_charge_measurement 						= output.stk_charge_measurement;
	evtfilter_xtrl_tight_cut 								= output.xtrl_tight_cut;
	evtfilter_xtrl_loose_cut 								= output.xtrl_loose_cut;
	evtfilter_all_cut 										= output.all_cut;
}

void mc_tuple::fill_simu_info(
	const std::shared_ptr<_tmp_simu> _simu_res,
	const std::shared_ptr<_tmp_energy> _energy_res)
{
	energy_corr_w 						= _energy_res->correct_w;
	simu_energy_w 						= _energy_res->simu_w;
	simu_position 						= _simu_res->position;
	simu_momentum 						= _simu_res->momentum;
	simu_energy 						= _energy_res->simu;
	simu_slope_x 						= _simu_res->momentum.Z() ? _simu_res->momentum.X() / _simu_res->momentum.Z() : -999;
	simu_slope_y 						= _simu_res->momentum.Z() ? _simu_res->momentum.Y() / _simu_res->momentum.Z() : -999;
	simu_intercept_x 					= _simu_res->position.X() - simu_slope_x * _simu_res->position.Z();
	simu_intercept_y 					= _simu_res->position.Y() - simu_slope_y * _simu_res->position.Z();
	simu_radius 						= _simu_res->radius;
	simu_theta 							= _simu_res->theta;
	simu_phi 							= _simu_res->phi;
	simu_flux_w 						= _simu_res->flux_w;
	simu_n_particle 					= _simu_res->n_particle;
	simu_cos_x 							= _simu_res->cos_x;
	simu_cos_y 							= _simu_res->cos_y;
	simu_cos_z 							= _simu_res->cos_z;
	simu_charge 						= _simu_res->charge;
	simu_zenith 						= _simu_res->zenith;
	simu_azimuth 						= _simu_res->azimuth;
	simu_w 								= _simu_res->w;
	simu_PDG 							= _simu_res->PDG;
	simu_geocut 						= _simu_res->geocut;
	simu_thruthtrajectory_x 			= _simu_res->thruthtrajectory_x;
	simu_thruthtrajectory_y 			= _simu_res->thruthtrajectory_y;
	simu_thruthtrajectory_z 			= _simu_res->thruthtrajectory_z;
	simu_thruthtrajectory_energy 		= _simu_res->thruthtrajectory_energy;
	simu_thruthtrajectory_start_x 		= _simu_res->thruthtrajectory_start_x;
	simu_thruthtrajectory_start_y 		= _simu_res->thruthtrajectory_start_y;
	simu_thruthtrajectory_start_z 		= _simu_res->thruthtrajectory_start_z;
	simu_thruthtrajectory_stop_x 		= _simu_res->thruthtrajectory_stop_x;
	simu_thruthtrajectory_stop_y 		= _simu_res->thruthtrajectory_stop_y;
	simu_thruthtrajectory_stop_z 		= _simu_res->thruthtrajectory_stop_z;
	simu_thruthtrajectory_trackID 		= _simu_res->thruthtrajectory_trackID;
	simu_thruthtrajectory_parentID 		= _simu_res->thruthtrajectory_parentID;
	simu_thruthtrajectory_charge 		= _simu_res->thruthtrajectory_charge;
	simu_thruthtrajectory_PDG 			= _simu_res->thruthtrajectory_PDG;
	simu_thruthtrajectory_stop_index 	= _simu_res->thruthtrajectory_stop_index;
}

void mc_tuple::Reset()
{
	core_reset();
	reset_simu_info();
}

void mc_tuple::reset_simu_info()
{
	simu_energy 						= -999;
	simu_energy_w 						= -999;
	energy_corr_w 						= -999;
	simu_position						.SetXYZ(-999, -999, -999);
	simu_momentum						.SetXYZ(-999, -999, -999);
	simu_slope_x 						= -999;
	simu_slope_y 						= -999;
	simu_intercept_x 					= -999;
	simu_intercept_y 					= -999;
	simu_radius 						= -999;
	simu_theta 							= -999;
	simu_phi 							= -999;
	simu_flux_w 						= -999;
	simu_n_particle 					= -999;
	simu_cos_x 							= -999;
	simu_cos_y 							= -999;
	simu_cos_z 							= -999;
	simu_charge 						= -999;
	simu_zenith 						= -999;
	simu_azimuth 						= -999;
	simu_w 								= -999;
	simu_PDG 							= -999;
	simu_geocut 						= -999;
	simu_thruthtrajectory_x 			= std::vector<double>();
	simu_thruthtrajectory_y 			= std::vector<double>();
	simu_thruthtrajectory_z 			= std::vector<double>();
	simu_thruthtrajectory_energy 		= std::vector<double>();
	simu_thruthtrajectory_start_x 		= std::vector<double>();
	simu_thruthtrajectory_start_y 		= std::vector<double>();
	simu_thruthtrajectory_start_z		= std::vector<double>();
	simu_thruthtrajectory_stop_x 		= std::vector<double>();
	simu_thruthtrajectory_stop_y 		= std::vector<double>();
	simu_thruthtrajectory_stop_z 		= std::vector<double>();
	simu_thruthtrajectory_trackID 		= std::vector<double>();
	simu_thruthtrajectory_parentID 		= std::vector<double>();
	simu_thruthtrajectory_charge 		= std::vector<double>();
	simu_thruthtrajectory_PDG 			= std::vector<double>();
	simu_thruthtrajectory_stop_index 	= std::vector<double>();
	evtfilter_geometric_before_trigger 	= false;
	evtfilter_trigger_check 			= false;
}