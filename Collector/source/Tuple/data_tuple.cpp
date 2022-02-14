#include "Tuple/data_tuple.h"

void data_tuple::init(const active_cuts &acuts)
{
	// Init Tree
	DmpNtupTree = std::make_unique<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	// Set active cuts
	set_active_cuts(acuts);
	// Branch tree
	branch_tree();
}

void data_tuple::branch_tree()
{
	// Trigger
	DmpNtupTree->Branch(
		"mip1_trigger",
		&mip1_trigger,
		"mip1_trigger/O");
	DmpNtupTree->Branch(
		"mip2_trigger",
		&mip2_trigger,
		"mip2_trigger/O");
	DmpNtupTree->Branch(
		"HET_trigger",
		&HET_trigger,
		"HET_trigger/O");
	DmpNtupTree->Branch(
		"LET_trigger",
		&LET_trigger,
		"LET_trigger/O");
	DmpNtupTree->Branch(
		"MIP_trigger",
		&MIP_trigger,
		"MIP_trigger/O");
	DmpNtupTree->Branch(
		"general_trigger",
		&general_trigger,
		"general_trigger/O");
	DmpNtupTree->Branch(
		"unbiased_trigger",
		&unbiased_trigger,
		"unbiased_trigger/O");
	// Event time
	DmpNtupTree->Branch(
		"second",
		&second,
		"second/i");
	// STK
	DmpNtupTree->Branch(
		"STK_bestTrack_npoints",
		&STK_bestTrack_npoints,
		"STK_bestTrack_npoints/I");
	DmpNtupTree->Branch(
		"STK_bestTrack_nholesX",
		&STK_bestTrack_nholesX,
		"STK_bestTrack_nholesX/I");
	DmpNtupTree->Branch(
		"STK_bestTrack_nholesY",
		&STK_bestTrack_nholesY,
		"STK_bestTrack_nholesY/I");
	DmpNtupTree->Branch(
		"STK_bestTrack_slopeX",
		&STK_bestTrack_slopeX,
		"STK_bestTrack_slopeX/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_slopeY",
		&STK_bestTrack_slopeY,
		"STK_bestTrack_slopeY/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_interceptX",
		&STK_bestTrack_interceptX,
		"STK_bestTrack_interceptX/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_interceptY",
		&STK_bestTrack_interceptY,
		"STK_bestTrack_interceptY/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_costheta",
		&STK_bestTrack_costheta,
		"STK_bestTrack_costheta/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_phi",
		&STK_bestTrack_phi,
		"STK_bestTrack_phi/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_extr_BGO_topX",
		&STK_bestTrack_extr_BGO_topX,
		"STK_bestTrack_extr_BGO_topX/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_extr_BGO_topY",
		&STK_bestTrack_extr_BGO_topY,
		"STK_bestTrack_extr_BGO_topY/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_STK_BGO_topX_distance",
		&STK_bestTrack_STK_BGO_topX_distance,
		"STK_bestTrack_STK_BGO_topX_distance/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_STK_BGO_topY_distance",
		&STK_bestTrack_STK_BGO_topY_distance,
		"STK_bestTrack_STK_BGO_topY_distance/D");
	DmpNtupTree->Branch(
		"STK_bestTrack_angular_distance_STK_BGO",
		&STK_bestTrack_angular_distance_STK_BGO,
		"STK_bestTrack_angular_distance_STK_BGO/D");
	DmpNtupTree->Branch(
		"STK_chargeX",
		&STK_chargeX,
		"STK_chargeX/D");
	DmpNtupTree->Branch(
		"STK_chargeY",
		&STK_chargeY,
		"STK_chargeY/D");
	DmpNtupTree->Branch(
		"STK_charge",
		&STK_charge,
		"STK_charge/D");
	DmpNtupTree->Branch(
		"STK_plane_clusters",
		&STK_plane_clusters);
	// BGO
	DmpNtupTree->Branch(
		"energy",
		&energy,
		"energy/D");
	DmpNtupTree->Branch(
		"energy_corr",
		&energy_corr,
		"energy_corr/D");
	DmpNtupTree->Branch(
		"eLayer",
		&eLayer);
	DmpNtupTree->Branch(
		"layerBarEnergy",
		&layerBarEnergy);
	DmpNtupTree->Branch(
		"BGOrec_slopeX",
		&BGOrec_slopeX,
		"BGOrec_slopeX/D");
	DmpNtupTree->Branch(
		"BGOrec_slopeY",
		&BGOrec_slopeY,
		"BGOrec_slopeY/D");
	DmpNtupTree->Branch(
		"BGOrec_interceptX",
		&BGOrec_interceptX,
		"BGOrec_interceptX/D");
	DmpNtupTree->Branch(
		"BGOrec_interceptY",
		&BGOrec_interceptY,
		"BGOrec_interceptY/D");
	DmpNtupTree->Branch(
		"BGOrec_trajectoryDirection2D", 
		"TVector3", 
		&trajectoryDirection2D);
	DmpNtupTree->Branch(
		"sumRms",
		&sumRms,
		"sumRms/D");
	DmpNtupTree->Branch(
		"rmsLayer",
		&rmsLayer);
	DmpNtupTree->Branch(
		"fracLayer",
		&fracLayer);
	DmpNtupTree->Branch(
		"fracLast",
		&fracLast,
		"fracLast/D");
	DmpNtupTree->Branch(
		"fracLast_13",
		&fracLast_13,
		"fracLast_13/D");
	DmpNtupTree->Branch(
		"lastBGOLayer",
		&lastBGOLayer,
		"lastBGOLayer/I");
	DmpNtupTree->Branch(
		"nBGOentries",
		&nBGOentries,
		"nBGOentries/I");
	// PSD
	DmpNtupTree->Branch(
		"PSD_chargeX",
		&PSD_chargeX,
		"PSD_chargeX/D");
	DmpNtupTree->Branch(
		"PSD_chargeY",
		&PSD_chargeY,
		"PSD_chargeY/D");
	DmpNtupTree->Branch(
		"PSD_charge",
		&PSD_charge,
		"PSD_charge/D");
	// NUD
	DmpNtupTree->Branch(
		"nud_adc",
		&nud_adc);
	DmpNtupTree->Branch(
		"nud_total_adc",
		&nud_total_adc,
		"nud_total_adc/I");
	DmpNtupTree->Branch(
		"nud_max_adc",
		&nud_max_adc,
		"nud_max_adc/I");
	DmpNtupTree->Branch(
		"nud_max_channel_id",
		&nud_max_channel_id,
		"nud_max_channel_id/I");
	// Classifiers
	DmpNtupTree->Branch(
		"xtr",
		&xtr,
		"xtr/D");
	DmpNtupTree->Branch(
		"xtrl",
		&xtrl,
		"xtrl/D");
	// Attitude
	DmpNtupTree->Branch(
		"glat",
		&glat,
		"glat/D");
	DmpNtupTree->Branch(
		"glon",
		&glon,
		"glon/D");
	DmpNtupTree->Branch(
		"geo_lat",
		&geo_lat,
		"geo_lat/D");
	DmpNtupTree->Branch(
		"geo_lon",
		&geo_lon,
		"geo_lon/D");
	DmpNtupTree->Branch(
		"ra_zenith",
		&ra_zenith,
		"ra_zenith/D");
	DmpNtupTree->Branch(
		"dec_zenith",
		&dec_zenith,
		"dec_zenith/D");
	DmpNtupTree->Branch(
		"ra_scz",
		&ra_scz,
		"ra_scz/D");
	DmpNtupTree->Branch(
		"dec_scz",
		&dec_scz,
		"dec_scz/D");
	DmpNtupTree->Branch(
		"ra_scx",
		&ra_scx,
		"ra_scx/D");
	DmpNtupTree->Branch(
		"dec_scx",
		&dec_scx,
		"dec_scx/D");
	DmpNtupTree->Branch(
		"ra_scy",
		&ra_scy,
		"ra_scy/D");
	DmpNtupTree->Branch(
		"dec_scy",
		&dec_scy,
		"dec_scy/D");
	DmpNtupTree->Branch(
		"cutoff",
		&cutoff,
		"cutoff/D");
	// Preselection Cuts
	DmpNtupTree->Branch(
		"cut_nBarLayer13",
		&cut_nBarLayer13,
		"cut_nBarLayer13/O");
	DmpNtupTree->Branch(
		"cut_maxRms",
		&cut_maxRms,
		"cut_maxRms/O");
	DmpNtupTree->Branch(
		"cut_track_selection",
		&cut_track_selection,
		"cut_track_selection/O");
	DmpNtupTree->Branch(
		"cut_psd_stk_match",
		&cut_psd_stk_match,
		"cut_psd_stk_match/O");
	DmpNtupTree->Branch(
		"cut_psd_charge",
		&cut_psd_charge,
		"cut_psd_charge/O");
	DmpNtupTree->Branch(
		"cut_stk_charge",
		&cut_stk_charge,
		"cut_stk_charge/O");
	DmpNtupTree->Branch(
		"nActiveCuts",
		&nActiveCuts,
		"nActiveCuts/I");
	// Filters
	DmpNtupTree->Branch(
		"evtfilter_out_energy_range",
		&evtfilter_out_energy_range,
		"evtfilter_out_energy_range/O");
	DmpNtupTree->Branch(
		"evtfilter_evt_in_saa",
		&evtfilter_evt_in_saa,
		"evtfilter_evt_in_saa/O");
	DmpNtupTree->Branch(
		"evtfilter_evt_triggered",
		&evtfilter_evt_triggered,
		"evtfilter_evt_triggered/O");
	DmpNtupTree->Branch(
		"evtfilter_correct_bgo_reco",
		&evtfilter_correct_bgo_reco,
		"evtfilter_correct_bgo_reco/O");
	DmpNtupTree->Branch(
		"evtfilter_good_event",
		&evtfilter_good_event,
		"evtfilter_good_event/O");
	DmpNtupTree->Branch(
		"evtfilter_geometric",
		&evtfilter_geometric,
		"evtfilter_geometric/O");
	DmpNtupTree->Branch(
		"evtfilter_BGO_fiducial",
		&evtfilter_BGO_fiducial,
		"evtfilter_BGO_fiducial/O");
	DmpNtupTree->Branch(
		"evtfilter_BGO_fiducial_maxElayer_cut",
		&evtfilter_BGO_fiducial_maxElayer_cut,
		"evtfilter_BGO_fiducial_maxElayer_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_BGO_fiducial_maxBarLayer_cut",
		&evtfilter_BGO_fiducial_maxBarLayer_cut,
		"evtfilter_BGO_fiducial_maxBarLayer_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_BGO_fiducial_BGOTrackContainment_cut",
		&evtfilter_BGO_fiducial_BGOTrackContainment_cut, "evtfilter_BGO_fiducial_BGOTrackContainment_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_nBarLayer13_cut",
		&evtfilter_nBarLayer13_cut,
		"evtfilter_nBarLayer13_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_maxRms_cut",
		&evtfilter_maxRms_cut,
		"evtfilter_maxRms_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_track_selection_cut",
		&evtfilter_track_selection_cut,
		"evtfilter_track_selection_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_psd_stk_match_cut",
		&evtfilter_psd_stk_match_cut,
		"evtfilter_psd_stk_match_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_psd_charge_cut",
		&evtfilter_psd_charge_cut,
		"evtfilter_psd_charge_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_stk_charge_cut",
		&evtfilter_stk_charge_cut,
		"evtfilter_stk_charge_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_psd_charge_measurement",
		&evtfilter_psd_charge_measurement,
		"evtfilter_psd_charge_measurement/O");
	DmpNtupTree->Branch(
		"evtfilter_stk_charge_measurement",
		&evtfilter_stk_charge_measurement,
		"evtfilter_stk_charge_measurement/O");
	DmpNtupTree->Branch(
		"evtfilter_xtrl_tight_cut",
		&evtfilter_xtrl_tight_cut,
		"evtfilter_xtrl_tight_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_xtrl_loose_cut",
		&evtfilter_xtrl_loose_cut,
		"evtfilter_xtrl_loose_cut/O");
	DmpNtupTree->Branch(
		"evtfilter_all_cut",
		&evtfilter_all_cut,
		"evtfilter_all_cut/O");
	// Preselection Filters
	DmpNtupTree->Branch(
		"preselection_maxelayer_cut",
		&preselection_maxelayer_cut,
		"preselection_maxelayer_cut/O");
	DmpNtupTree->Branch(
		"preselection_maxbarlayer_cut",
		&preselection_maxbarlayer_cut,
		"preselection_maxbarlayer_cut/O");
	DmpNtupTree->Branch(
		"preselection_bgotrack_cut",
		&preselection_bgotrack_cut,
		"preselection_bgotrack_cut/O");
	DmpNtupTree->Branch(
		"preselection_bgofiducial_cut",
		&preselection_bgofiducial_cut,
		"preselection_bgofiducial_cut/O");
	DmpNtupTree->Branch(
		"preselection_nbarlayer13_cut",
		&preselection_nbarlayer13_cut,
		"preselection_nbarlayer13_cut/O");
	DmpNtupTree->Branch(
		"preselection_maxrms_cut",
		&preselection_maxrms_cut,
		"preselection_maxrms_cut/O");
	DmpNtupTree->Branch(
		"preselection_trackselection_cut",
		&preselection_trackselection_cut,
		"preselection_trackselection_cut/O");
	DmpNtupTree->Branch(
		"preselection_psdstkmatch_cut",
		&preselection_psdstkmatch_cut,
		"preselection_psdstkmatch_cut/O");
	DmpNtupTree->Branch(
		"preselection_psdcharge_cut",
		&preselection_psdcharge_cut,
		"preselection_psdcharge_cut/O");
	DmpNtupTree->Branch(
		"preselection_stkcharge_cut",
		&preselection_stkcharge_cut,
		"preselection_stkcharge_cut/O");
}

void data_tuple::Fill(
	const std::shared_ptr<_tmp_filter> _filter_res,
	const std::vector<int> _stk_clusters_on_plane,
	const std::shared_ptr<_tmp_bgo> _bgo_res,
	const std::shared_ptr<_tmp_energy_data> _energy_res,
	const std::shared_ptr<DmpEvtAttitude> attitude,
	const std::shared_ptr<DmpEvtHeader> header,
	const std::shared_ptr<_tmp_nud> _nud_res)
{
	fill_trigger_info(_filter_res->evt_trigger_info);
	fill_filter_info(_filter_res->output);
	fill_preselectionfilters_info(_filter_res->preselection_output);
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
	fill_attitude_info(attitude);
	fill_psdcharge_info(_filter_res->evt_psd_charge);
	fill_stkcharge_info(_filter_res->evt_stk_charge);
	fill_classifier_info(_filter_res->evt_bgo_classifier);
	fill_nud_info(
		_nud_res->adc,
		_nud_res->total_adc,
		_nud_res->max_adc,
		_nud_res->max_channel_ID);
	second = (unsigned int)header->GetSecond();
	DmpNtupTree->Fill();
}

void data_tuple::fill_filter_info(const filter_output &output)
{
	evtfilter_out_energy_range = output.out_energy_range;
	evtfilter_evt_in_saa = output.evt_in_saa;
	evtfilter_evt_triggered = output.evt_triggered;
	evtfilter_correct_bgo_reco = output.correct_bgo_reco;
	evtfilter_good_event = output.good_event;
	evtfilter_geometric = output.geometric;
	evtfilter_BGO_fiducial = output.BGO_fiducial;
	evtfilter_BGO_fiducial_maxElayer_cut = output.BGO_fiducial_maxElayer_cut;
	evtfilter_BGO_fiducial_maxBarLayer_cut = output.BGO_fiducial_maxBarLayer_cut;
	evtfilter_BGO_fiducial_BGOTrackContainment_cut = output.BGO_fiducial_BGOTrackContainment_cut;
	evtfilter_nBarLayer13_cut = output.nBarLayer13_cut;
	evtfilter_maxRms_cut = output.maxRms_cut;
	evtfilter_track_selection_cut = output.track_selection_cut;
	evtfilter_psd_stk_match_cut = output.psd_stk_match_cut;
	evtfilter_psd_charge_cut = output.psd_charge_cut;
	evtfilter_stk_charge_cut = output.stk_charge_cut;
	evtfilter_psd_charge_measurement = output.psd_charge_measurement;
	evtfilter_stk_charge_measurement = output.stk_charge_measurement;
	evtfilter_xtrl_tight_cut = output.xtrl_tight_cut;
	evtfilter_xtrl_loose_cut = output.xtrl_loose_cut;
	evtfilter_all_cut = output.all_cut;
}

void data_tuple::fill_attitude_info(const std::shared_ptr<DmpEvtAttitude> &attitude)
{
	glat = attitude->glat;
	glon = attitude->glon;
	geo_lat = attitude->lat_geo;
	geo_lon = attitude->lon_geo;
	ra_zenith = attitude->ra_zenith;
	dec_zenith = attitude->dec_zenith;
	ra_scz = attitude->ra_scz;
	dec_scz = attitude->dec_scz;
	ra_scx = attitude->ra_scx;
	dec_scx = attitude->dec_scx;
	ra_scy = attitude->ra_scy;
	dec_scy = attitude->dec_scy;
	cutoff = attitude->verticalRigidityCutoff;
}

void data_tuple::Reset()
{
	core_reset();
	// Event time
	second = 0;
	// Attitude
	glat = -999;
	glon = -999;
	geo_lat = -999;
	geo_lon = -999;
	ra_zenith = -999;
	dec_zenith = -999;
	ra_scz = -999;
	dec_scz = -999;
	ra_scx = -999;
	dec_scx = -999;
	ra_scy = -999;
	dec_scy = -999;
	cutoff = -999;
	// Filters
	evtfilter_evt_in_saa = false;
}