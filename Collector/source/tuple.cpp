#include "tuple.h"

#include "DmpStkTrack.h"

void ntuple::init(const active_cuts &acuts)
{
	// Init Tree
	DmpNtupTree = std::make_unique<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	// Set active cuts
	set_active_cuts(acuts);
	// Branch tree
	branch_tree();
}

void ntuple::set_active_cuts(const active_cuts &acuts)
{
	cut_nBarLayer13 = acuts.nBarLayer13;
	cut_maxRms = acuts.maxRms;
	cut_track_selection = acuts.track_selection;
	cut_psd_stk_match = acuts.psd_stk_match;
	cut_psd_charge = acuts.psd_charge;
	cut_stk_charge = acuts.stk_charge;
	nActiveCuts = acuts.nActiveCuts;
}

void ntuple::branch_tree()
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
	// Event time
	DmpNtupTree->Branch(
		"second",
		&second,
		"second/i");
	// STK
	DmpNtupTree->Branch(
		"STK_bestTrack_npoints",
		&STK_bestTrack_npoints,
		"STK_bestTrack_npoints/i");
	DmpNtupTree->Branch(
		"STK_bestTrack_nholesX",
		&STK_bestTrack_nholesX,
		"STK_bestTrack_nholesX/i");
	DmpNtupTree->Branch(
		"STK_bestTrack_nholesY",
		&STK_bestTrack_nholesY,
		"STK_bestTrack_nholesY/i");
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
		"sumRms",
		&sumRms,
		"sumRms/D");
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
		"lastBGOLayer/i");
	DmpNtupTree->Branch(
		"nBGOentries",
		&nBGOentries,
		"nBGOentries/i");
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
		"nActiveCuts/i");
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
}

void ntuple::Fill(
	const filter_output &output,
	const std::shared_ptr<DmpEvtAttitude> &attitude,
	const best_track &event_best_track,
	const double raw_energy,
	const double corr_energy,
	const std::vector<double> &bgoRec_slope,
	const std::vector<double> &bgoRec_intercept,
	const double sumRMS,
	const std::vector<double> bgo_fracLayer,
	const double lastFracLayer,
	const double frac_layer_13,
	const unsigned int last_bgo_layer,
	const unsigned int bgo_entries,
	const psd_charge &extracted_psd_charge,
	const stk_charge &extracted_stk_charge,
	const bgo_classifiers &classifier,
	const trigger_info &evt_trigger)
{
	fill_trigger_info(evt_trigger);
	fill_filter_info(output);
	fill_attitude_info(attitude);
	fill_stk_info(event_best_track);
	fill_bgo_info(
		raw_energy,
		corr_energy,
		bgoRec_slope,
		bgoRec_intercept,
		sumRMS,
		bgo_fracLayer,
		lastFracLayer,
		frac_layer_13,
		last_bgo_layer,
		bgo_entries);
	fill_psdcharge_info(extracted_psd_charge);
	fill_stkcharge_info(extracted_stk_charge);
	fill_classifier_info(classifier);

	DmpNtupTree->Fill();
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

void ntuple::fill_filter_info(const filter_output &output)
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

void ntuple::fill_attitude_info(const std::shared_ptr<DmpEvtAttitude> &attitude)
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
	const std::vector<double> &bgoRec_slope,
	const std::vector<double> &bgoRec_intercept,
	const double sumRMS,
	const std::vector<double> bgo_fracLayer,
	const double lastFracLayer,
	const double frac_layer_13,
	const unsigned int last_bgo_layer,
	const unsigned int bgo_entries)
{
	if (evtfilter_correct_bgo_reco)
	{
		energy = raw_energy;
		energy_corr = corr_energy;
		BGOrec_slopeX = bgoRec_slope[0];
		BGOrec_slopeY = bgoRec_slope[1];
		BGOrec_interceptX = bgoRec_intercept[0];
		BGOrec_interceptY = bgoRec_intercept[1];
		sumRms = sumRMS;
		fracLayer = bgo_fracLayer;
		fracLast = lastFracLayer;
		fracLast_13 = frac_layer_13;
		lastBGOLayer = last_bgo_layer;
		nBGOentries = bgo_entries;
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