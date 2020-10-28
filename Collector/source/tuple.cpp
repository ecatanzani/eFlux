#include "tuple.h"

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