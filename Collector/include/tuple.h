#ifndef TUPLE_H
#define TUPLE_H

#include "config.h"
#include "DmpFilterContainer.h"
#include "DAMPE_geo_structure.h"

#include "DmpEvtAttitude.h"

#include "TTree.h"
#include "TFile.h"

#include <memory>
#include <vector>

class ntuple
{
public:
	ntuple(const active_cuts &acuts)	:	fracLayer(DAMPE_bgo_nLayers, -999)
	{
		init(acuts);
	};
	~ntuple(){};
	void Fill(
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
		const trigger_info &evt_trigger);
	void Write(TFile &outfile);

private:
	void init(const active_cuts &acuts);
	void set_active_cuts(const active_cuts &acuts);
	void branch_tree();
	void fill_trigger_info(const trigger_info &evt_trigger);
	void fill_filter_info(const filter_output &output);
	void fill_attitude_info(const std::shared_ptr<DmpEvtAttitude> &attitude);
	void fill_stk_info(const best_track &event_best_track);
	void fill_bgo_info(
		const double raw_energy,
		const double corr_energy,
		const std::vector<double> &bgoRec_slope,
		const std::vector<double> &bgoRec_intercept,
		const double sumRMS,
		const std::vector<double> bgo_fracLayer,
		const double lastFracLayer,
		const double frac_layer_13,
		const unsigned int last_bgo_layer,
		const unsigned int bgo_entries);
	void fill_psdcharge_info(const psd_charge &extracted_psd_charge);
	void fill_stkcharge_info(const stk_charge &extracted_stk_charge);
	void fill_classifier_info(const bgo_classifiers &classifier);

	// Tree
	std::unique_ptr<TTree> DmpNtupTree;
	// Trigger
	bool mip1_trigger = false;
	bool mip2_trigger = false;
	bool HET_trigger = false;
	bool LET_trigger = false;
	bool MIP_trigger = false;
	bool general_trigger = false;
	// Event time
	unsigned int second = 0;
	// STK
	unsigned int STK_bestTrack_npoints = 0;
	unsigned int STK_bestTrack_nholesX = 0;
	unsigned int STK_bestTrack_nholesY = 0;
	double STK_bestTrack_slopeX = -999;
	double STK_bestTrack_slopeY = -999;
	double STK_bestTrack_interceptX = -999;
	double STK_bestTrack_interceptY = -999;
	double STK_bestTrack_costheta = -999;
	double STK_bestTrack_phi = -999;
	double STK_bestTrack_extr_BGO_topX = -999;
	double STK_bestTrack_extr_BGO_topY = -999;
	double STK_bestTrack_STK_BGO_topX_distance = -999;
	double STK_bestTrack_STK_BGO_topY_distance = -999;
	double STK_bestTrack_angular_distance_STK_BGO = -999;
	double STK_chargeX = -999;
	double STK_chargeY = -999;
	double STK_charge = -999;
	// BGO
	double energy = -999;
	double energy_corr = -999;
	double BGOrec_slopeX = -999;
	double BGOrec_slopeY = -999;
	double BGOrec_interceptX = -999;
	double BGOrec_interceptY = -999;
	double sumRms = -999;
	std::vector<double> fracLayer;
	double fracLast = -999;
	double fracLast_13 = -999;
	unsigned int lastBGOLayer = 0;
	unsigned int nBGOentries = 0;
	// PSD
	double PSD_chargeX = -999;
	double PSD_chargeY = -999;
	double PSD_charge = -999;
	// Classifiers
	double xtr = -999;
	double xtrl = -999;
	// Attitude
	double glat = -999;
	double glon = -999;
	double geo_lat = -999;
	double geo_lon = -999;
	double ra_zenith = -999;
	double dec_zenith = -999;
	double ra_scz = -999;
	double dec_scz = -999;
	double ra_scx = -999;
	double dec_scx = -999;
	double ra_scy = -999;
	double dec_scy = -999;
	double cutoff = -999;
	// Preselection Cuts
	bool cut_nBarLayer13 = false;
	bool cut_maxRms = false;
	bool cut_track_selection = false;
	bool cut_psd_stk_match = false;
	bool cut_psd_charge = false;
	bool cut_stk_charge = false;
	unsigned int nActiveCuts = 0;
	// Filters
	bool evtfilter_out_energy_range = false;
	bool evtfilter_evt_in_saa = false;
	bool evtfilter_evt_triggered = false;
	bool evtfilter_correct_bgo_reco = false;
	bool evtfilter_good_event = false;
	bool evtfilter_geometric = false;
	bool evtfilter_BGO_fiducial = false;
	bool evtfilter_BGO_fiducial_maxElayer_cut = false;
	bool evtfilter_BGO_fiducial_maxBarLayer_cut = false;
	bool evtfilter_BGO_fiducial_BGOTrackContainment_cut = false;
	bool evtfilter_nBarLayer13_cut = false;
	bool evtfilter_maxRms_cut = false;
	bool evtfilter_track_selection_cut = false;
	bool evtfilter_psd_stk_match_cut = false;
	bool evtfilter_psd_charge_cut = false;
	bool evtfilter_stk_charge_cut = false;
	bool evtfilter_psd_charge_measurement = false;
	bool evtfilter_stk_charge_measurement = false;
	bool evtfilter_xtrl_tight_cut = false;
	bool evtfilter_xtrl_loose_cut = false;
	bool evtfilter_all_cut = false;
};

#endif