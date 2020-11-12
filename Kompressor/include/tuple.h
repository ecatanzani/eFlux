#ifndef TUPLE_H
#define TUPLE_H

#include "DAMPE_geo_structure.h"

#include "TChain.h"

#include <memory>
#include <vector>

class ntuple
{
public:
	ntuple(){};
	~ntuple(){};
protected:
	void set_core_address();
	const std::vector<double> fit_shower_profile(const double costheta);
	const double compute_bgoreco_costheta();
	// Tree
	std::shared_ptr<TChain> evtch;
	// Trigger
	bool mip1_trigger;
	bool mip2_trigger;
	bool HET_trigger;
	bool LET_trigger;
	bool MIP_trigger;
	bool general_trigger;
	// STK
	int STK_bestTrack_npoints;
	int STK_bestTrack_nholesX;
	int STK_bestTrack_nholesY;
	double STK_bestTrack_slopeX;
	double STK_bestTrack_slopeY;
	double STK_bestTrack_interceptX;
	double STK_bestTrack_interceptY;
	double STK_bestTrack_costheta;
	double STK_bestTrack_phi;
	double STK_bestTrack_extr_BGO_topX;
	double STK_bestTrack_extr_BGO_topY;
	double STK_bestTrack_STK_BGO_topX_distance;
	double STK_bestTrack_STK_BGO_topY_distance;
	double STK_bestTrack_angular_distance_STK_BGO;
	double STK_chargeX;
	double STK_chargeY;
	double STK_charge;
	// BGO
	double energy;
	double energy_corr;
	std::vector<double> *eLayer = nullptr;
	double BGOrec_slopeX;
	double BGOrec_slopeY;
	double BGOrec_interceptX;
	double BGOrec_interceptY;
	double sumRms;
	std::vector<double> *rmsLayer = nullptr;
	std::vector<double> *fracLayer = nullptr;
	double fracLast;
	double fracLast_13;
	int lastBGOLayer;
	int nBGOentries;
	std::vector<double> *energy_1R_radius = nullptr;
	std::vector<double> *energy_2R_radius = nullptr;
	std::vector<double> *energy_3R_radius = nullptr;
	std::vector<double> *energy_5R_radius = nullptr;
	// PSD
	double PSD_chargeX;
	double PSD_chargeY;
	double PSD_charge;
	// NUD
	std::vector<double> *nud_adc = nullptr;
	double nud_total_adc;
	double nud_max_adc;
	int nud_max_channel_id;
	// Classifiers
	double xtr;
	double xtrl;
	// Filters
	bool evtfilter_out_energy_range;
	bool evtfilter_evt_triggered;
	bool evtfilter_correct_bgo_reco;
	bool evtfilter_good_event;
	bool evtfilter_geometric;
	bool evtfilter_BGO_fiducial;
	bool evtfilter_BGO_fiducial_maxElayer_cut;
	bool evtfilter_BGO_fiducial_maxBarLayer_cut;
	bool evtfilter_BGO_fiducial_BGOTrackContainment_cut;
	bool evtfilter_nBarLayer13_cut;
	bool evtfilter_maxRms_cut;
	bool evtfilter_track_selection_cut;
	bool evtfilter_psd_stk_match_cut;
	bool evtfilter_psd_charge_cut;
	bool evtfilter_stk_charge_cut;
	bool evtfilter_psd_charge_measurement;
	bool evtfilter_stk_charge_measurement;
	bool evtfilter_xtrl_tight_cut;
	bool evtfilter_xtrl_loose_cut;
	bool evtfilter_all_cut;
	// Preselection Cuts
	bool cut_nBarLayer13;
	bool cut_maxRms;
	bool cut_track_selection;
	bool cut_psd_stk_match;
	bool cut_psd_charge;
	bool cut_stk_charge;
	int nActiveCuts;
};

#endif