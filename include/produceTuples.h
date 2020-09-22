#ifndef PRODUCETUPLES_H
#define PRODUCETUPLES_H

#include "data_loop.h"
#include "aggregate_events.h"
#include "read_sets_config_file.h"
#include "data_cuts.h"
#include "BGO_energy_cuts.h"
#include "flux.h"
#include "wtsydp.h"
#include "binning.h"
#include "charge.h"
#include "mc_ancillary.h"
#include "fill_event_histo.h"
#include "myHeader.h"
#include "DAMPE_geo_structure.h"

#include "TEfficiency.h"
#include "TTree.h"
#include "TChain.h"

#include "DmpEvtAttitude.h"
#include "DmpFilterOrbit.h"
#include "DmpEvtHeader.h"
#include "DmpIOSvc.h"
#include "DmpCore.h"

#include <memory>

struct t_variables
{
	// Trigger
	bool unbiased_trigger;
	bool mip1_trigger;
	bool mip2_trigger;
	bool HET_trigger;
	bool LET_trigger;
	bool MIP_trigger;
	bool general_trigger;

	// Event time
	unsigned int second;
	unsigned int msecond;

	// STK
	unsigned int STK_bestTrack_npoints;
	unsigned int STK_bestTrack_nholesX;
	unsigned int STK_bestTrack_nholesY;
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
	double BGOrec_slopeX;
	double BGOrec_slopeY;
	double BGOrec_interceptX;
	double BGOrec_interceptY;
	double sumRms;
	std::vector<double> fracLayer;
	std::vector<double> *fracLayer_ptr;
	double fracLast;
	double fracLast_13;
	unsigned int lastBGOLayer;
	unsigned int nBGOentries;

	// PSD
	double PSD_chargeX;
	double PSD_chargeY;
	double PSD_charge;

	// Classifiers
	double xtr;
	double xtrl;

	// Attitude
	double glat;
	double glon;
	double geo_lat;
	double geo_lon;
	double ra_zenith;
	double dec_zenith;
	double ra_scz;
	double dec_scz;
	double ra_scx;
	double dec_scx;
	double ra_scy;
	double dec_scy;
	double verticalRigidityCutoff;

	// Cuts
	bool cut_nBarLayer13;
	bool cut_maxRms;
	bool cut_track_selection;
	bool cut_psd_stk_match;
	bool cut_psd_charge;
	bool cut_stk_charge;
	bool cut_xtrl;
	unsigned int nActiveCuts;

	bool evtfilter_geometric;
	bool evtfilter_BGO_fiducial;
	bool evtfilter_all_cut;
	bool evtfilter_all_cut_no_xtrl;
	bool evtfilter_BGO_fiducial_maxElayer_cut;
	bool evtfilter_BGO_fiducial_maxBarLayer_cut;
	bool evtfilter_BGO_fiducial_BGOTrackContainment_cut;
	bool evtfilter_nBarLayer13_cut;
	bool evtfilter_maxRms_cut;
	bool evtfilter_track_selection_cut;
	bool evtfilter_psd_stk_match_cut;
	bool evtfilter_psd_charge_cut;
	bool evtfilter_stk_charge_cut;
	bool evtfilter_xtrl_cut;
	bool evtfilter_psd_charge_measurement;
	bool evtfilter_stk_charge_measurement;
};

extern void produceTuples(
	AnyOption &opt,
	const std::string inputPath,
	const bool verbose,
	const std::string wd);

template <typename linkedStruct> extern void linker(
    std::shared_ptr<linkedStruct> myStruct, 
    t_variables &vars);

extern template void linker(
	std::shared_ptr<TTree> tree,
	t_variables &vars);

extern template void linker(
	std::shared_ptr<TChain> chain,
	t_variables &vars);

#endif