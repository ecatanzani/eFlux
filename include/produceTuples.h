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
	bool unbiased_trigger = false;
	bool mip1_trigger = false;
	bool mip2_trigger = false;
	bool HET_trigger = false;
	bool LET_trigger = false;
	bool MIP_trigger = false;
	bool general_trigger = false;

	// Event time
	unsigned int second = 0;
	unsigned int msecond = 0;

	// STK
	unsigned int STK_bestTrack_npoints = 0;
	unsigned int STK_bestTrack_nholesX = 0;
	unsigned int STK_bestTrack_nholesY = 0;
	double STK_bestTrack_slopeX = -1;
	double STK_bestTrack_slopeY = -1;
	double STK_bestTrack_interceptX = -1;
	double STK_bestTrack_interceptY = -1;
	double STK_bestTrack_costheta = -1;
	double STK_bestTrack_phi = -1;
	double STK_bestTrack_extr_BGO_topX = -1;
	double STK_bestTrack_extr_BGO_topY = -1;
	double STK_bestTrack_STK_BGO_topX_distance = -1;
	double STK_bestTrack_STK_BGO_topY_distance = -1;
	double STK_bestTrack_angular_distance_STK_BGO = -2;
	double STK_chargeX = -1;
	double STK_chargeY = -1;
	double STK_charge = -1;

	// BGO

	double energy = -1;
	double energy_corr = -1;
	double BGOrec_slopeX = -1;
	double BGOrec_slopeY = -1;
	double BGOrec_interceptX = -1;
	double BGOrec_interceptY = -1;
	double sumRms = -1;
	std::vector<double> fracLayer;
	std::vector<double> *fracLayer_ptr = nullptr;
	double fracLast = -1;
	double fracLast_13 = -1;
	unsigned int lastBGOLayer = 0;
	unsigned int nBGOentries = 0;

	// PSD
	double PSD_chargeX = -1;
	double PSD_chargeY = -1;
	double PSD_charge = -1;

	// Classifiers
	double xtr = -1;
	double xtrl = -1;

	// Attitude
	double glat = -1;
	double glon = -1;
	double geo_lat = -1;
	double geo_lon = -1;
	double ra_zenith = -1;
	double dec_zenith = -1;
	double ra_scz = -1;
	double dec_scz = -1;
	double ra_scx = -1;
	double dec_scx = -1;
	double ra_scy = -1;
	double dec_scy = -1;
	double verticalRigidityCutoff = -1;

	// Cuts
	bool cut_nBarLayer13 = false;
	bool cut_maxRms = false;
	bool cut_track_selection = false;
	bool cut_psd_stk_match = false;
	bool cut_psd_charge = false;
	bool cut_stk_charge = false;
	bool cut_xtrl = false;
	unsigned int nActiveCuts = 0;

	bool evtfilter_geometric = false;
	bool evtfilter_BGO_fiducial = false;
	bool evtfilter_all_cut = false;
	bool evtfilter_all_cut_no_xtrl = false;
	bool evtfilter_BGO_fiducial_maxElayer_cut = false;
	bool evtfilter_BGO_fiducial_maxBarLayer_cut = false;
	bool evtfilter_BGO_fiducial_BGOTrackContainment_cut = false;
	bool evtfilter_nBarLayer13_cut = false;
	bool evtfilter_maxRms_cut = false;
	bool evtfilter_track_selection_cut = false;
	bool evtfilter_psd_stk_match_cut = false;
	bool evtfilter_psd_charge_cut = false;
	bool evtfilter_stk_charge_cut = false;
	bool evtfilter_xtrl_cut = false;
	bool evtfilter_psd_charge_measurement = false;
	bool evtfilter_stk_charge_measurement = false;

	t_variables() :	fracLayer (DAMPE_bgo_nLayers, 0)
					{
						
					}
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