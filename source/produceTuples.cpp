#include "produceTuples.h"

inline void updateProcessStatus(const int evIdx, int &kStep, const int nevents)
{
	auto percentage = ((evIdx + 1) / (double)nevents) * 100;
	if (floor(percentage) != 0 && ((int)floor(percentage) % kStep) == 0)
	{
		std::cout << "\n"
				  << floor(percentage) << " %\t | \tProcessed " << evIdx + 1 << " events / " << nevents;
		kStep += 10;
	}
}

inline void branchTree(
	std::shared_ptr<TTree> tree, 
	t_variables &vars)
{
	tree->Branch("unbiased_trigger", &vars.unbiased_trigger, "unbiased_trigger/O");
	tree->Branch("mip1_trigger", &vars.mip1_trigger, "mip1_trigger/O");
	tree->Branch("mip2_trigger", &vars.mip2_trigger, "mip2_trigger/O");
	tree->Branch("HET_trigger", &vars.HET_trigger, "HET_trigger/O");
	tree->Branch("LET_trigger", &vars.LET_trigger, "LET_trigger/O");
	tree->Branch("MIP_trigger", &vars.MIP_trigger, "MIP_trigger/O");
	tree->Branch("general_trigger", &vars.general_trigger, "general_trigger/O");
	tree->Branch("second", &vars.second, "second/i");
	tree->Branch("msecond", &vars.msecond, "msecond/i");
	tree->Branch("STK_bestTrack_npoints", &vars.STK_bestTrack_npoints, "STK_bestTrack_npoints/i");
	tree->Branch("STK_bestTrack_nholesX", &vars.STK_bestTrack_nholesX, "STK_bestTrack_nholesX/i");
	tree->Branch("STK_bestTrack_nholesY", &vars.STK_bestTrack_nholesY, "STK_bestTrack_nholesY/i");
	tree->Branch("STK_bestTrack_slopeX", &vars.STK_bestTrack_slopeX, "STK_bestTrack_slopeX/D");
	tree->Branch("STK_bestTrack_slopeY", &vars.STK_bestTrack_slopeY, "STK_bestTrack_slopeY/D");
	tree->Branch("STK_bestTrack_interceptX", &vars.STK_bestTrack_interceptX, "STK_bestTrack_interceptX/D");
	tree->Branch("STK_bestTrack_interceptY", &vars.STK_bestTrack_interceptY, "STK_bestTrack_interceptY/D");
	tree->Branch("STK_bestTrack_costheta", &vars.STK_bestTrack_costheta, "STK_bestTrack_costheta/D");
	tree->Branch("STK_bestTrack_phi", &vars.STK_bestTrack_phi, "STK_bestTrack_phi/D");
	tree->Branch("STK_bestTrack_extr_BGO_topX", &vars.STK_bestTrack_extr_BGO_topX, "STK_bestTrack_extr_BGO_topX/D");
	tree->Branch("STK_bestTrack_extr_BGO_topY", &vars.STK_bestTrack_extr_BGO_topY, "STK_bestTrack_extr_BGO_topY/D");
	tree->Branch("STK_bestTrack_STK_BGO_topX_distance", &vars.STK_bestTrack_STK_BGO_topX_distance, "STK_bestTrack_STK_BGO_topX_distance/D");
	tree->Branch("STK_bestTrack_STK_BGO_topY_distance", &vars.STK_bestTrack_STK_BGO_topY_distance, "STK_bestTrack_STK_BGO_topY_distance/D");
	tree->Branch("STK_bestTrack_angular_distance_STK_BGO", &vars.STK_bestTrack_angular_distance_STK_BGO, "STK_bestTrack_angular_distance_STK_BGO/D");
	tree->Branch("STK_chargeX", &vars.STK_chargeX, "STK_chargeX/D");
	tree->Branch("STK_chargeY", &vars.STK_chargeY, "STK_chargeY/D");
	tree->Branch("STK_charge", &vars.STK_charge, "STK_charge/D");
	tree->Branch("energy", &vars.energy, "energy/D");
	tree->Branch("energy_corr", &vars.energy_corr, "energy_corr/D");
	tree->Branch("BGOrec_slopeX", &vars.BGOrec_slopeX, "BGOrec_slopeX/D");
	tree->Branch("BGOrec_slopeY", &vars.BGOrec_slopeY, "BGOrec_slopeY/D");
	tree->Branch("BGOrec_interceptX", &vars.BGOrec_interceptX, "BGOrec_interceptX/D");
	tree->Branch("BGOrec_interceptY", &vars.BGOrec_interceptY, "BGOrec_interceptY/D");
	tree->Branch("sumRms", &vars.sumRms, "sumRms/D");
	tree->Branch("fracLayer", &vars.fracLayer);
	tree->Branch("fracLast", &vars.fracLast, "fracLast/D");
	tree->Branch("fracLast_13", &vars.fracLast_13, "fracLast_13/D");
	tree->Branch("lastBGOLayer", &vars.lastBGOLayer, "lastBGOLayer/i");
	tree->Branch("nBGOentries", &vars.nBGOentries, "nBGOentries/i");
	tree->Branch("PSD_chargeX", &vars.PSD_chargeX, "PSD_chargeX/D");
	tree->Branch("PSD_chargeY", &vars.PSD_chargeY, "PSD_chargeY/D");
	tree->Branch("PSD_charge", &vars.PSD_charge, "PSD_charge/D");
	tree->Branch("xtr", &vars.xtr, "xtr/D");
	tree->Branch("xtrl", &vars.xtrl, "xtrl/D");
	tree->Branch("glat", &vars.glat, "glat/D");
	tree->Branch("glon", &vars.glon, "glon/D");
	tree->Branch("cut_nBarLayer13", &vars.cut_nBarLayer13, "cut_nBarLayer13/O");
	tree->Branch("cut_maxRms", &vars.cut_maxRms, "cut_maxRms/O");
	tree->Branch("cut_track_selection", &vars.cut_track_selection, "cut_track_selection/O");
	tree->Branch("cut_psd_stk_match", &vars.cut_psd_stk_match, "cut_psd_stk_match/O");
	tree->Branch("cut_psd_charge", &vars.cut_psd_charge, "cut_psd_charge/O");
	tree->Branch("cut_stk_charge", &vars.cut_stk_charge, "cut_stk_charge/O");
	tree->Branch("cut_xtrl", &vars.cut_xtrl, "cut_xtrl/O");
	tree->Branch("nActiveCuts", &vars.nActiveCuts, "nActiveCuts/i");
	tree->Branch("evtfilter_geometric", &vars.evtfilter_geometric, "evtfilter_geometric/O");
	tree->Branch("evtfilter_BGO_fiducial", &vars.evtfilter_BGO_fiducial, "evtfilter_BGO_fiducial/O");
	tree->Branch("evtfilter_all_cut", &vars.evtfilter_all_cut, "evtfilter_all_cut/O");
	tree->Branch("evtfilter_all_cut_no_xtrl", &vars.evtfilter_all_cut_no_xtrl, "evtfilter_all_cut_no_xtrl/O");
	tree->Branch("evtfilter_BGO_fiducial_maxElayer_cut", &vars.evtfilter_BGO_fiducial_maxElayer_cut, "evtfilter_BGO_fiducial_maxElayer_cut/O");
	tree->Branch("evtfilter_BGO_fiducial_maxBarLayer_cut", &vars.evtfilter_BGO_fiducial_maxBarLayer_cut, "evtfilter_BGO_fiducial_maxBarLayer_cut/O");
	tree->Branch("evtfilter_BGO_fiducial_BGOTrackContainment_cut", &vars.evtfilter_BGO_fiducial_BGOTrackContainment_cut, "evtfilter_BGO_fiducial_BGOTrackContainment_cut/O");
	tree->Branch("evtfilter_nBarLayer13_cut", &vars.evtfilter_nBarLayer13_cut, "evtfilter_nBarLayer13_cut/O");
	tree->Branch("evtfilter_maxRms_cut", &vars.evtfilter_maxRms_cut, "evtfilter_maxRms_cut/O");
	tree->Branch("evtfilter_track_selection_cut", &vars.evtfilter_track_selection_cut, "evtfilter_track_selection_cut/O");
	tree->Branch("evtfilter_psd_stk_match_cut", &vars.evtfilter_psd_stk_match_cut, "evtfilter_psd_stk_match_cut/O");
	tree->Branch("evtfilter_psd_charge_cut", &vars.evtfilter_psd_charge_cut, "evtfilter_psd_charge_cut/O");
	tree->Branch("evtfilter_stk_charge_cut", &vars.evtfilter_stk_charge_cut, "evtfilter_stk_charge_cut/O");
	tree->Branch("evtfilter_xtrl_cut", &vars.evtfilter_xtrl_cut, "evtfilter_xtrl_cut/O");
	tree->Branch("evtfilter_psd_charge_measurement", &vars.evtfilter_psd_charge_measurement, "evtfilter_psd_charge_measurement/O");
	tree->Branch("evtfilter_stk_charge_measurement", &vars.evtfilter_stk_charge_measurement, "evtfilter_stk_charge_measurement/O");
}

template <typename linkedStruct> void linker(
    std::shared_ptr<linkedStruct> myStruct, 
    t_variables &vars)
{
	myStruct->SetBranchAddress("unbiased_trigger", &vars.unbiased_trigger);
	myStruct->SetBranchAddress("mip1_trigger", &vars.mip1_trigger);
	myStruct->SetBranchAddress("mip2_trigger", &vars.mip2_trigger);
	myStruct->SetBranchAddress("HET_trigger", &vars.HET_trigger);
	myStruct->SetBranchAddress("LET_trigger", &vars.LET_trigger);
	myStruct->SetBranchAddress("MIP_trigger", &vars.MIP_trigger);
	myStruct->SetBranchAddress("general_trigger", &vars.general_trigger);
	myStruct->SetBranchAddress("second", &vars.second);
	myStruct->SetBranchAddress("msecond", &vars.msecond);
	myStruct->SetBranchAddress("STK_bestTrack_npoints", &vars.STK_bestTrack_npoints);
	myStruct->SetBranchAddress("STK_bestTrack_nholesX", &vars.STK_bestTrack_nholesX);
	myStruct->SetBranchAddress("STK_bestTrack_nholesY", &vars.STK_bestTrack_nholesY);
	myStruct->SetBranchAddress("STK_bestTrack_slopeX", &vars.STK_bestTrack_slopeX);
	myStruct->SetBranchAddress("STK_bestTrack_slopeY", &vars.STK_bestTrack_slopeY);
	myStruct->SetBranchAddress("STK_bestTrack_interceptX", &vars.STK_bestTrack_interceptX);
	myStruct->SetBranchAddress("STK_bestTrack_interceptY", &vars.STK_bestTrack_interceptY);
	myStruct->SetBranchAddress("STK_bestTrack_costheta", &vars.STK_bestTrack_costheta);
	myStruct->SetBranchAddress("STK_bestTrack_phi", &vars.STK_bestTrack_phi);
	myStruct->SetBranchAddress("STK_bestTrack_extr_BGO_topX", &vars.STK_bestTrack_extr_BGO_topX);
	myStruct->SetBranchAddress("STK_bestTrack_extr_BGO_topY", &vars.STK_bestTrack_extr_BGO_topY);
	myStruct->SetBranchAddress("STK_bestTrack_STK_BGO_topX_distance", &vars.STK_bestTrack_STK_BGO_topX_distance);
	myStruct->SetBranchAddress("STK_bestTrack_STK_BGO_topY_distance", &vars.STK_bestTrack_STK_BGO_topY_distance);
	myStruct->SetBranchAddress("STK_bestTrack_angular_distance_STK_BGO", &vars.STK_bestTrack_angular_distance_STK_BGO);
	myStruct->SetBranchAddress("STK_chargeX", &vars.STK_chargeX);
	myStruct->SetBranchAddress("STK_chargeY", &vars.STK_chargeY);
	myStruct->SetBranchAddress("STK_charge", &vars.STK_charge);
	myStruct->SetBranchAddress("energy", &vars.energy);
	myStruct->SetBranchAddress("energy_corr", &vars.energy_corr);
	myStruct->SetBranchAddress("BGOrec_slopeX", &vars.BGOrec_slopeX);
	myStruct->SetBranchAddress("BGOrec_slopeY", &vars.BGOrec_slopeY);
	myStruct->SetBranchAddress("BGOrec_interceptX", &vars.BGOrec_interceptX);
	myStruct->SetBranchAddress("BGOrec_interceptY", &vars.BGOrec_interceptY);
	myStruct->SetBranchAddress("sumRms", &vars.sumRms);
	myStruct->SetBranchAddress("fracLayer", &vars.fracLayer_ptr);
	myStruct->SetBranchAddress("fracLast", &vars.fracLast);
	myStruct->SetBranchAddress("fracLast_13", &vars.fracLast_13);
	myStruct->SetBranchAddress("lastBGOLayer", &vars.lastBGOLayer);
	myStruct->SetBranchAddress("nBGOentries", &vars.nBGOentries);
	myStruct->SetBranchAddress("PSD_chargeX", &vars.PSD_chargeX);
	myStruct->SetBranchAddress("PSD_chargeY", &vars.PSD_chargeY);
	myStruct->SetBranchAddress("PSD_charge", &vars.PSD_charge);
	myStruct->SetBranchAddress("xtr", &vars.xtr);
	myStruct->SetBranchAddress("xtrl", &vars.xtrl);
	myStruct->SetBranchAddress("glat", &vars.glat);
	myStruct->SetBranchAddress("glon", &vars.glon);
	myStruct->SetBranchAddress("cut_nBarLayer13", &vars.cut_nBarLayer13);
	myStruct->SetBranchAddress("cut_maxRms", &vars.cut_maxRms);
	myStruct->SetBranchAddress("cut_track_selection", &vars.cut_track_selection);
	myStruct->SetBranchAddress("cut_psd_stk_match", &vars.cut_psd_stk_match);
	myStruct->SetBranchAddress("cut_psd_charge", &vars.cut_psd_charge);
	myStruct->SetBranchAddress("cut_stk_charge", &vars.cut_stk_charge);
	myStruct->SetBranchAddress("cut_xtrl", &vars.cut_xtrl);
	myStruct->SetBranchAddress("nActiveCuts", &vars.nActiveCuts);
	myStruct->SetBranchAddress("evtfilter_geometric", &vars.evtfilter_geometric);
	myStruct->SetBranchAddress("evtfilter_BGO_fiducial", &vars.evtfilter_BGO_fiducial);
	myStruct->SetBranchAddress("evtfilter_all_cut", &vars.evtfilter_all_cut);
	myStruct->SetBranchAddress("evtfilter_all_cut_no_xtrl", &vars.evtfilter_all_cut_no_xtrl);
	myStruct->SetBranchAddress("evtfilter_BGO_fiducial_maxElayer_cut", &vars.evtfilter_BGO_fiducial_maxElayer_cut);
	myStruct->SetBranchAddress("evtfilter_BGO_fiducial_maxBarLayer_cut", &vars.evtfilter_BGO_fiducial_maxBarLayer_cut);
	myStruct->SetBranchAddress("evtfilter_BGO_fiducial_BGOTrackContainment_cut", &vars.evtfilter_BGO_fiducial_BGOTrackContainment_cut);
	myStruct->SetBranchAddress("evtfilter_nBarLayer13_cut", &vars.evtfilter_nBarLayer13_cut);
	myStruct->SetBranchAddress("evtfilter_maxRms_cut", &vars.evtfilter_maxRms_cut);
	myStruct->SetBranchAddress("evtfilter_track_selection_cut", &vars.evtfilter_track_selection_cut);
	myStruct->SetBranchAddress("evtfilter_psd_stk_match_cut", &vars.evtfilter_psd_stk_match_cut);
	myStruct->SetBranchAddress("evtfilter_psd_charge_cut", &vars.evtfilter_psd_charge_cut);
	myStruct->SetBranchAddress("evtfilter_stk_charge_cut", &vars.evtfilter_stk_charge_cut);
	myStruct->SetBranchAddress("evtfilter_xtrl_cut", &vars.evtfilter_xtrl_cut);
	myStruct->SetBranchAddress("evtfilter_psd_charge_measurement", &vars.evtfilter_psd_charge_measurement);
	myStruct->SetBranchAddress("evtfilter_stk_charge_measurement", &vars.evtfilter_stk_charge_measurement);
}


inline void fill_STK_bestTrack_info(best_track event_best_track, t_variables &vars)
{
	vars.STK_bestTrack_npoints = event_best_track.n_points;
	vars.STK_bestTrack_nholesX = event_best_track.n_holes[0];
	vars.STK_bestTrack_nholesY = event_best_track.n_holes[1];
	vars.STK_bestTrack_slopeX = event_best_track.track_slope[0];
	vars.STK_bestTrack_slopeY = event_best_track.track_slope[1];
	vars.STK_bestTrack_interceptX = event_best_track.track_intercept[0];
	vars.STK_bestTrack_interceptY = event_best_track.track_intercept[1];
	vars.STK_bestTrack_costheta = event_best_track.myBestTrack.getDirection().CosTheta();
	vars.STK_bestTrack_phi = event_best_track.myBestTrack.getDirection().Phi();
	vars.STK_bestTrack_extr_BGO_topX = event_best_track.extr_BGO_topX;
	vars.STK_bestTrack_extr_BGO_topY = event_best_track.extr_BGO_topY;
	vars.STK_bestTrack_STK_BGO_topX_distance = event_best_track.STK_BGO_topX_distance;
	vars.STK_bestTrack_STK_BGO_topY_distance = event_best_track.STK_BGO_topY_distance;
	vars.STK_bestTrack_angular_distance_STK_BGO = event_best_track.angular_distance_STK_BGO;
}

inline void fill_attitude_info(const std::shared_ptr<DmpEvtAttitude> attitude, t_variables &vars)
{
	vars.glat = attitude->glat;
	vars.glon = attitude->glon;
	vars.geo_lat = attitude->lat_geo;
	vars.geo_lon = attitude->lon_geo;
	vars.ra_zenith = attitude->ra_zenith;
	vars.dec_zenith = attitude->dec_zenith;
	vars.ra_scz = attitude->ra_scz;
	vars.dec_scz = attitude->dec_scz;
	vars.ra_scx = attitude->ra_scx;
	vars.dec_scx = attitude->dec_scx;
	vars.ra_scy = attitude->ra_scy;
	vars.dec_scy = attitude->dec_scy;
	vars.verticalRigidityCutoff = attitude->verticalRigidityCutoff;
}

inline void setActiveCuts(t_variables &vars, const data_active_cuts active_cuts)
{
	vars.cut_nBarLayer13 = active_cuts.nBarLayer13;
	vars.cut_maxRms = active_cuts.maxRms;
	vars.cut_track_selection = active_cuts.track_selection;
	vars.cut_psd_stk_match = active_cuts.psd_stk_match;
	vars.cut_psd_charge = active_cuts.psd_charge;
	vars.cut_stk_charge = active_cuts.stk_charge;
	vars.cut_xtrl = active_cuts.xtrl;
	vars.nActiveCuts = active_cuts.nActiveCuts;
}

inline void setEvtFilter(t_variables &vars, const event_filter filter)
{
	vars.evtfilter_geometric = filter.geometric;
	vars.evtfilter_BGO_fiducial = filter.BGO_fiducial;
	vars.evtfilter_all_cut = filter.all_cut;
	vars.evtfilter_all_cut_no_xtrl = filter.all_cut_no_xtrl;
	vars.evtfilter_BGO_fiducial_maxElayer_cut = filter.BGO_fiducial_maxElayer_cut;
	vars.evtfilter_BGO_fiducial_maxBarLayer_cut = filter.BGO_fiducial_maxBarLayer_cut;
	vars.evtfilter_BGO_fiducial_BGOTrackContainment_cut = filter.BGO_fiducial_BGOTrackContainment_cut;
	vars.evtfilter_nBarLayer13_cut = filter.nBarLayer13_cut;
	vars.evtfilter_maxRms_cut = filter.maxRms_cut;
	vars.evtfilter_track_selection_cut = filter.track_selection_cut;
	vars.evtfilter_psd_stk_match_cut = filter.psd_stk_match_cut;
	vars.evtfilter_psd_charge_cut = filter.psd_charge_cut;
	vars.evtfilter_stk_charge_cut = filter.stk_charge_cut;
	vars.evtfilter_xtrl_cut = filter.xtrl_cut;
	vars.evtfilter_psd_charge_measurement = filter.psd_charge_measurement;
	vars.evtfilter_stk_charge_measurement = filter.stk_charge_measurement;
}

void produceTuples(
	AnyOption &opt,
	const std::string inputPath,
	const bool verbose,
	const std::string wd)
{	
	int data_year;
	int data_month;

	auto dmpch = aggregateTupleDataEventsTChain(
		inputPath, 
		data_year, 
		data_month, 
		verbose);

	// Register Header container
	std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
	dmpch->SetBranchAddress("EventHeader", &evt_header);

	// Register BGO constainer
	std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
	dmpch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

	// Register BGO REC constainer
	std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
	dmpch->SetBranchAddress("DmpEvtBgoRec", &bgorec);

	// Register STK container
	std::shared_ptr<TClonesArray> stkclusters = std::make_shared<TClonesArray>("DmpStkSiCluster");
	dmpch->SetBranchAddress("StkClusterCollection", &stkclusters);

	// Check if STK tracks collection exists
	bool fStkKalmanTracksFound = false;
	for (int brIdx = 0; brIdx < dmpch->GetListOfBranches()->GetEntries(); ++brIdx)
		if (strcmp(dmpch->GetListOfBranches()->At(brIdx)->GetName(), "StkKalmanTracks"))
		{
			fStkKalmanTracksFound = true;
			break;
		}

	// Register STK tracks collection
	std::shared_ptr<TClonesArray> stktracks = std::make_shared<TClonesArray>("DmpStkTrack", 200);
	if (fStkKalmanTracksFound)
		dmpch->SetBranchAddress("StkKalmanTracks", &stktracks);

	// Register PSD container
	std::shared_ptr<DmpEvtPsdHits> psdhits = std::make_shared<DmpEvtPsdHits>();
	dmpch->SetBranchAddress("DmpPsdHits", &psdhits);

	// Register attitude container
	std::shared_ptr<DmpEvtAttitude> attitude = std::make_shared<DmpEvtAttitude>();
	dmpch->SetBranchAddress("EvtAttitudeContainer", &attitude);

	// Orbit filter
	// Set gIOSvc
	gIOSvc->Set("OutData/NoOutput", "True");
	gIOSvc->Initialize();
	// Create orbit filter
	std::unique_ptr<DmpFilterOrbit> pFilter = std::make_unique<DmpFilterOrbit>("EventHeader");
	// Activate orbit filter
	pFilter->ActiveMe(); // Call this function to calculate SAA through House Keeping Data

	// Create flux cuts struct
	cuts_conf flux_cuts;
	// Create active cuts struct
	data_active_cuts active_cuts;

	// Load structs reading config file
	auto logEBins = load_flux_struct(
		flux_cuts,
		active_cuts,
		wd);

	// Load MC statistics structure
	data_statistics data_selection;

	// Print filter status
	if (verbose)
		print_filter_status(active_cuts);

	// Event loop
	auto nevents = dmpch->GetEntries();
	if (verbose)
		std::cout << "Total number of events: " << nevents << "\n\n";

	double _GeV = 0.001;
	int kStep = 10;

	// Create TTrees
	std::shared_ptr<TTree> DmpNtupTree_20_100 = std::make_shared<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	std::shared_ptr<TTree> DmpNtupTree_100_250 = std::make_shared<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	std::shared_ptr<TTree> DmpNtupTree_250_500 = std::make_shared<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	std::shared_ptr<TTree> DmpNtupTree_500_1000 = std::make_shared<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	std::shared_ptr<TTree> DmpNtupTree_1000_5000 = std::make_shared<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");
	std::shared_ptr<TTree> DmpNtupTree_5000_10000 = std::make_shared<TTree>("DmpEvtNtup", "DAMPE Event nTuple Tree");

	t_variables vars;

	setActiveCuts(vars, active_cuts);

	branchTree(DmpNtupTree_20_100, vars);
	branchTree(DmpNtupTree_100_250, vars);
	branchTree(DmpNtupTree_250_500, vars);
	branchTree(DmpNtupTree_500_1000, vars);
	branchTree(DmpNtupTree_1000_5000, vars);
	branchTree(DmpNtupTree_5000_10000, vars);

	for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
	{
		// Get chain event
		dmpch->GetEvent(evIdx);

		// Update event counter
		++data_selection.event_counter;

		vars.second = evt_header->GetSecond();
		vars.msecond = evt_header->GetMillisecond();

		if (pFilter->IsInSAA(vars.second))
		{
			continue;
			++data_selection.events_in_saa;
		}

		// Event printout
		if (verbose)
			updateProcessStatus(evIdx, kStep, nevents);

		// Get event total energy
		double bgoTotalE_raw = bgorec->GetTotalEnergy();
		double bgoTotalE_corr = bgorec->GetElectronEcor();

		// Don't accept events outside the selected energy window
		if (bgoTotalE_raw * _GeV < flux_cuts.min_event_energy || bgoTotalE_raw * _GeV > flux_cuts.max_event_energy)
		{
			++data_selection.events_out_range;
			continue;
		}

		++data_selection.events_in_range;

		// Read trigger status
		// For MC events triggers 1 and 2 are always disabled
		bool unbiased_tr = evt_header->GeneratedTrigger(0) && evt_header->EnabledTrigger(0);
		bool mip1_tr = evt_header->GeneratedTrigger(1);
		bool mip2_tr = evt_header->GeneratedTrigger(2);
		bool HET_tr = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
		bool LET_tr = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
		bool MIP_tr = mip1_tr || mip2_tr;

		bool general_trigger = MIP_tr || HET_tr || LET_tr;

		vars.unbiased_trigger = unbiased_tr;
		vars.mip1_trigger = mip1_tr;
		vars.mip2_trigger = mip2_tr;
		vars.HET_trigger = HET_tr;
		vars.LET_trigger = LET_tr;
		vars.MIP_trigger = MIP_tr;
		vars.general_trigger = general_trigger;

		// Check if the event has been triggered or not
		if (general_trigger)
		{
			++data_selection.triggered_events;
			if (!checkBGOreco_data(bgorec))
				continue;
		}
		else
			continue;

		// Load BGO event class
		DmpBgoContainer bgoVault;

		// Load PSD event class
		DmpPsdContainer psdVault;

		bgoVault.scanBGOHits(
			bgohits,
			bgoTotalE_raw,
			flux_cuts);

		psdVault.scanPSDHits(
			psdhits,
			flux_cuts);

		// Create best_track struct
		best_track event_best_track;

		// Create PSD-STK match struct
		psd_cluster_match clu_matching;

		// Create PSD charge struct
		psd_charge extracted_psd_charge;

		// Create STK charge struct
		stk_charge extracted_stk_charge;

		// Event filter
		event_filter filter;
		if (filter_this_data_event(
				filter,
				bgorec,
				bgohits,
				flux_cuts,
				bgoTotalE_raw,
				bgoVault,
				psdVault,
				event_best_track,
				clu_matching,
				extracted_psd_charge,
				extracted_stk_charge,
				stkclusters,
				stktracks,
				active_cuts))
			++data_selection.selected_events;
		setEvtFilter(vars, filter);

		if (filter.all_cut_no_xtrl)
		{
			// Filling STK best track info
			fill_STK_bestTrack_info(event_best_track, vars);

			// Filling BGO info
			vars.energy = bgoTotalE_raw;
			vars.energy_corr = bgoTotalE_corr;
			vars.BGOrec_slopeX = bgorec->GetSlopeXZ();
			vars.BGOrec_slopeY = bgorec->GetSlopeYZ();
			vars.BGOrec_interceptX = bgorec->GetInterceptXZ();
			vars.BGOrec_interceptY = bgorec->GetInterceptYZ();
			vars.sumRms = bgoVault.GetSumRMS();
			vars.fracLayer = bgoVault.GetFracLayer();
			vars.fracLast = bgoVault.GetLastFFracLayer();
			vars.fracLast_13 = bgoVault.GetSingleFracLayer(13);
			vars.lastBGOLayer = bgoVault.GetFracIdxLastLayer();
			vars.nBGOentries = bgoVault.GetNhits();

			// Filling charges info
			vars.STK_chargeX = extracted_stk_charge.chargeX;
			vars.STK_chargeY = extracted_stk_charge.chargeY;
			vars.STK_charge = 0.5 * (extracted_stk_charge.chargeX + extracted_stk_charge.chargeY);
			vars.PSD_chargeX = extracted_psd_charge.chargeX;
			vars.PSD_chargeY = extracted_psd_charge.chargeY;
			vars.PSD_charge = 0.5 * (extracted_psd_charge.chargeX + extracted_psd_charge.chargeY);

			// Filling classifiers info
			if (bgoVault.GetLastFFracLayer() != -1)
			{
				vars.xtr = 0.125e-6 * pow(bgoVault.GetSumRMS(), 4) * bgoVault.GetSingleFracLayer(13);
				vars.xtrl = 0.125e-6 * pow(bgoVault.GetSumRMS(), 4) * bgoVault.GetLastFFracLayer();
			}
			else
			{
				vars.xtr = -1;
				vars.xtrl = -1;
			}

			// Filling attitude info
			fill_attitude_info(attitude, vars);

			// Fill tree
			if (bgoTotalE_raw >=20 && bgoTotalE_raw * _GeV <=100)
				DmpNtupTree_20_100->Fill();
			if (bgoTotalE_raw >100 && bgoTotalE_raw * _GeV <=250)
				DmpNtupTree_100_250->Fill();
			if (bgoTotalE_raw >250 && bgoTotalE_raw * _GeV <=500)
				DmpNtupTree_250_500->Fill();
			if (bgoTotalE_raw >500 && bgoTotalE_raw * _GeV <=1000)
				DmpNtupTree_500_1000->Fill();
			if (bgoTotalE_raw >1000 && bgoTotalE_raw * _GeV <=5000)
				DmpNtupTree_1000_5000->Fill();
			if (bgoTotalE_raw >5000 && bgoTotalE_raw * _GeV <=10000)
				DmpNtupTree_5000_10000->Fill();
		}
	}
	
	auto outFilePath_20_100 = uniqueTupleOutFile(opt, data_year, data_month, 20, 100);
	auto outFilePath_100_250 = uniqueTupleOutFile(opt, data_year, data_month, 100, 250);
	auto outFilePath_250_500 = uniqueTupleOutFile(opt, data_year, data_month, 250, 500);
	auto outFilePath_500_1000 = uniqueTupleOutFile(opt, data_year, data_month, 500, 1000);
	auto outFilePath_1000_5000 = uniqueTupleOutFile(opt, data_year, data_month, 1000, 5000);
	auto outFilePath_5000_10000 = uniqueTupleOutFile(opt, data_year, data_month, 5000, 10000);
	
    TFile outFile_20_100(outFilePath_20_100.c_str(), "NEW", "Analysis Output File");
    if (!outFile_20_100.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_20_100 << std::endl;
        exit(123);
    }
	DmpNtupTree_20_100->Write();
	outFile_20_100.Close();

	TFile outFile_100_250(outFilePath_100_250.c_str(), "NEW", "Analysis Output File");
    if (!outFile_100_250.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_100_250 << std::endl;
        exit(123);
    }
	DmpNtupTree_100_250->Write();
	outFile_100_250.Close();

	TFile outFile_250_500(outFilePath_250_500.c_str(), "NEW", "Analysis Output File");
    if (!outFile_250_500.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_250_500 << std::endl;
        exit(123);
    }
	DmpNtupTree_250_500->Write();
	outFile_250_500.Close();

	TFile outFile_500_1000(outFilePath_500_1000.c_str(), "NEW", "Analysis Output File");
    if (!outFile_500_1000.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_500_1000 << std::endl;
        exit(123);
    }
	DmpNtupTree_500_1000->Write();
	outFile_500_1000.Close();

	TFile outFile_1000_5000(outFilePath_1000_5000.c_str(), "NEW", "Analysis Output File");
    if (!outFile_1000_5000.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_1000_5000 << std::endl;
        exit(123);
    }
	DmpNtupTree_1000_5000->Write();
	outFile_1000_5000.Close();

	TFile outFile_5000_10000(outFilePath_5000_10000.c_str(), "NEW", "Analysis Output File");
    if (!outFile_5000_10000.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_5000_10000 << std::endl;
        exit(123);
    }
	DmpNtupTree_5000_10000->Write();
	outFile_5000_10000.Close();
	
}

template void linker(
	std::shared_ptr<TTree> tree,
	t_variables &vars);

template void linker(
	std::shared_ptr<TChain> chain,
	t_variables &vars);