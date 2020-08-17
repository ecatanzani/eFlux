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

#include "TEfficiency.h"
#include "TTree.h"

#include "DmpEvtAttitude.h"
#include "DmpFilterOrbit.h"
#include "DmpEvtHeader.h"
#include "DmpIOSvc.h"
#include "DmpCore.h"

struct t_variables
{
	// Trigger
	bool unbiased_trigger;
	bool mip1_trigger;
	bool mip2_trigger;
	bool HET_trigger;
	bool LET_trigger;
	bool MIP_trigger;

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
	double STK_charge;

	// BGO

	double energy;
	double energy_corr;
	double BGOrec_slopeX;
	double BGOrec_slopeY;
	double BGOrec_interceptX;
	double BGOrec_interceptY;
	double sumRms;
	double fracLast;
	unsigned int lastBGOLayer;
	unsigned int nBGOentries;

	// PSD
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
};

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

inline void branchTree(TTree &tree, t_variables &vars)
{
	tree.Branch("unbiased_trigger", &vars.unbiased_trigger, "unbiased_trigger/O");
	tree.Branch("mip1_trigger", &vars.mip1_trigger, "mip1_trigger/O");
	tree.Branch("mip2_trigger", &vars.mip2_trigger, "mip2_trigger/O");
	tree.Branch("HET_trigger", &vars.HET_trigger, "HET_trigger/O");
	tree.Branch("LET_trigger", &vars.LET_trigger, "LET_trigger/O");
	tree.Branch("MIP_trigger", &vars.MIP_trigger, "MIP_trigger/O");
	tree.Branch("second", &vars.second, "second/i");
	tree.Branch("msecond", &vars.msecond, "msecond/i");
	tree.Branch("STK_bestTrack_npoints", &vars.STK_bestTrack_npoints, "STK_bestTrack_npoints/i");
	tree.Branch("STK_bestTrack_nholesX", &vars.STK_bestTrack_nholesX, "STK_bestTrack_nholesX/i");
	tree.Branch("STK_bestTrack_nholesY", &vars.STK_bestTrack_nholesY, "STK_bestTrack_nholesY/i");
	tree.Branch("STK_bestTrack_slopeX", &vars.STK_bestTrack_slopeX, "STK_bestTrack_slopeX/D");
	tree.Branch("STK_bestTrack_slopeY", &vars.STK_bestTrack_slopeY, "STK_bestTrack_slopeY/D");
	tree.Branch("STK_bestTrack_interceptX", &vars.STK_bestTrack_interceptX, "STK_bestTrack_interceptX/D");
	tree.Branch("STK_bestTrack_interceptY", &vars.STK_bestTrack_interceptY, "STK_bestTrack_interceptY/D");
	tree.Branch("STK_bestTrack_costheta", &vars.STK_bestTrack_costheta, "STK_bestTrack_costheta/D");
	tree.Branch("STK_bestTrack_phi", &vars.STK_bestTrack_phi, "STK_bestTrack_phi/D");
	tree.Branch("STK_bestTrack_extr_BGO_topX", &vars.STK_bestTrack_extr_BGO_topX, "STK_bestTrack_extr_BGO_topX/D");
	tree.Branch("STK_bestTrack_extr_BGO_topY", &vars.STK_bestTrack_extr_BGO_topY, "STK_bestTrack_extr_BGO_topY/D");
	tree.Branch("STK_bestTrack_STK_BGO_topX_distance", &vars.STK_bestTrack_STK_BGO_topX_distance, "STK_bestTrack_STK_BGO_topX_distance/D");
	tree.Branch("STK_bestTrack_STK_BGO_topY_distance", &vars.STK_bestTrack_STK_BGO_topY_distance, "STK_bestTrack_STK_BGO_topY_distance/D");
	tree.Branch("STK_bestTrack_angular_distance_STK_BGO", &vars.STK_bestTrack_angular_distance_STK_BGO, "STK_bestTrack_angular_distance_STK_BGO/D");
	tree.Branch("STK_charge", &vars.STK_charge, "STK_charge/D");
	tree.Branch("energy", &vars.energy, "energy/D");
	tree.Branch("energy_corr", &vars.energy_corr, "energy_corr/D");
	tree.Branch("BGOrec_slopeX", &vars.BGOrec_slopeX, "BGOrec_slopeX/D");
	tree.Branch("BGOrec_slopeY", &vars.BGOrec_slopeY, "BGOrec_slopeY/D");
	tree.Branch("BGOrec_interceptX", &vars.BGOrec_interceptX, "BGOrec_interceptX/D");
	tree.Branch("BGOrec_interceptY", &vars.BGOrec_interceptY, "BGOrec_interceptY/D");
	tree.Branch("sumRms", &vars.sumRms, "sumRms/D");
	tree.Branch("fracLast", &vars.fracLast, "fracLast/D");
	tree.Branch("lastBGOLayer", &vars.lastBGOLayer, "lastBGOLayer/i");
	tree.Branch("nBGOentries", &vars.nBGOentries, "nBGOentries/i");
	tree.Branch("PSD_charge", &vars.PSD_charge, "PSD_charge/D");
	tree.Branch("xtr", &vars.xtr, "xtr/D");
	tree.Branch("xtrl", &vars.xtrl, "xtrl/D");
	tree.Branch("glat", &vars.glat, "glat/D");
	tree.Branch("glon", &vars.glon, "glon/D");
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
	TTree DmpNtupTree_20_100("DmpEvtNtup", "DAMPE Event nTuple Tree");
	TTree DmpNtupTree_100_250("DmpEvtNtup", "DAMPE Event nTuple Tree");
	TTree DmpNtupTree_250_500("DmpEvtNtup", "DAMPE Event nTuple Tree");
	TTree DmpNtupTree_500_1000("DmpEvtNtup", "DAMPE Event nTuple Tree");
	TTree DmpNtupTree_1000_5000("DmpEvtNtup", "DAMPE Event nTuple Tree");
	TTree DmpNtupTree_5000_10000("DmpEvtNtup", "DAMPE Event nTuple Tree");

	t_variables vars;

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

		if (pFilter->IsInSAA(evt_header->GetSecond()))
		{
			continue;
			++data_selection.events_in_saa;
		}

		// Event printout
		if (verbose)
			updateProcessStatus(evIdx, kStep, nevents);

		// Get event total energy
		double bgoTotalE = bgorec->GetTotalEnergy();
		double bgoTotalE_corr = bgorec->GetElectronEcor();

		// Don't accept events outside the selected energy window
		if (bgoTotalE * _GeV < flux_cuts.min_event_energy || bgoTotalE * _GeV > flux_cuts.max_event_energy)
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
			bgoTotalE,
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
				bgoTotalE,
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

		if (filter.all_cut_no_xtrl)
		{
			// Filling STK best track info
			fill_STK_bestTrack_info(event_best_track, vars);

			// Filling BGO info
			vars.energy = bgoTotalE;
			vars.energy_corr = bgoTotalE_corr;
			vars.BGOrec_slopeX = bgorec->GetSlopeXZ();
			vars.BGOrec_slopeY = bgorec->GetSlopeYZ();
			vars.BGOrec_interceptX = bgorec->GetInterceptXZ();
			vars.BGOrec_interceptY = bgorec->GetInterceptYZ();
			vars.sumRms = bgoVault.GetSumRMS();
			vars.fracLast = bgoVault.GetLastFFracLayer();
			vars.lastBGOLayer = bgoVault.GetFracIdxLastLayer();
			vars.nBGOentries = bgoVault.GetNhits();

			// Filling charges info
			vars.STK_charge = 0.5 * (extracted_stk_charge.chargeX + extracted_stk_charge.chargeY);
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
			if (bgoTotalE >=20 && bgoTotalE * _GeV <=100)
				DmpNtupTree_20_100.Fill();
			if (bgoTotalE >100 && bgoTotalE * _GeV <=250)
				DmpNtupTree_100_250.Fill();
			if (bgoTotalE >250 && bgoTotalE * _GeV <=500)
				DmpNtupTree_250_500.Fill();
			if (bgoTotalE >500 && bgoTotalE * _GeV <=1000)
				DmpNtupTree_500_1000.Fill();
			if (bgoTotalE >1000 && bgoTotalE * _GeV <=5000)
				DmpNtupTree_1000_5000.Fill();
			if (bgoTotalE >5000 && bgoTotalE * _GeV <=10000)
				DmpNtupTree_5000_10000.Fill();
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
	DmpNtupTree_20_100.Write();
	outFile_20_100.Close();

	TFile outFile_100_250(outFilePath_100_250.c_str(), "NEW", "Analysis Output File");
    if (!outFile_100_250.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_100_250 << std::endl;
        exit(123);
    }
	DmpNtupTree_100_250.Write();
	outFile_100_250.Close();

	TFile outFile_250_500(outFilePath_250_500.c_str(), "NEW", "Analysis Output File");
    if (!outFile_250_500.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_250_500 << std::endl;
        exit(123);
    }
	DmpNtupTree_250_500.Write();
	outFile_250_500.Close();

	TFile outFile_500_1000(outFilePath_500_1000.c_str(), "NEW", "Analysis Output File");
    if (!outFile_500_1000.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_500_1000 << std::endl;
        exit(123);
    }
	DmpNtupTree_500_1000.Write();
	outFile_500_1000.Close();

	TFile outFile_1000_5000(outFilePath_1000_5000.c_str(), "NEW", "Analysis Output File");
    if (!outFile_1000_5000.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_1000_5000 << std::endl;
        exit(123);
    }
	DmpNtupTree_1000_5000.Write();
	outFile_1000_5000.Close();

	TFile outFile_5000_10000(outFilePath_5000_10000.c_str(), "NEW", "Analysis Output File");
    if (!outFile_5000_10000.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath_5000_10000 << std::endl;
        exit(123);
    }
	DmpNtupTree_5000_10000.Write();
	outFile_5000_10000.Close();
	
}