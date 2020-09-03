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

#include "TEfficiency.h"

#include "DmpFilterOrbit.h"
#include "DmpEvtHeader.h"
#include "DmpIOSvc.h"
#include "DmpCore.h"

inline void init_BGO_histos(std::vector<TH1D> &h_layer_energy_ratio)
{
	h_layer_energy_ratio.resize(DAMPE_bgo_nLayers);

	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
	{
		TString h_ratio_name = "h_layer_energy_ratio_";
		TString h_ratio_title = "Energy Ratio - BGO layer ";

		h_ratio_name += lIdx;
		h_ratio_title += lIdx;

		h_layer_energy_ratio[lIdx] = TH1D(h_ratio_name.Data(), h_ratio_title.Data(), 100, 0, 1);
		h_layer_energy_ratio[lIdx].Sumw2();
	}
}

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

std::vector<TH1D> evLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const std::string wd)
{
	auto dmpch = aggregateDataEventsTChain(inputPath, verbose);

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
	
	data_evt_time evt_time;
	dmpch->GetEntry(0);
	evt_time.start_second = evt_header->GetSecond();
	evt_time.start_msecond = evt_header->GetMillisecond();
	dmpch->GetEntry(nevents -1);
	evt_time.end_second = evt_header->GetSecond();
	evt_time.end_msecond = evt_header->GetMillisecond();
	
	if (verbose)
	{
		std::cout << "RUN (info) \tTotal number of events: " << nevents << std::endl;
		std::cout << "RUN (info) \tStart Time (second): " << evt_time.start_second << std::endl;
		std::cout << "RUN (info) \tStart Time (millisecond): " << evt_time.start_msecond << std::endl;
		std::cout << "RUN (info) \tEnd Time (second): " << evt_time.end_second << std::endl;
		std::cout << "RUN (info) \tEnd Time (millisecond): " << evt_time.end_msecond << "\n\n";
	}

	// Cut histos
	TH1D h_trigger("h_trigger", "Energy Distribution of the triggered particles; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_cut("h_geometric_cut", "Energy Distribution - geometric (trigger selection) cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_maxElayer_cut("h_maxElayer_cut", "Energy Distribution - maxElayer cut; Raw Energy (GeV); counts ", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_maxBarLayer_cut("h_maxBarLayer_cut", "Energy Distribution - maxBarLayer cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOTrackContainment_cut("h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGO_fiducial_cut("h_BGO_fiducial_cut", "Energy Distibution - BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_all_cut("h_all_cut", "Electron counts energy distribution; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_all_cut_ce("h_all_cut_ce", "Electron counts energy distribution; Corrected Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// Cuts && Geometric Cut
	TH1D h_geometric_maxElayer_cut("h_geometric_maxElayer_cut", "Energy Distribution - maxElayer + geometric (trigger selection) cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_maxBarLayer_cut("h_geometric_maxBarLayer_cut", "Energy Distribution - maxBarLayer + geometric (trigger selection) cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_BGOTrackContainment_cut("h_geometric_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment + geometric (trigger selection) cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_BGO_fiducial_cut("h_geometric_BGO_fiducial_cut", "Energy Distibution - BGO fiducial + geometric (trigger selection) cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_all_cut("h_geometric_all_cut", "Energy Distribution - All + geometric (trigger selection) cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_all_cut_ce("h_geometric_all_cut_ce", "Energy Distribution - All + geometric (trigger selection) cut; Corrected Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// Cuts && BGO fiducial volume cut
	TH1D h_BGOfiducial_nBarLayer13_cut("h_BGOfiducial_nBarLayer13_cut", "Energy Distribution - nBarLayer13 + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_maxRms_cut("h_BGOfiducial_maxRms_cut", "Energy Distribution - maxRms  + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_track_selection_cut("h_BGOfiducial_track_selection_cut", "Energy Distribution - track selection + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_xtrl_cut("h_BGOfiducial_xtrl_cut", "Energy Distribution - xtrl + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_psd_stk_match_cut("h_BGOfiducial_psd_stk_match_cut", "Energy Distribution - psd-stk match + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_psd_charge_cut("h_BGOfiducial_psd_charge_cut", "Energy Distribution - psd charge + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_stk_charge_cut("h_BGOfiducial_stk_charge_cut", "Energy Distribution - stk charge + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_all_cut("h_BGOfiducial_all_cut", "Energy Distribution - All + BGO fiducial cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_all_cut_ce("h_BGOfiducial_all_cut_ce", "Energy Distribution - All + BGO fiducial cut; Corrected Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// Analysis histos - reco energy of incoming events
	TH1D h_BGOrec_energy("h_BGOrec_energy", "BGO Energy: Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// After Geometric Cut
	// Slope X and Y
	TH1D h_geo_BGOrec_slopeX("h_geo_BGOrec_slopeX", "BGOrec Slope X", 1000, -90, 90);
	TH1D h_geo_BGOrec_slopeY("h_geo_BGOrec_slopeY", "BGOrec Slope Y", 1000, -90, 90);

	// Intercept X and Y
	TH1D h_geo_BGOrec_interceptX("h_geo_BGOrec_interceptX", "BGOrec Intercept X; X(mm); counts", 500, -500, 500);
	TH1D h_geo_BGOrec_interceptY("h_geo_BGOrec_interceptY", "BGOrec Intercept Y; Y(mm); counts", 500, -500, 500);

	// Top Maps
	TH2D h_geo_BGOreco_topMap("h_geo_BGOreco_topMap", "BGOreco TOP Map; X(mm); Y(mm)", 500, -500, 500, 500, -500, 500);

	// Bottom Maps
	TH2D h_geo_BGOreco_bottomMap("h_geo_BGOreco_bottomMap", "BGOreco BOTTOM Map; X(mm); Y(mm)", 500, -500, 500, 500, -500, 500);

	// sumRms - cosine correlation
	std::vector<std::shared_ptr<TH2D>> sumRms_cosine (logEBins.size() - 1);
	auto cosine_bins = createLinearBinning(0, 1, 1e+2);
	auto sumRms_bins = createLogBinning(10, 2e+3, 1e+3);
	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
	{
		std::string histo_name = "sumRms_cosine_" + to_string(std::distance(sumRms_cosine.begin(), it));
		(*it) = std::make_shared<TH2D>(histo_name.c_str(), "sumRms - cos(#theta) correlation; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));
	}

	TH2D sumRms_cosine_20_100("sumRms_cosine_20_100", "sumRms - cos(#theta) correlation 20 GeV - 100 GeV; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));
	TH2D sumRms_cosine_100_250("sumRms_cosine_100_250", "sumRms - cos(#theta) correlation 100 GeV - 250 GeV; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));
	TH2D sumRms_cosine_250_500("sumRms_cosine_250_500", "sumRms - cos(#theta) correlation 250 GeV - 500 GeV; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));
	TH2D sumRms_cosine_500_1000("sumRms_cosine_500_1000", "sumRms - cos(#theta) correlation 500 GeV - 1 TeV; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));
	TH2D sumRms_cosine_1000_3000("sumRms_cosine_1000_3000", "sumRms - cos(#theta) correlation 1 TeV - 3 TeV; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));
	TH2D sumRms_cosine_3000_10000("sumRms_cosine_3000_10000", "sumRms - cos(#theta) correlation 3 TeV - 10 TeV; cos(#theta); sumRms [mm]", cosine_bins.size() -1, &(cosine_bins[0]),  sumRms_bins.size() -1, &(sumRms_bins[0]));


	// Ratio of layer energy respect to total BGO energy
	TH1D h_layer_max_energy_ratio("h_layer_max_energy_ratio", "Layer Energy Ratio", 100, 0, 1);
	std::vector<TH1D> h_layer_energy_ratio;
	init_BGO_histos(h_layer_energy_ratio);

	// XTRL histos
	auto xtrl_bins = createLinearBinning(0, 150, 1e+3);
	auto flast_binning = createLogBinning(1e-5, 2e-1, 1e+3);

	// Energy integrated XTRL distribution
	TH1D h_xtrl_energy_int("h_xtrl_energy_int", "Energy integrated XTRL distribution", xtrl_bins.size() - 1, &(xtrl_bins[0]));
	// 2D XTRL energy distribution
	TH2D h_xtrl("h_xtrl", "XTRL energy Distribution", logEBins.size() - 1, &(logEBins[0]), xtrl_bins.size() - 1, &(xtrl_bins[0]));
	TH2D e_discrimination_last("e_discrimination_last", "Electron Discrimination; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_last_20_100("e_discrimination_last_20_100", "Electron Discrimination 20 GeV - 100 GeV; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_last_100_250("e_discrimination_last_100_250", "Electron Discrimination 100 GeV - 250 GeV; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_last_250_500("e_discrimination_last_250_500", "Electron Discrimination 250 GeV - 500 GeV; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_last_500_1000("e_discrimination_last_500_1000", "Electron Discrimination 500 GeV - 1 TeV; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_last_1000_3000("e_discrimination_last_1000_3000", "Electron Discrimination 1 TeV - 3 TeV; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_last_3000_10000("e_discrimination_last_3000_10000", "Electron Discrimination 3 TeV - 10 TeV; sumRms [mm]; F_{last}", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));

	TH2D e_discrimination("e_discrimination", "Electron Discrimination; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_20_100("e_discrimination_20_100", "Electron Discrimination 20 GeV - 100 GeV; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_100_250("e_discrimination_100_250", "Electron Discrimination 100 GeV - 250 GeV; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_250_500("e_discrimination_250_500", "Electron Discrimination 250 GeV - 500 GeV; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_500_1000("e_discrimination_500_1000", "Electron Discrimination 500 GeV - 1 TeV; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_1000_3000("e_discrimination_1000_3000", "Electron Discrimination 1 TeV - 3 TeV; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));
	TH2D e_discrimination_3000_10000("e_discrimination_3000_10000", "Electron Discrimination 3 TeV - 10 TeV; sumRms [mm]; F", sumRms_bins.size() -1, &(sumRms_bins[0]), flast_binning.size() -1, &(flast_binning[0]));

	// Bin energy integrated XTRL distributions
	std::vector<std::shared_ptr<TH1D>> bin_xtrl(logEBins.size() - 1);
	for (auto it = bin_xtrl.begin(); it != bin_xtrl.end(); ++it)
	{
		std::string histo_name = "h_xtrl_bin_" + to_string(std::distance(bin_xtrl.begin(), it));
		(*it) = std::make_shared<TH1D>(histo_name.c_str(), "XTRL bin distribution; xtrl; counts", 100, 0, 150);
	}

	// PSD charge histos
	TH1D h_psd_chargeX("h_psd_chargeX", "Charge distribution X; X charge", 1000, 0, 100);
	TH1D h_psd_chargeY("h_psd_chargeY", "Charge distribution Y; Y charge", 1000, 0, 100);
	TH2D h_psd_charge2D("h_psd_charge2D", "PSD charge; X charge; Y charge", 1000, 0, 100, 1000, 0, 100);
	TH1D h_psd_charge("h_psd_charge", "Mean PSD charge; Mean charge", 1000, 0, 100);

	TH1D h_psd_selected_chargeX("h_psd_selected_chargeX", "Charge distribution X; X charge", 1000, 0, 100);
	TH1D h_psd_selected_chargeY("h_psd_selected_chargeY", "Charge distribution Y; Y charge", 1000, 0, 100);
	TH2D h_psd_selected_charge2D("h_psd_selected_charge2D", "PSD charge; X charge; Y charge", 1000, 0, 100, 1000, 0, 100);
	TH1D h_psd_selected_charge("h_psd_selected_charge", "Mean PSD charge; Mean charge", 1000, 0, 1000);

	// STK charge histos
	TH1D h_stk_chargeX("h_stk_chargeX", "Charge distribution X; X charge", 1000, 0, 100);
	TH1D h_stk_chargeY("h_stk_chargeY", "Charge distribution Y; Y charge", 1000, 0, 100);
	TH2D h_stk_charge2D("h_stk_charge2D", "STK charge; X charge; Y charge", 1000, 0, 100, 1000, 0, 100);
	TH1D h_stk_charge("h_stk_charge", "Mean STK charge; Mean charge", 1000, 0, 100);

	TH1D h_stk_selected_chargeX("h_stk_selected_chargeX", "Charge distribution X; X charge", 1000, 0, 100);
	TH1D h_stk_selected_chargeY("h_stk_selected_chargeY", "Charge distribution Y; Y charge", 1000, 0, 100);
	TH2D h_stk_selected_charge2D("h_stk_selected_charge2D", "STK charge; X charge; Y charge", 1000, 0, 100, 1000, 0, 100);
	TH1D h_stk_selected_charge("h_stk_selected_charge", "Mean STK charge; Mean charge", 1000, 0, 100);

	// Proton background histos
	TH1D h_background_under_xtrl_cut("h_background_under_xtrl_cut", "Proton background - XTRL < cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_background_over_xtrl_cut("h_background_over_xtrl_cut", "Proton background - XTRL > 20 and XTRL < 100 cut; Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	TH1D h_second("h_second", "second", 100, evt_time.start_second, evt_time.end_second);
	TH1D h_msecond("h_msecond", "millisecond", 100, evt_time.start_msecond, evt_time.end_msecond);

	// Sumw2 - First-Cut histos
	h_trigger.Sumw2();
	h_geometric_cut.Sumw2();
	h_maxElayer_cut.Sumw2();
	h_maxBarLayer_cut.Sumw2();
	h_BGOTrackContainment_cut.Sumw2();
	h_BGO_fiducial_cut.Sumw2();
	h_all_cut.Sumw2();
	h_all_cut_ce.Sumw2();

	// Sumw2 - Cuts && Geometric Cut
	h_geometric_maxElayer_cut.Sumw2();
	h_geometric_maxBarLayer_cut.Sumw2();
	h_geometric_BGOTrackContainment_cut.Sumw2();
	h_geometric_BGO_fiducial_cut.Sumw2();
	h_geometric_all_cut.Sumw2();
	h_geometric_all_cut_ce.Sumw2();

	// Sumw2 - Cuts && BGO fiducial volume cut
	h_BGOfiducial_nBarLayer13_cut.Sumw2();
	h_BGOfiducial_maxRms_cut.Sumw2();
	h_BGOfiducial_track_selection_cut.Sumw2();
	h_BGOfiducial_psd_stk_match_cut.Sumw2();
	h_BGOfiducial_psd_charge_cut.Sumw2();
	h_BGOfiducial_stk_charge_cut.Sumw2();
	h_BGOfiducial_xtrl_cut.Sumw2();
	h_BGOfiducial_all_cut.Sumw2();
	h_BGOfiducial_all_cut_ce.Sumw2();

	// Sumw2 Analysis histos - reco energy of incoming events
	h_BGOrec_energy.Sumw2();

	// Sumw2 Analysis histos - Geo
	h_geo_BGOrec_slopeX.Sumw2();
	h_geo_BGOrec_slopeY.Sumw2();
	h_geo_BGOrec_interceptX.Sumw2();
	h_geo_BGOrec_interceptY.Sumw2();
	h_geo_BGOreco_topMap.Sumw2();
	h_geo_BGOreco_bottomMap.Sumw2();
	h_layer_max_energy_ratio.Sumw2();

	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
		(*it)->Sumw2();
	
	sumRms_cosine_20_100.Sumw2();
	sumRms_cosine_100_250.Sumw2();
	sumRms_cosine_250_500.Sumw2();
	sumRms_cosine_500_1000.Sumw2();
	sumRms_cosine_1000_3000.Sumw2();
	sumRms_cosine_3000_10000.Sumw2();

	// Sumw2 XTRL histos
	h_xtrl_energy_int.Sumw2();
	h_xtrl.Sumw2();
	for (auto it = bin_xtrl.begin(); it != bin_xtrl.end(); ++it)
		(*it)->Sumw2();
	e_discrimination.Sumw2();
	e_discrimination_20_100.Sumw2();
	e_discrimination_100_250.Sumw2();
	e_discrimination_250_500.Sumw2();
	e_discrimination_500_1000.Sumw2();
	e_discrimination_1000_3000.Sumw2();
	e_discrimination_3000_10000.Sumw2();

	e_discrimination_last.Sumw2();
	e_discrimination_last_20_100.Sumw2();
	e_discrimination_last_100_250.Sumw2();
	e_discrimination_last_250_500.Sumw2();
	e_discrimination_last_500_1000.Sumw2();
	e_discrimination_last_1000_3000.Sumw2();
	e_discrimination_last_3000_10000.Sumw2();

	// Sumw2 PSD charge histos
	h_psd_chargeX.Sumw2();
	h_psd_chargeY.Sumw2();
	h_psd_charge.Sumw2();
	h_psd_charge2D.Sumw2();

	h_psd_selected_chargeX.Sumw2();
	h_psd_selected_chargeY.Sumw2();
	h_psd_selected_charge.Sumw2();
	h_psd_selected_charge2D.Sumw2();

	// Sumw2 STK charge histos
	h_stk_chargeX.Sumw2();
	h_stk_chargeY.Sumw2();
	h_stk_charge.Sumw2();
	h_stk_charge2D.Sumw2();

	h_stk_selected_chargeX.Sumw2();
	h_stk_selected_chargeY.Sumw2();
	h_stk_selected_charge.Sumw2();
	h_stk_selected_charge2D.Sumw2();

	// Proton background histos
	h_background_under_xtrl_cut.Sumw2();
	h_background_over_xtrl_cut.Sumw2();
	
	double _GeV = 0.001;
	int kStep = 10;

	if (verbose)
		std::cout << "Starting analysing data...\n\n";

	for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
	{
		// Get chain event
		dmpch->GetEvent(evIdx);

		// Update event counter
		++data_selection.event_counter;
		

		// Extract data time information
		int second = evt_header->GetSecond();
		short msecond = evt_header->GetMillisecond();
		if (pFilter->IsInSAA(second))
		{
			continue;
			++data_selection.events_in_saa;
		}
		h_second.Fill(second);
		h_msecond.Fill(msecond);

		// Event printout
		if (verbose)
			updateProcessStatus(evIdx, kStep, nevents);

		// Get event total energy
		double bgoTotalE = bgorec->GetTotalEnergy();

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

		// Check if the event has been triggered or not
		if (general_trigger)
		{
			++data_selection.triggered_events;
			h_trigger.Fill(bgoTotalE * _GeV);
			if (!checkBGOreco_data(bgorec))
				continue;
		}
		else
			continue;

		// Fill the energy histos only for good reco events
		h_BGOrec_energy.Fill(bgoTotalE * _GeV);

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

		// evaluate the energy raio on each single layer of the BGO
		evaluateEnergyRatio(
			bgorec,
			flux_cuts,
			bgoTotalE,
			h_layer_max_energy_ratio,
			h_layer_energy_ratio);

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
		
		// Fill cuts histos

		// Fill geometric cut histos
		if (filter.geometric)
		{
			h_geometric_cut.Fill(bgoTotalE * _GeV);

			// Evaluate the position on the First BGO layer (after geometric cut)
			evaluateTopBottomPosition_data(
				bgorec,
				h_geo_BGOrec_slopeX,
				h_geo_BGOrec_slopeY,
				h_geo_BGOrec_interceptX,
				h_geo_BGOrec_interceptY,
				h_geo_BGOreco_topMap,
				h_geo_BGOreco_bottomMap);

			// Geometric cut && maxElayer cut
			if (filter.BGO_fiducial_maxElayer_cut)
				h_geometric_maxElayer_cut.Fill(bgoTotalE * _GeV);

			// Geometric cut && maxBarLayer cut
			if (filter.BGO_fiducial_maxBarLayer_cut)
				h_geometric_maxBarLayer_cut.Fill(bgoTotalE * _GeV);

			// Geometric cut && BGOTrackContainment cut
			if (filter.BGO_fiducial_BGOTrackContainment_cut)
				h_geometric_BGOTrackContainment_cut.Fill(bgoTotalE * _GeV);

			// Geometric cut && BGO fiducial cut
			if (filter.BGO_fiducial)
				h_geometric_BGO_fiducial_cut.Fill(bgoTotalE * _GeV);

			// Geometric cut and all cuts
			if (filter.all_cut)
			{
				h_geometric_all_cut.Fill(bgoTotalE * _GeV);
				h_geometric_all_cut_ce.Fill(bgorec->GetElectronEcor() * _GeV);
			}
		}

		// Fill BGO_fiducial_maxElayer cut histo
		if (filter.BGO_fiducial_maxElayer_cut)
			h_maxElayer_cut.Fill(bgoTotalE * _GeV);

		// Fill BGO_fiducial_maxBarLayer cut histo
		if (filter.BGO_fiducial_maxBarLayer_cut)
			h_maxBarLayer_cut.Fill(bgoTotalE * _GeV);

		// Fill BGO_fiducial_BGOTrackContainment cut histo
		if (filter.BGO_fiducial_BGOTrackContainment_cut)
			h_BGOTrackContainment_cut.Fill(bgoTotalE * _GeV);

		// Fill BGO fiducial volume cut
		if (filter.BGO_fiducial)
		{
			h_BGO_fiducial_cut.Fill(bgoTotalE * _GeV);

			// BGO fiducial cut && nBarLayer13 cut
			if (filter.nBarLayer13_cut)
				h_BGOfiducial_nBarLayer13_cut.Fill(bgoTotalE * _GeV);

			// BGO fiducial cut && maxRms cut
			if (filter.maxRms_cut)
				h_BGOfiducial_maxRms_cut.Fill(bgoTotalE * _GeV);

			// BGO fiducial cut && track selection cut
			if (filter.track_selection_cut)
				h_BGOfiducial_track_selection_cut.Fill(bgoTotalE * _GeV);
			
			// BGO fiducial cut && PSD-STK match cut
			if (filter.psd_stk_match_cut)
				h_BGOfiducial_psd_stk_match_cut.Fill(bgoTotalE * _GeV);
			
			// BGO fiducial cut && PSD charge cut
			if (filter.psd_charge_cut)
				h_BGOfiducial_psd_charge_cut.Fill(bgoTotalE * _GeV);

			// BGO fiducial cut && STK charge cut
			if (filter.stk_charge_cut)
				h_BGOfiducial_stk_charge_cut.Fill(bgoTotalE * _GeV);

			// BGO fiducial cut && XTRL cut
			if (filter.xtrl_cut)
				h_BGOfiducial_xtrl_cut.Fill(bgoTotalE * _GeV);

			// BGO fiducial cut && all cut
			if (filter.all_cut)
			{
				h_BGOfiducial_all_cut.Fill(bgoTotalE * _GeV);
				h_BGOfiducial_all_cut_ce.Fill(bgorec->GetElectronEcor() * _GeV);
			}
		}

		// Collect PSD particle charge
		if (filter.psd_charge_measurement)
			fillChargeHistos<psd_charge>(
				h_psd_chargeX,
				h_psd_chargeY,
				h_psd_charge,
				h_psd_charge2D,
				extracted_psd_charge);
		
		// Collect selected STK particle charge
		if (filter.psd_charge_cut)
			fillChargeHistos<psd_charge>(
				h_psd_selected_chargeX,
				h_psd_selected_chargeY,
				h_psd_selected_charge,
				h_psd_selected_charge2D,
				extracted_psd_charge);

		// Collect STK particle charge
		if (filter.stk_charge_measurement)
			fillChargeHistos<stk_charge>(
				h_stk_chargeX,
				h_stk_chargeY,
				h_stk_charge,
				h_stk_charge2D,
				extracted_stk_charge);

		// Collect selected STK particle charge
		if (filter.stk_charge_cut)
			fillChargeHistos<stk_charge>(
				h_stk_selected_chargeX,
				h_stk_selected_chargeY,
				h_stk_selected_charge,
				h_stk_selected_charge2D,
				extracted_stk_charge);

		// Fill all cut histo
		if (active_cuts.nActiveCuts)
		{
			if (filter.all_cut)
			{
				h_all_cut.Fill(bgoTotalE * _GeV);
				h_all_cut_ce.Fill(bgorec->GetElectronEcor() * _GeV);

				// Fill sumRms - cosine correlation histo
				fill_sumRms_cosine_histo(
					bgoVault.GetSumRMS(),
					event_best_track.myBestTrack.getDirection().CosTheta(),
					bgorec->GetElectronEcor() * _GeV,
					logEBins,
					sumRms_cosine,
					sumRms_cosine_20_100,
					sumRms_cosine_100_250,
					sumRms_cosine_250_500,
					sumRms_cosine_500_1000,
					sumRms_cosine_1000_3000,
					sumRms_cosine_3000_10000);
			}

			if (filter.all_cut_no_xtrl)
			{
				// Fill xtrl histos
				fill_XTRL_histo(
					bgoVault.GetSumRMS(),
					bgoVault.GetLastFFracLayer(),
					bgorec->GetElectronEcor(),
					bin_xtrl,
					h_xtrl_energy_int,
					h_xtrl);

				// Fill e/p discrimination histos
				fill_ep_histos(
					bgoVault.GetSumRMS(),
					bgoVault.GetLastFFracLayer(),
					bgoVault.GetSingleFracLayer(13),
					e_discrimination,
					e_discrimination_20_100,
					e_discrimination_100_250,
					e_discrimination_250_500,
					e_discrimination_500_1000,
					e_discrimination_1000_3000,
					e_discrimination_3000_10000,
					e_discrimination_last,
					e_discrimination_last_20_100,
					e_discrimination_last_100_250,
					e_discrimination_last_250_500,
					e_discrimination_last_500_1000,
					e_discrimination_last_1000_3000,
					e_discrimination_last_3000_10000,
					//bgoTotalE * _GeV);
					bgorec->GetElectronEcor() * _GeV);
					
				// Compute proton background
				compute_proton_background(
					bgoVault.GetSumRMS(),
					bgoVault.GetLastFFracLayer(),
					bgorec->GetElectronEcor(),
					flux_cuts,
					h_background_under_xtrl_cut,
					h_background_over_xtrl_cut);
			}
		}
	}

	if (verbose)
	{
		std::cout << "\n\n ****** \n\n";
		std::cout << "Triggered events: " << data_selection.triggered_events << std::endl;
		std::cout << "Triggered events in energy range: " << data_selection.events_in_range << std::endl;
		std::cout << "Triggered events out of energy range: " << data_selection.events_out_range << std::endl;
		std::cout << "Events in SAA: " << data_selection.events_in_saa << std::endl;
		std::cout << "\n\n**** Filter result ****\n";
		std::cout << "***********************\n\n";
		std::cout << "Particles surviving the selection cuts: " << data_selection.selected_events << "\t | " << ((double)data_selection.selected_events / data_selection.triggered_events) * 100 << "%";
		std::cout << "\n\n ***********************\n\n";
	}
	
	outFile.cd();
	
	// Write histos to file
	// First-Cut histos
	h_trigger.Write();
	h_geometric_cut.Write();
	h_maxElayer_cut.Write();
	h_maxBarLayer_cut.Write();
	h_BGOTrackContainment_cut.Write();
	h_BGO_fiducial_cut.Write();
	h_all_cut.Write();
	h_all_cut_ce.Write();
	// Cuts && Geometric Cut
	h_geometric_maxElayer_cut.Write();
	h_geometric_maxBarLayer_cut.Write();
	h_geometric_BGOTrackContainment_cut.Write();
	h_geometric_BGO_fiducial_cut.Write();
	h_geometric_all_cut.Write();
	h_geometric_all_cut_ce.Write();
	// Cuts && BGO fiducial volume cut
	h_BGOfiducial_nBarLayer13_cut.Write();
	h_BGOfiducial_maxRms_cut.Write();
	h_BGOfiducial_track_selection_cut.Write();
	h_BGOfiducial_psd_stk_match_cut.Write();
	h_BGOfiducial_psd_charge_cut.Write();
	h_BGOfiducial_stk_charge_cut.Write();
	h_BGOfiducial_xtrl_cut.Write();
	h_BGOfiducial_all_cut.Write();
	h_BGOfiducial_all_cut_ce.Write();

	// Create output ratio dir in the output TFile
	auto ratioDir = outFile.mkdir("Efficiency");

	// Create trigger folder
	auto trigger_dir = ratioDir->mkdir("Trigger");
	trigger_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> tr_eff_gometric_cut;
	std::shared_ptr<TEfficiency> tr_eff_maxElayer_cut;
	std::shared_ptr<TEfficiency> tr_eff_maxBarLayer_cut;
	std::shared_ptr<TEfficiency> tr_eff_BGOTrackContainment_cut;
	std::shared_ptr<TEfficiency> tr_eff_BGO_fiducial_cut;
	std::shared_ptr<TEfficiency> tr_eff_all_cut;

	if (TEfficiency::CheckConsistency(h_geometric_cut, h_trigger))
		tr_eff_gometric_cut = std::make_shared<TEfficiency>(h_geometric_cut, h_trigger);

	if (TEfficiency::CheckConsistency(h_maxElayer_cut, h_trigger))
		tr_eff_maxElayer_cut = std::make_shared<TEfficiency>(h_maxElayer_cut, h_trigger);

	if (TEfficiency::CheckConsistency(h_maxBarLayer_cut, h_trigger))
		tr_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(h_maxBarLayer_cut, h_trigger);

	if (TEfficiency::CheckConsistency(h_BGOTrackContainment_cut, h_trigger))
		tr_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(h_BGOTrackContainment_cut, h_trigger);

	if (TEfficiency::CheckConsistency(h_BGO_fiducial_cut, h_trigger))
		tr_eff_BGO_fiducial_cut = std::make_shared<TEfficiency>(h_BGO_fiducial_cut, h_trigger);

	if (TEfficiency::CheckConsistency(h_all_cut, h_trigger))
		tr_eff_all_cut = std::make_shared<TEfficiency>(h_all_cut, h_trigger);

	// Set uniform statistic option
	tr_eff_gometric_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_BGO_fiducial_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);
	
	tr_eff_gometric_cut->SetName("tr_eff_gometric_cut");
	tr_eff_maxElayer_cut->SetName("tr_eff_maxElayer_cut");
	tr_eff_maxBarLayer_cut->SetName("tr_eff_maxBarLayer_cut");
	tr_eff_BGOTrackContainment_cut->SetName("tr_eff_BGOTrackContainment_cut");
	tr_eff_BGO_fiducial_cut->SetName("tr_eff_BGO_fiducial_cut");
	tr_eff_all_cut->SetName("tr_eff_all_cut");
	
	tr_eff_gometric_cut->SetTitle("Gometric cut efficiency");
	tr_eff_maxElayer_cut->SetTitle("maxElayer cut efficiency");
	tr_eff_maxBarLayer_cut->SetTitle("maxBarLayer cut efficiency");
	tr_eff_BGOTrackContainment_cut->SetTitle("BGOTrackContainment cut efficiency");
	tr_eff_all_cut->SetTitle("all cut efficiency");

	// Write histos to disk
	tr_eff_gometric_cut->Write();
	tr_eff_maxElayer_cut->Write();
	tr_eff_maxBarLayer_cut->Write();
	tr_eff_BGOTrackContainment_cut->Write();
	tr_eff_BGO_fiducial_cut->Write();
	tr_eff_all_cut->Write();

	// Create geometric folder
	auto geometric_dir = ratioDir->mkdir("Geometric");
	geometric_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> geo_eff_maxElayer_cut;
	std::shared_ptr<TEfficiency> geo_eff_maxBarLayer_cut;
	std::shared_ptr<TEfficiency> geo_eff_BGOTrackContainment_cut;
	std::shared_ptr<TEfficiency> geo_eff_BGO_fiducial;
	std::shared_ptr<TEfficiency> geo_eff_all_cut;

	if (TEfficiency::CheckConsistency(h_geometric_maxElayer_cut, h_geometric_cut))
		geo_eff_maxElayer_cut = std::make_shared<TEfficiency>(h_geometric_maxElayer_cut, h_geometric_cut);

	if (TEfficiency::CheckConsistency(h_geometric_maxBarLayer_cut, h_geometric_cut))
		geo_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(h_geometric_maxBarLayer_cut, h_geometric_cut);

	if (TEfficiency::CheckConsistency(h_geometric_BGOTrackContainment_cut, h_geometric_cut))
		geo_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(h_geometric_BGOTrackContainment_cut, h_geometric_cut);

	if (TEfficiency::CheckConsistency(h_geometric_BGO_fiducial_cut, h_geometric_cut))
		geo_eff_BGO_fiducial = std::make_shared<TEfficiency>(h_geometric_BGO_fiducial_cut, h_geometric_cut);

	if (TEfficiency::CheckConsistency(h_geometric_all_cut, h_geometric_cut))
		geo_eff_all_cut = std::make_shared<TEfficiency>(h_geometric_all_cut, h_geometric_cut);

	// Set uniform statistic option
	geo_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_BGO_fiducial->SetStatisticOption(TEfficiency::kBUniform);
	geo_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

	geo_eff_maxElayer_cut->SetName("geo_eff_maxElayer_cut");
	geo_eff_maxBarLayer_cut->SetName("geo_eff_maxBarLayer_cut");
	geo_eff_BGOTrackContainment_cut->SetName("geo_eff_BGOTrackContainment_cut");
	geo_eff_BGO_fiducial->SetName("geo_eff_BGO_fiducial");
	geo_eff_all_cut->SetName("geo_eff_all_cut");

	geo_eff_maxElayer_cut->SetTitle("geometic maxElayer cut efficiency");
	geo_eff_maxBarLayer_cut->SetTitle("geometric maxBarLayer cut efficiency");
	geo_eff_BGOTrackContainment_cut->SetTitle("geometric BGOTrackContainment cut efficiency");
	geo_eff_BGO_fiducial->SetTitle("geometric BGO fiducial cut efficiency");
	geo_eff_all_cut->SetTitle("geometric all cut efficiency");

	//Write histos to disk
	geo_eff_maxElayer_cut->Write();
	geo_eff_maxBarLayer_cut->Write();
	geo_eff_BGOTrackContainment_cut->Write();
	geo_eff_BGO_fiducial->Write();
	geo_eff_all_cut->Write();

	// Create BGO_fiducial_volume folder
	auto BGOfiducial_dir = ratioDir->mkdir("BGO_fiducial_volume");
	BGOfiducial_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> BGOfiducial_eff_nBarLayer13_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_maxRms_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_track_selection_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_psd_stk_match_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_psd_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_stk_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_xtrl_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_all_cut;

	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_maxRms_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_track_selection_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut;
	std::shared_ptr<TEfficiency> BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut;


	if (TEfficiency::CheckConsistency(h_BGOfiducial_nBarLayer13_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_nBarLayer13_cut = std::make_shared<TEfficiency>(h_BGOfiducial_nBarLayer13_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_maxRms_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_maxRms_cut = std::make_shared<TEfficiency>(h_BGOfiducial_maxRms_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_track_selection_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_track_selection_cut = std::make_shared<TEfficiency>(h_BGOfiducial_track_selection_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_stk_match_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_stk_match_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_charge_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_charge_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_stk_charge_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_stk_charge_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_xtrl_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut = std::make_shared<TEfficiency>(h_BGOfiducial_xtrl_cut, h_BGO_fiducial_cut);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_all_cut, h_BGO_fiducial_cut))
		BGOfiducial_eff_all_cut = std::make_shared<TEfficiency>(h_BGOfiducial_all_cut, h_BGO_fiducial_cut);

	if (active_cuts.maxRms && active_cuts.nBarLayer13)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_maxRms_cut, h_BGOfiducial_nBarLayer13_cut))
			BGOfiducial_eff_maxRms_cut = std::make_shared<TEfficiency>(h_BGOfiducial_maxRms_cut, h_BGOfiducial_nBarLayer13_cut);

	if (active_cuts.track_selection && active_cuts.maxRms)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_track_selection_cut, h_BGOfiducial_maxRms_cut))
			BGOfiducial_eff_track_selection_cut = std::make_shared<TEfficiency>(h_BGOfiducial_track_selection_cut, h_BGOfiducial_maxRms_cut);

	if (active_cuts.psd_stk_match && active_cuts.track_selection)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_stk_match_cut, h_BGOfiducial_track_selection_cut))
			BGOfiducial_eff_psd_stk_match_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_stk_match_cut, h_BGOfiducial_track_selection_cut);

	if (active_cuts.psd_charge && active_cuts.psd_stk_match)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_charge_cut, h_BGOfiducial_psd_stk_match_cut))
			BGOfiducial_eff_psd_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_charge_cut, h_BGOfiducial_psd_stk_match_cut);

	if (active_cuts.stk_charge && active_cuts.psd_charge)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_stk_charge_cut, h_BGOfiducial_psd_stk_match_cut))
			BGOfiducial_eff_stk_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_stk_charge_cut, h_BGOfiducial_psd_stk_match_cut);

	if (active_cuts.xtrl && active_cuts.stk_charge)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_xtrl_cut, h_BGOfiducial_stk_charge_cut))
			BGOfiducial_eff_xtrl_cut = std::make_shared<TEfficiency>(h_BGOfiducial_xtrl_cut, h_BGOfiducial_stk_charge_cut);

	// Set uniform statistic option
	BGOfiducial_eff_l13_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);

	BGOfiducial_eff_nBarLayer13_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_maxRms_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_track_selection_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_psd_stk_match_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_psd_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_stk_charge_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_xtrl_cut->SetStatisticOption(TEfficiency::kBUniform);
	BGOfiducial_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

	BGOfiducial_eff_l13_maxRms_cut->SetName("BGOfiducial_eff_l13_maxRms_cut");
	BGOfiducial_eff_l13_rms_track_selection_cut->SetName("BGOfiducial_eff_l13_rms_track_selection_cut");
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->SetName("BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut");
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->SetName("BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->SetName("BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->SetName("BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut");

	BGOfiducial_eff_nBarLayer13_cut->SetName("BGOfiducial_eff_nBarLayer13_cut");
	BGOfiducial_eff_maxRms_cut->SetName("BGOfiducial_eff_maxRms_cut");
	BGOfiducial_eff_track_selection_cut->SetName("BGOfiducial_eff_track_selection_cut");
	BGOfiducial_eff_psd_stk_match_cut->SetName("BGOfiducial_eff_psd_stk_match_cut");
	BGOfiducial_eff_psd_charge_cut->SetName("BGOfiducial_eff_psd_charge_cut");
	BGOfiducial_eff_stk_charge_cut->SetName("BGOfiducial_eff_stk_charge_cut");
	BGOfiducial_eff_xtrl_cut->SetName("BGOfiducial_eff_xtrl_cut");
	BGOfiducial_eff_all_cut->SetName("BGOfiducial_eff_all_cut");

	BGOfiducial_eff_l13_maxRms_cut->SetTitle("BGOfiducial l13 + maxRms cut efficiency");
	BGOfiducial_eff_l13_rms_track_selection_cut->SetTitle("BGOfiducial l13 + rms + track_selection cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->SetTitle("BGOfiducial l13 + rms + ts + PSD-STK match cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->SetTitle("BGOfiducial l13 + rms + ts + psdstk + PSD charge cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->SetTitle("BGOfiducial l13 + rms + ts + psdstk + pc + STK charge cut efficiency");
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->SetTitle("BGOfiducial l13 + rms + ts + psdstk + pc + sc + xtrl cut efficiency");

	BGOfiducial_eff_nBarLayer13_cut->SetTitle("BGOfiducial nBarLayer13 cut efficiency");
	BGOfiducial_eff_maxRms_cut->SetTitle("BGOfiducial maxRms cut efficiency");
	BGOfiducial_eff_track_selection_cut->SetTitle("BGOfiducial track selection cut efficiency");
	BGOfiducial_eff_psd_stk_match_cut->SetTitle("BGOfiducial PSD-STK match cut efficiency");
	BGOfiducial_eff_psd_charge_cut->SetTitle("BGOfiducial PSD charge cut efficiency");
	BGOfiducial_eff_stk_charge_cut->SetTitle("BGOfiducial STK charge cut efficiency");
	BGOfiducial_eff_xtrl_cut->SetTitle("BGOfiducial xtrl cut efficiency");
	BGOfiducial_eff_all_cut->SetTitle("BGOfiducial all cut efficiency");

	// Write histos to disk
	BGOfiducial_eff_l13_maxRms_cut->Write();
	BGOfiducial_eff_l13_rms_track_selection_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut->Write();
	BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut->Write();

	BGOfiducial_eff_nBarLayer13_cut->Write();
	BGOfiducial_eff_maxRms_cut->Write();
	BGOfiducial_eff_track_selection_cut->Write();
	BGOfiducial_eff_psd_stk_match_cut->Write();
	BGOfiducial_eff_psd_charge_cut->Write();
	BGOfiducial_eff_stk_charge_cut->Write();
	BGOfiducial_eff_xtrl_cut->Write();
	BGOfiducial_eff_all_cut->Write();

	auto geo_analysisDir = outFile.mkdir("Analysis_GeoCut");
	geo_analysisDir->cd();
	
	h_geo_BGOrec_slopeX.Write();
	h_geo_BGOrec_slopeY.Write();
	h_geo_BGOrec_interceptX.Write();
	h_geo_BGOrec_interceptY.Write();
	h_geo_BGOreco_topMap.Write();
	h_geo_BGOreco_bottomMap.Write();

	auto BGOdir = outFile.mkdir("BGO_Energy");
	BGOdir->cd();

	h_BGOrec_energy.Write();
	h_layer_max_energy_ratio.Write();

	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
		h_layer_energy_ratio[lIdx].Write();

	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
		(*it)->Write();

	sumRms_cosine_20_100.Write();
	sumRms_cosine_100_250.Write();
	sumRms_cosine_250_500.Write();
	sumRms_cosine_500_1000.Write();
	sumRms_cosine_1000_3000.Write();
	sumRms_cosine_3000_10000.Write();

	auto XTRLdir = outFile.mkdir("xtrl");
	XTRLdir->cd();

	h_xtrl_energy_int.Write();
	h_xtrl.Write();
	e_discrimination.Write();
	e_discrimination_20_100.Write();
	e_discrimination_100_250.Write();
	e_discrimination_250_500.Write();
	e_discrimination_500_1000.Write();
	e_discrimination_1000_3000.Write();
	e_discrimination_3000_10000.Write();

	e_discrimination_last.Write();
	e_discrimination_last_20_100.Write();
	e_discrimination_last_100_250.Write();
	e_discrimination_last_250_500.Write();
	e_discrimination_last_500_1000.Write();
	e_discrimination_last_1000_3000.Write();
	e_discrimination_last_3000_10000.Write();

	for (auto it = bin_xtrl.begin(); it != bin_xtrl.end(); ++it)
		(*it)->Write();

	auto psd_chargeDir = outFile.mkdir("PSDcharge");
	psd_chargeDir->cd();

	h_psd_chargeX.Write();
	h_psd_chargeY.Write();
	h_psd_charge.Write();
	h_psd_charge2D.Write();

	h_psd_selected_chargeX.Write();
	h_psd_selected_chargeY.Write();
	h_psd_selected_charge.Write();
	h_psd_selected_charge2D.Write();

	auto stk_chargeDir = outFile.mkdir("STKcharge");
	stk_chargeDir->cd();

	h_stk_chargeX.Write();
	h_stk_chargeY.Write();
	h_stk_charge.Write();
	h_stk_charge2D.Write();

	h_stk_selected_chargeX.Write();
	h_stk_selected_chargeY.Write();
	h_stk_selected_charge.Write();
	h_stk_selected_charge2D.Write();

	// Create ancillary output folder
	auto pDir = outFile.mkdir("proton_background");
	pDir->cd();

	h_background_under_xtrl_cut.Write();
	h_background_over_xtrl_cut.Write();

	// Create proton background ratio
	auto proton_background_ratio = static_cast<TH1D *>(h_background_under_xtrl_cut.Clone("proton_background_ratio"));
	proton_background_ratio->SetTitle("Proton background ratio");
	proton_background_ratio->Divide(&h_background_over_xtrl_cut);

	proton_background_ratio->Write();

	// Fill final vector
	std::vector<TH1D> data_selection_result (2, TH1D());
	data_selection_result[0] = h_BGOfiducial_all_cut_ce;
	data_selection_result[1] = h_background_over_xtrl_cut;

	// Create time output folder
	auto timeDir = outFile.mkdir("Time");
	timeDir->cd();
	
	h_second.Write();
	h_msecond.Write();

	return data_selection_result;
}