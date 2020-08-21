#include "acceptance.h"
#include "acceptance_cuts.h"
#include "data_cuts.h"
#include "energy_match.h"
#include "aggregate_events.h"
#include "wtsydp.h"
#include "BGO_energy_cuts.h"
#include "DAMPE_geo_structure.h"
#include "DmpBgoContainer.h"
#include "DmpPsdContainer.h"
#include "read_sets_config_file.h"
#include "binning.h"
#include "charge.h"
#include "mc_ancillary.h"
#include "fill_event_histo.h"

#include "TGraphErrors.h"
#include "TEfficiency.h"

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

inline std::shared_ptr<TH1D> buildHistoFromVector(
	const std::vector<double> &energyValues,
	const std::vector<double> &consgFactor)
{
	std::shared_ptr<TH1D> histo = std::make_shared<TH1D>("histo", "histoTitle", consgFactor.size(), energyValues[0], energyValues[energyValues.size() - 1]);
	for (auto idx = 1; idx <= histo->GetNbinsX(); ++idx)
		histo->SetBinContent(idx, consgFactor[idx - 1]);

	return histo;
}

inline void updateProcessStatus(const int evIdx, int &kStep, const int nevents)
{
	auto percentage = ((evIdx + 1) / (double)nevents) * 100;
	if (floor(percentage) != 0 && ((int)floor(percentage) % kStep) == 0)
	{
		std::cout << "\n"
				  << percentage << " %\t | \tProcessed " << evIdx + 1 << " events / " << nevents;
		kStep += 10;
	}
}

inline void load_particle_struct(
    p_type &simu_particle,
    const std::shared_ptr<DmpEvtSimuPrimaries> &simu_primaries)
{
	switch(simu_primaries->pvpart_pdg)
	{
		case 11:
			simu_particle.is_electron = true;
			break;
		case 2212:
			simu_particle.is_proton = true;
			break;
		default:
			std::cerr << "\n\nWrong particle ID\n\n";
	}	
}

void buildAcceptance(
	const std::string accInputPath,
	const bool verbose,
	TFile &outFile,
	const std::string wd)
{
	auto dmpch = aggregateEventsTChain(accInputPath, verbose);

	// Register Header container
	std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
	dmpch->SetBranchAddress("EventHeader", &evt_header);

	// Register SimuPrimaries container
	std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::make_shared<DmpEvtSimuPrimaries>();
	dmpch->SetBranchAddress("DmpEvtSimuPrimaries", &simu_primaries);

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

	// Create acceptance cuts struct
	cuts_conf acceptance_cuts;
	// Create active cuts struct
	data_active_cuts active_cuts;
	// Create ancillary cuts struct
	mc_ancillary_cuts ancillary_cuts;
	// Create particle type struct
	p_type simu_particle;

	// Load structs reading config file
	auto logEBins = load_acceptance_struct(
		acceptance_cuts,
		active_cuts,
		ancillary_cuts,
		wd);

	// Read dataSets connfig file
	data_set_conf input_sets;
	load_input_dsets_config(input_sets, wd);

	// Load MC statistics structure
	mc_statistics mc_selection;

	// Print filter status
	if (verbose)
		print_filter_status(active_cuts);
	
	// Event loop
	auto nevents = dmpch->GetEntries();
	if (verbose)
		std::cout << "Total number of events: " << nevents << "\n\n";

	// Cut histos
	TH1D h_geo_factor("h_geo_factor", "Energy Distribution of the geometric factor; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_incoming("h_incoming", "Energy Distribution of the incoming particles; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_trigger("h_trigger", "Energy Distribution of the triggered particles; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_cut("h_geometric_cut", "Energy Distribution - geometric (trigger selection) cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_maxElayer_cut("h_maxElayer_cut", "Energy Distribution - maxElayer cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_maxBarLayer_cut("h_maxBarLayer_cut", "Energy Distribution - maxBarLayer cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOTrackContainment_cut("h_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGO_fiducial_cut("h_BGO_fiducial_cut", "Energy Distibution - BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_all_cut("h_all_cut", "Energy Distribution - All cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	TH1D h_geo_factor_nw("h_geo_factor_nw", "Energy Distribution of the geometric factor; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_incoming_nw("h_incoming_nw", "Energy Distribution of the incoming particles; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_trigger_nw("h_trigger_nw", "Energy Distribution of the triggered particles; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_cut_nw("h_geometric_cut_nw", "Energy Distribution - geometric (trigger selection) cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_maxElayer_cut_nw("h_maxElayer_cut_nw", "Energy Distribution - maxElayer cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_maxBarLayer_cut_nw("h_maxBarLayer_cut_nw", "Energy Distribution - maxBarLayer cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOTrackContainment_cut_nw("h_BGOTrackContainment_cut_nw", "Energy Distribution - BGOTrackContainment cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGO_fiducial_cut_nw("h_BGO_fiducial_cut_nw", "Energy Distibution - BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_all_cut_nw("h_all_cut_nw", "Energy Distribution - All cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// Cuts && Geometric Cut
	TH1D h_geometric_maxElayer_cut("h_geometric_maxElayer_cut", "Energy Distribution - maxElayer + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_maxBarLayer_cut("h_geometric_maxBarLayer_cut", "Energy Distribution - maxBarLayer + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_BGOTrackContainment_cut("h_geometric_BGOTrackContainment_cut", "Energy Distribution - BGOTrackContainment + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_BGO_fiducial_cut("h_geometric_BGO_fiducial_cut", "Energy Distibution - BGO fiducial + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_all_cut("h_geometric_all_cut", "Energy Distribution - All + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	TH1D h_geometric_maxElayer_cut_nw("h_geometric_maxElayer_cut_nw", "Energy Distribution - maxElayer + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_maxBarLayer_cut_nw("h_geometric_maxBarLayer_cut_nw", "Energy Distribution - maxBarLayer + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_BGOTrackContainment_cut_nw("h_geometric_BGOTrackContainment_cut_nw", "Energy Distribution - BGOTrackContainment + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_BGO_fiducial_cut_nw("h_geometric_BGO_fiducial_cut_nw", "Energy Distibution - BGO fiducial + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_geometric_all_cut_nw("h_geometric_all_cut_nw", "Energy Distribution - All + geometric cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// Cuts && BGO fiducial volume cut
	TH1D h_BGOfiducial_nBarLayer13_cut("h_BGOfiducial_nBarLayer13_cut", "Energy Distribution - nBarLayer13 + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_maxRms_cut("h_BGOfiducial_maxRms_cut", "Energy Distribution - maxRms  + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_track_selection_cut("h_BGOfiducial_track_selection_cut", "Energy Distribution - track selection + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_psd_stk_match_cut("h_BGOfiducial_psd_stk_match_cut", "Energy Distribution - psd-stk match + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_psd_charge_cut("h_BGOfiducial_psd_charge_cut", "Energy Distribution - psd charge + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_stk_charge_cut("h_BGOfiducial_stk_charge_cut", "Energy Distribution - stk charge + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_xtrl_cut("h_BGOfiducial_xtrl_cut", "Energy Distribution - xtrl + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_all_cut("h_BGOfiducial_all_cut", "Energy Distribution - All + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	TH1D h_BGOfiducial_nBarLayer13_cut_nw("h_BGOfiducial_nBarLayer13_cut_nw", "Energy Distribution - nBarLayer13 + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_maxRms_cut_nw("h_BGOfiducial_maxRms_cut_nw", "Energy Distribution - maxRms  + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_track_selection_cut_nw("h_BGOfiducial_track_selection_cut_nw", "Energy Distribution - track selection + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_psd_stk_match_cut_nw("h_BGOfiducial_psd_stk_match_cut_nw", "Energy Distribution - psd-stk match + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_psd_charge_cut_nw("h_BGOfiducial_psd_charge_cut_nw", "Energy Distribution - psd charge + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_stk_charge_cut_nw("h_BGOfiducial_stk_charge_cut_nw", "Energy Distribution - stk charge + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_xtrl_cut_nw("h_BGOfiducial_xtrl_cut_nw", "Energy Distribution - xtrl + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_BGOfiducial_all_cut_nw("h_BGOfiducial_all_cut_nw", "Energy Distribution - All + BGO fiducial cut; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));

	// Analysis histos - simu and reco energy of incoming events
	auto energy_ratio_line_binning = createLinearBinning(0, 1, 100);

	TH1D h_BGOrec_energy("h_BGOrec_energy", "BGO Energy: Raw Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_simu_energy("h_simu_energy", "Simu Energy; Real Energy (GeV); counts", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_energy_diff("h_energy_diff", "Simu vs Corrected Reco BGO energy: Real Energy - Corrected Energy (GeV); counts", 100, 0, 1);
	TH2D h_energy_diff2D("h_energy_diff2D", "Energy Ratio; Real Energy (GeV); (Real - Raw)/Raw", logEBins.size() - 1, &(logEBins[0]), energy_ratio_line_binning.size() - 1, &(energy_ratio_line_binning[0]));
	TH2D h_energy_unfold("h_energy_unfold", "Energy Unfolding Matrix; Real Energy (GeV); Raw Energy (GeV)", logEBins.size() - 1, &(logEBins[0]), logEBins.size() - 1, &(logEBins[0]));

	// Pre Geometric Cut
	// Top X and Y
	TH1D h_preGeo_BGOrec_topX_vs_realX("h_preGeo_BGOrec_topX_vs_realX", "Real X - BGOrec TOP X", 100, -100, 100);
	TH1D h_preGeo_BGOrec_topY_vs_realY("h_preGeo_BGOrec_topY_vs_realY", "Real Y - BGOrec TOP Y", 100, -100, 100);

	// Slope X and Y
	TH1D h_preGeo_real_slopeX("h_preGeo_real_slopeX", "Real Slope X", 1000, -90, 90);
	TH1D h_preGeo_real_slopeY("h_preGeo_real_slopeY", "Real Slope Y", 1000, -90, 90);
	TH1D h_preGeo_BGOrec_slopeX("h_preGeo_BGOrec_slopeX", "BGOrec Slope X", 1000, -90, 90);
	TH1D h_preGeo_BGOrec_slopeY("h_preGeo_BGOrec_slopeY", "BGOrec Slope Y", 1000, -90, 90);

	// Intercept X and Y
	TH1D h_preGeo_real_interceptX("h_preGeo_real_interceptX", "Real Intercept X", 500, -500, 500);
	TH1D h_preGeo_real_interceptY("h_preGeo_real_interceptY", "Real Intercept Y", 500, -500, 500);
	TH1D h_preGeo_BGOrec_interceptX("h_preGeo_BGOrec_interceptX", "BGOrec Intercept X", 500, -500, 500);
	TH1D h_preGeo_BGOrec_interceptY("h_preGeo_BGOrec_interceptY", "BGOrec Intercept Y", 500, -500, 500);

	// Top Maps
	TH2D h_preGeo_real_topMap("h_preGeo_real_topMap", "Real BGO TOP Map", 500, -500, 500, 500, -500, 500);
	TH2D h_preGeo_BGOreco_topMap("h_preGeo_BGOreco_topMap", "BGOreco TOP Map", 500, -500, 500, 500, -500, 500);

	// Bottom Maps
	TH2D h_preGeo_real_bottomMap("h_preGeo_real_bottomMap", "Real BGO BOTTOM Map", 500, -500, 500, 500, -500, 500);
	TH2D h_preGeo_BGOreco_bottomMap("h_preGeo_BGOreco_bottomMap", "BGOreco BOTTOM Map", 500, -500, 500, 500, -500, 500);

	// Map of events outside the "real" first BGO layer
	TH2D h_noBGOenergy_real_topMap("h_noBGOenergy_real_topMap", "Real BGO TOP Map", 500, -500, 500, 500, -500, 500);

	// After Geometric Cut
	// Top X and Y
	TH1D h_geo_BGOrec_topX_vs_realX("h_geo_BGOrec_topX_vs_realX", "Real X - BGOrec TOP X", 100, -100, 100);
	TH1D h_geo_BGOrec_topY_vs_realY("h_geo_BGOrec_topY_vs_realY", "Real Y - BGOrec TOP Y", 100, -100, 100);

	// Slope X and Y
	TH1D h_geo_real_slopeX("h_geo_real_slopeX", "Real Slope X", 1000, -90, 90);
	TH1D h_geo_real_slopeY("h_geo_real_slopeY", "Real Slope Y", 1000, -90, 90);
	TH1D h_geo_BGOrec_slopeX("h_geo_BGOrec_slopeX", "BGOrec Slope X", 1000, -90, 90);
	TH1D h_geo_BGOrec_slopeY("h_geo_BGOrec_slopeY", "BGOrec Slope Y", 1000, -90, 90);

	// Intercept X and Y
	TH1D h_geo_real_interceptX("h_geo_real_interceptX", "Real Intercept X", 500, -500, 500);
	TH1D h_geo_real_interceptY("h_geo_real_interceptY", "Real Intercept Y", 500, -500, 500);
	TH1D h_geo_BGOrec_interceptX("h_geo_BGOrec_interceptX", "BGOrec Intercept X", 500, -500, 500);
	TH1D h_geo_BGOrec_interceptY("h_geo_BGOrec_interceptY", "BGOrec Intercept Y", 500, -500, 500);

	// Top Maps
	TH2D h_geo_real_topMap("h_geo_real_topMap", "Real BGO TOP Map", 500, -500, 500, 500, -500, 500);
	TH2D h_geo_BGOreco_topMap("h_geo_BGOreco_topMap", "BGOreco TOP Map", 500, -500, 500, 500, -500, 500);

	// Bottom Maps
	TH2D h_geo_real_bottomMap("h_geo_real_bottomMap", "Real BGO BOTTOM Map", 500, -500, 500, 500, -500, 500);
	TH2D h_geo_BGOreco_bottomMap("h_geo_BGOreco_bottomMap", "BGOreco BOTTOM Map", 500, -500, 500, 500, -500, 500);

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
	TH1D h_background_under_xtrl_cut("h_background_under_xtrl_cut", "Proton background - XTRL < cut", logEBins.size() - 1, &(logEBins[0]));
	TH1D h_background_over_xtrl_cut("h_background_over_xtrl_cut", "Proton background - XTRL > cut", logEBins.size() - 1, &(logEBins[0]));

	// Sumw2 Acceptance - First-Cut histos
	h_geo_factor.Sumw2();
	h_incoming.Sumw2();
	h_trigger.Sumw2();
	h_geometric_cut.Sumw2();
	h_maxElayer_cut.Sumw2();
	h_maxBarLayer_cut.Sumw2();
	h_BGOTrackContainment_cut.Sumw2();
	h_BGO_fiducial_cut.Sumw2();
	h_all_cut.Sumw2();

	h_geo_factor_nw.Sumw2();
	h_incoming_nw.Sumw2();
	h_trigger_nw.Sumw2();
	h_geometric_cut_nw.Sumw2();
	h_maxElayer_cut_nw.Sumw2();
	h_maxBarLayer_cut_nw.Sumw2();
	h_BGOTrackContainment_cut_nw.Sumw2();
	h_BGO_fiducial_cut_nw.Sumw2();
	h_all_cut_nw.Sumw2();

	// Sumw2 Acceptance - Cuts && Geometric Cut
	h_geometric_maxElayer_cut.Sumw2();
	h_geometric_maxBarLayer_cut.Sumw2();
	h_geometric_BGOTrackContainment_cut.Sumw2();
	h_geometric_BGO_fiducial_cut.Sumw2();
	h_geometric_all_cut.Sumw2();

	h_geometric_maxElayer_cut_nw.Sumw2();
	h_geometric_maxBarLayer_cut_nw.Sumw2();
	h_geometric_BGOTrackContainment_cut_nw.Sumw2();
	h_geometric_BGO_fiducial_cut_nw.Sumw2();
	h_geometric_all_cut_nw.Sumw2();

	// Sumw2 Acceptance - Cuts && BGO fiducial volume cut
	h_BGOfiducial_nBarLayer13_cut.Sumw2();
	h_BGOfiducial_maxRms_cut.Sumw2();
	h_BGOfiducial_track_selection_cut.Sumw2();
	h_BGOfiducial_psd_stk_match_cut.Sumw2();
	h_BGOfiducial_psd_charge_cut.Sumw2();
	h_BGOfiducial_stk_charge_cut.Sumw2();
	h_BGOfiducial_xtrl_cut.Sumw2();
	h_BGOfiducial_all_cut.Sumw2();

	h_BGOfiducial_nBarLayer13_cut_nw.Sumw2();
	h_BGOfiducial_maxRms_cut_nw.Sumw2();
	h_BGOfiducial_track_selection_cut_nw.Sumw2();
	h_BGOfiducial_psd_stk_match_cut_nw.Sumw2();
	h_BGOfiducial_psd_charge_cut_nw.Sumw2();
	h_BGOfiducial_stk_charge_cut_nw.Sumw2();
	h_BGOfiducial_xtrl_cut_nw.Sumw2();
	h_BGOfiducial_all_cut_nw.Sumw2();

	// Sumw2 Analysis histos - simu and reco energy of incoming events
	h_BGOrec_energy.Sumw2();
	h_simu_energy.Sumw2();
	h_energy_diff.Sumw2();
	h_energy_diff2D.Sumw2();
	h_energy_unfold.Sumw2();

	// Sumw2 Analysis histos - preGeo
	h_preGeo_BGOrec_topX_vs_realX.Sumw2();
	h_preGeo_BGOrec_topY_vs_realY.Sumw2();
	h_preGeo_real_slopeX.Sumw2();
	h_preGeo_real_slopeY.Sumw2();
	h_preGeo_BGOrec_slopeX.Sumw2();
	h_preGeo_BGOrec_slopeY.Sumw2();
	h_preGeo_real_interceptX.Sumw2();
	h_preGeo_real_interceptY.Sumw2();
	h_preGeo_BGOrec_interceptX.Sumw2();
	h_preGeo_BGOrec_interceptY.Sumw2();
	h_preGeo_real_topMap.Sumw2();
	h_preGeo_BGOreco_topMap.Sumw2();
	h_preGeo_real_bottomMap.Sumw2();
	h_preGeo_BGOreco_bottomMap.Sumw2();
	h_noBGOenergy_real_topMap.Sumw2();

	// Sumw2 Analysis histos - Geo
	h_geo_BGOrec_topX_vs_realX.Sumw2();
	h_geo_BGOrec_topY_vs_realY.Sumw2();
	h_geo_real_slopeX.Sumw2();
	h_geo_real_slopeY.Sumw2();
	h_geo_BGOrec_slopeX.Sumw2();
	h_geo_BGOrec_slopeY.Sumw2();
	h_geo_real_interceptX.Sumw2();
	h_geo_real_interceptY.Sumw2();
	h_geo_BGOrec_interceptX.Sumw2();
	h_geo_BGOrec_interceptY.Sumw2();
	h_geo_real_topMap.Sumw2();
	h_geo_BGOreco_topMap.Sumw2();
	h_geo_real_bottomMap.Sumw2();
	h_geo_BGOreco_bottomMap.Sumw2();
	h_layer_max_energy_ratio.Sumw2();

	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
		(*it)->Sumw2();

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

	for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
	{
		// Get chain event
		dmpch->GetEvent(evIdx);

		// Read the particle ID
		load_particle_struct(
			simu_particle, 
			simu_primaries);

		// Update event counter
		++mc_selection.event_counter;

		// Event printout
		if (verbose)
			updateProcessStatus(evIdx, kStep, nevents);
		
		// Get event total energy
		double bgoTotalE_raw = bgorec->GetTotalEnergy(); 	// Energy in MeV - not correcte
		double simuEnergy = simu_primaries->pvpart_ekin; 	//Energy of simu primaries particle (MeV)

		double energy_w;
		if (simu_particle.is_electron)
			energy_w = pow(simuEnergy * _GeV, -2);
		else if (simu_particle.is_proton)
			energy_w = pow(simuEnergy * _GeV, -1.7);

		// Don't accept events outside the selected energy window
		if (simuEnergy * _GeV < acceptance_cuts.min_event_energy || simuEnergy * _GeV > acceptance_cuts.max_event_energy)
		{
			++mc_selection.generated_events_out_range;
			continue;
		}

		++mc_selection.generated_events_in_range;

		if (geometric_cut(simu_primaries))
		{
			h_geo_factor.Fill(simuEnergy * _GeV, energy_w);
			h_geo_factor_nw.Fill(simuEnergy * _GeV);
		}

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
			++mc_selection.triggered_events;
			h_trigger.Fill(simuEnergy * _GeV, energy_w);
			h_trigger_nw.Fill(simuEnergy * _GeV);

			if (checkBGOreco(bgorec, simu_primaries))
			{
				h_incoming.Fill(simuEnergy * _GeV, energy_w);
				h_incoming_nw.Fill(simuEnergy * _GeV);
			}
			else
				continue;
		}
		else
		{
			// Increase the incoming events
			h_incoming.Fill(simuEnergy * _GeV, energy_w);
			h_incoming_nw.Fill(simuEnergy * _GeV);

			// Evaluate the position on the First BGO layer for the non - triggered events
			evaluateTopBottomPosition(
				simu_primaries,
				bgorec,
				h_preGeo_BGOrec_topX_vs_realX,
				h_preGeo_BGOrec_topY_vs_realY,
				h_preGeo_real_slopeX,
				h_preGeo_real_slopeY,
				h_preGeo_BGOrec_slopeX,
				h_preGeo_BGOrec_slopeY,
				h_preGeo_real_interceptX,
				h_preGeo_real_interceptY,
				h_preGeo_BGOrec_interceptX,
				h_preGeo_BGOrec_interceptY,
				h_preGeo_real_topMap,
				h_preGeo_BGOreco_topMap,
				h_preGeo_real_bottomMap,
				h_preGeo_BGOreco_bottomMap,
				energy_w);

			continue;
		}

		// Fill the energy histos only for good reco events
		h_BGOrec_energy.Fill(bgoTotalE_raw * _GeV, energy_w);
		h_simu_energy.Fill(simuEnergy * _GeV, energy_w);
		h_energy_diff.Fill((simuEnergy - bgoTotalE_raw) / simuEnergy, energy_w);
		h_energy_diff2D.Fill(simuEnergy * _GeV, (simuEnergy - bgoTotalE_raw) / simuEnergy, energy_w);
		h_energy_unfold.Fill(simuEnergy * _GeV, bgoTotalE_raw * _GeV, energy_w);

		// Create best_track struct
		best_track event_best_track;
		
		// Create PSD-STK match struct
		psd_cluster_match clu_matching;

		// Load BGO event class
		DmpBgoContainer bgoVault;

		// Load PSD event class
		DmpPsdContainer psdVault;
		
		bgoVault.scanBGOHits(
			bgohits,
			bgoTotalE_raw,
			acceptance_cuts);

		psdVault.scanPSDHits(
			psdhits,
			acceptance_cuts);

		// Create PSD charge struct
		psd_charge extracted_psd_charge;

		// Create STK charge struct
		stk_charge extracted_stk_charge;

		// evaluate the energy raio on each single layer of the BGO
		evaluateEnergyRatio(
			bgorec,
			acceptance_cuts,
			bgoTotalE_raw,
			h_layer_max_energy_ratio,
			h_layer_energy_ratio,
			energy_w);

		// Event filter
		event_filter filter;
		if (filter_this_mc_event(
				filter,
				simu_primaries,
				bgorec,
				bgohits,
				acceptance_cuts,
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
			++mc_selection.selected_events;

		// Fill cuts histos

		// Fill geometric cut histos
		if (filter.geometric)
		{
			h_geometric_cut.Fill(simuEnergy * _GeV, energy_w);
			h_geometric_cut_nw.Fill(simuEnergy * _GeV);
			
			// Evaluate the position on the First BGO layer (after geometric cut)
			evaluateTopBottomPosition_data(
				bgorec,
				h_geo_BGOrec_slopeX,
				h_geo_BGOrec_slopeY,
				h_geo_BGOrec_interceptX,
				h_geo_BGOrec_interceptY,
				h_geo_BGOreco_topMap,
				h_geo_BGOreco_bottomMap,
				energy_w);

			// Geometric cut && maxElayer cut
			if (filter.BGO_fiducial_maxElayer_cut)
			{
				h_geometric_maxElayer_cut.Fill(simuEnergy * _GeV, energy_w);
				h_geometric_maxElayer_cut_nw.Fill(simuEnergy * _GeV);
			}

			// Geometric cut && maxBarLayer cut
			if (filter.BGO_fiducial_maxBarLayer_cut)
			{
				h_geometric_maxBarLayer_cut.Fill(simuEnergy * _GeV, energy_w);
				h_geometric_maxBarLayer_cut_nw.Fill(simuEnergy * _GeV);
			}

			// Geometric cut && BGOTrackContainment cut
			if (filter.BGO_fiducial_BGOTrackContainment_cut)
			{
				h_geometric_BGOTrackContainment_cut.Fill(simuEnergy * _GeV, energy_w);
				h_geometric_BGOTrackContainment_cut_nw.Fill(simuEnergy * _GeV);
			}

			// Geometric cut && BGO fiducial cut
			if (filter.BGO_fiducial)
			{
				h_geometric_BGO_fiducial_cut.Fill(simuEnergy * _GeV, energy_w);
				h_geometric_BGO_fiducial_cut_nw.Fill(simuEnergy * _GeV);
			}
			
			// Geometric cut and all cuts
			if (filter.all_cut)
			{
				h_geometric_all_cut.Fill(simuEnergy * _GeV, energy_w);
				h_geometric_all_cut_nw.Fill(simuEnergy * _GeV);
			}
		}

		// Fill BGO_fiducial_maxElayer cut histo
		if (filter.BGO_fiducial_maxElayer_cut)
		{
			h_maxElayer_cut.Fill(simuEnergy * _GeV, energy_w);
			h_maxElayer_cut_nw.Fill(simuEnergy * _GeV);
		}
		
		// Fill BGO_fiducial_maxBarLayer cut histo
		if (filter.BGO_fiducial_maxBarLayer_cut)
		{
			h_maxBarLayer_cut.Fill(simuEnergy * _GeV, energy_w);
			h_maxBarLayer_cut_nw.Fill(simuEnergy * _GeV);
		}
		
		// Fill BGO_fiducial_BGOTrackContainment cut histo
		if (filter.BGO_fiducial_BGOTrackContainment_cut)
		{
			h_BGOTrackContainment_cut.Fill(simuEnergy * _GeV, energy_w);
			h_BGOTrackContainment_cut_nw.Fill(simuEnergy * _GeV);
		}

		// Fill BGO fiducial volume cut
		if (filter.BGO_fiducial)
		{
			h_BGO_fiducial_cut.Fill(simuEnergy * _GeV, energy_w);
			h_BGO_fiducial_cut_nw.Fill(simuEnergy * _GeV);

			// BGO fiducial cut && nBarLayer13 cut
			if (filter.nBarLayer13_cut)
			{
				h_BGOfiducial_nBarLayer13_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_nBarLayer13_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && maxRms cut
			if (filter.maxRms_cut)
			{
				h_BGOfiducial_maxRms_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_maxRms_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && track selection cut
			if (filter.track_selection_cut)
			{
				h_BGOfiducial_track_selection_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_track_selection_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && PSD-STK match cut
			if (filter.psd_stk_match_cut)
			{
				h_BGOfiducial_psd_stk_match_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_psd_stk_match_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && PSD charge cut
			if (filter.psd_charge_cut)
			{
				h_BGOfiducial_psd_charge_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_psd_charge_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && STK charge cut
			if (filter.stk_charge_cut)
			{
				h_BGOfiducial_stk_charge_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_stk_charge_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && XTRL cut
			if (filter.xtrl_cut)
			{
				h_BGOfiducial_xtrl_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_xtrl_cut_nw.Fill(simuEnergy * _GeV);
			}

			// BGO fiducial cut && all cut
			if (filter.all_cut)
			{
				h_BGOfiducial_all_cut.Fill(simuEnergy * _GeV, energy_w);
				h_BGOfiducial_all_cut_nw.Fill(simuEnergy * _GeV);
			}
		}

		// Collect PSD particle charge
		if (filter.psd_charge_measurement)
			fillChargeHistos<psd_charge>(
				h_psd_chargeX,
				h_psd_chargeY,
				h_psd_charge,
				h_psd_charge2D,
				extracted_psd_charge,
				energy_w);
		
		// Collect selected STK particle charge
		if (filter.psd_charge_cut)
			fillChargeHistos<psd_charge>(
				h_psd_selected_chargeX,
				h_psd_selected_chargeY,
				h_psd_selected_charge,
				h_psd_selected_charge2D,
				extracted_psd_charge,
				energy_w);

		// Collect STK particle charge
		if (filter.stk_charge_measurement)
			fillChargeHistos<stk_charge>(
				h_stk_chargeX,
				h_stk_chargeY,
				h_stk_charge,
				h_stk_charge2D,
				extracted_stk_charge,
				energy_w);

		// Collect selected STK particle charge
		if (filter.stk_charge_cut)
			fillChargeHistos<stk_charge>(
				h_stk_selected_chargeX,
				h_stk_selected_chargeY,
				h_stk_selected_charge,
				h_stk_selected_charge2D,
				extracted_stk_charge,
				energy_w);

		// Fill all cut histo
		if (active_cuts.nActiveCuts)
		{
			if (filter.all_cut)
			{
				h_all_cut.Fill(simuEnergy * _GeV, energy_w);
				
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
					simuEnergy,
					bin_xtrl,
					h_xtrl_energy_int,
					h_xtrl,
					energy_w);

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
					simuEnergy);

				// Compute proton background
				if (ancillary_cuts.compute_proton_background)
					compute_proton_background(
						bgoVault.GetSumRMS(),
						bgoVault.GetLastFFracLayer(),
						simuEnergy,
						acceptance_cuts,
						h_background_under_xtrl_cut,
						h_background_over_xtrl_cut,
						energy_w);
			}
		}
	}
	
	if (verbose)
	{
		std::cout << "\n\n ****** \n\n";
		std::cout << "Generated events: " << mc_selection.event_counter << std::endl;
		std::cout << "Generated events in energy range: " << mc_selection.generated_events_in_range << std::endl;
		std::cout << "Generated events out of energy range: " << mc_selection.generated_events_out_range << std::endl;
		std::cout << "Triggered events: " << mc_selection.triggered_events << std::endl;
		std::cout << "\n\n**** Filter result ****\n";
		std::cout << "***********************\n\n";
		std::cout << "Particles surviving the selection cuts: " << mc_selection.selected_events << "\t | " << ((double)mc_selection.selected_events / mc_selection.triggered_events) * 100 << "%";
		std::cout << "\n\n ***********************\n\n";
	}

	double genSurface = 4 * TMath::Pi() * pow(acceptance_cuts.vertex_radius, 2) / 2;
	double scaleFactor = TMath::Pi() * genSurface;

	// Building acceptance histos
	auto h_acceptance_geometric_factor = static_cast<TH1D *>(h_geo_factor.Clone("h_acceptance_geometric_factor"));
	auto h_acceptance_gometric_cut = static_cast<TH1D *>(h_geometric_cut.Clone("h_acceptance_gometric_cut"));
	auto h_acceptance_maxElayer_cut = static_cast<TH1D *>(h_maxElayer_cut.Clone("h_acceptance_maxElayer_cut"));
	auto h_acceptance_maxBarLayer_cut = static_cast<TH1D *>(h_maxBarLayer_cut.Clone("h_acceptance_maxBarLayer_cut"));
	auto h_acceptance_BGOTrackContainment_cut = static_cast<TH1D *>(h_BGOTrackContainment_cut.Clone("h_acceptance_BGOTrackContainment_cut"));
	auto h_acceptance_BGO_fiducial_cut = static_cast<TH1D *>(h_BGO_fiducial_cut.Clone("h_acceptance_BGO_fiducial_cut"));
	auto h_acceptance_BGO_fiducial_nBarLayer13_cut = static_cast<TH1D *>(h_BGOfiducial_nBarLayer13_cut.Clone("h_acceptance_BGO_fiducial_nBarLayer13_cut"));
	auto h_acceptance_BGO_fiducial_maxRms_cut = static_cast<TH1D *>(h_BGOfiducial_maxRms_cut.Clone("h_acceptance_BGO_fiducial_maxRms_cut"));
	auto h_acceptance_BGO_fiducial_track_selection_cut = static_cast<TH1D *>(h_BGOfiducial_track_selection_cut.Clone("h_acceptance_BGO_fiducial_track_selection_cut"));
	auto h_acceptance_BGO_fiducial_psd_stk_match_cut = static_cast<TH1D *>(h_BGOfiducial_psd_stk_match_cut.Clone("h_acceptance_BGO_fiducial_psd_stk_match_cut"));
	auto h_acceptance_BGO_fiducial_psd_charge_cut = static_cast<TH1D *>(h_BGOfiducial_psd_charge_cut.Clone("h_acceptance_BGO_fiducial_psd_charge_cut"));
	auto h_acceptance_BGO_fiducial_stk_charge_cut = static_cast<TH1D *>(h_BGOfiducial_stk_charge_cut.Clone("h_acceptance_BGO_fiducial_stk_charge_cut"));
	auto h_acceptance_BGO_fiducial_xtrl_cut = static_cast<TH1D *>(h_BGOfiducial_xtrl_cut.Clone("h_acceptance_BGO_fiducial_xtrl_cut"));
	auto h_acceptance_all_cut = static_cast<TH1D *>(h_all_cut.Clone("h_acceptance_all_cut"));

	h_acceptance_geometric_factor->Divide(&h_incoming);
	h_acceptance_gometric_cut->Divide(&h_incoming);
	h_acceptance_maxElayer_cut->Divide(&h_incoming);
	h_acceptance_maxBarLayer_cut->Divide(&h_incoming);
	h_acceptance_BGOTrackContainment_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_nBarLayer13_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_maxRms_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_track_selection_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_psd_stk_match_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_psd_charge_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_stk_charge_cut->Divide(&h_incoming);
	h_acceptance_BGO_fiducial_xtrl_cut->Divide(&h_incoming);
	h_acceptance_all_cut->Divide(&h_incoming);

	h_acceptance_geometric_factor->Scale(scaleFactor);
	h_acceptance_gometric_cut->Scale(scaleFactor);
	h_acceptance_maxElayer_cut->Scale(scaleFactor);
	h_acceptance_maxBarLayer_cut->Scale(scaleFactor);
	h_acceptance_BGOTrackContainment_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_nBarLayer13_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_maxRms_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_track_selection_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_psd_stk_match_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_psd_charge_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_stk_charge_cut->Scale(scaleFactor);
	h_acceptance_BGO_fiducial_xtrl_cut->Scale(scaleFactor);
	h_acceptance_all_cut->Scale(scaleFactor);

	// Builing vectors
	std::vector<double> energyValues(h_incoming.GetXaxis()->GetNbins(), 0);

	std::vector<double> acceptanceValues_geometric_factor(energyValues.size(), 0);
	std::vector<double> acceptanceValues_gometric_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_maxElayer_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_maxBarLayer_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGOTrackContainment_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_nBarLayer13_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_maxRms_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_track_selection_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_psd_stk_match_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_psd_charge_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_stk_charge_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_BGO_fiducial_xtrl_cut(energyValues.size(), 0);
	std::vector<double> acceptanceValues_all_cut(energyValues.size(), 0);

	//Building histo errors on energy and
	std::vector<double> acceptanceError_geometric_factor(h_incoming.GetXaxis()->GetNbins(), 0);
	std::vector<double> acceptanceError_gometric_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_maxElayer_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_maxBarLayer_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGOTrackContainment_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_nBarLayer13_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_maxRms_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_track_selection_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_psd_stk_match_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_psd_charge_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_stk_charge_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_BGO_fiducial_xtrl_cut(acceptanceError_geometric_factor.size(), 0);
	std::vector<double> acceptanceError_all_cut(acceptanceError_geometric_factor.size(), 0);

	std::vector<double> energyError(energyValues.size(), 0);

	for (auto it = logEBins.begin(); it != (logEBins.end() - 1); ++it)
	{
		auto index = std::distance(logEBins.begin(), it);
		energyValues[index] = wtsydp(*it, *(it + 1), getInputPowerLawIndex(*it, *(it + 1), input_sets));

		acceptanceValues_geometric_factor[index] = h_acceptance_geometric_factor->GetBinContent(index + 1);
		acceptanceValues_gometric_cut[index] = h_acceptance_gometric_cut->GetBinContent(index + 1);
		acceptanceValues_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinContent(index + 1);
		acceptanceValues_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinContent(index + 1);
		acceptanceValues_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_nBarLayer13_cut[index] = h_acceptance_BGO_fiducial_nBarLayer13_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_maxRms_cut[index] = h_acceptance_BGO_fiducial_maxRms_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_track_selection_cut[index] = h_acceptance_BGO_fiducial_track_selection_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_psd_stk_match_cut[index] = h_acceptance_BGO_fiducial_psd_stk_match_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_psd_charge_cut[index] = h_acceptance_BGO_fiducial_psd_charge_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_stk_charge_cut[index] = h_acceptance_BGO_fiducial_stk_charge_cut->GetBinContent(index + 1);
		acceptanceValues_BGO_fiducial_xtrl_cut[index] = h_acceptance_BGO_fiducial_xtrl_cut->GetBinContent(index + 1);
		acceptanceValues_all_cut[index] = h_acceptance_all_cut->GetBinContent(index + 1);

		acceptanceError_geometric_factor[index] = h_acceptance_geometric_factor->GetBinError(index + 1);
		acceptanceError_gometric_cut[index] = h_acceptance_gometric_cut->GetBinError(index + 1);
		acceptanceError_maxElayer_cut[index] = h_acceptance_maxElayer_cut->GetBinError(index + 1);
		acceptanceError_maxBarLayer_cut[index] = h_acceptance_maxBarLayer_cut->GetBinError(index + 1);
		acceptanceError_BGOTrackContainment_cut[index] = h_acceptance_BGOTrackContainment_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_cut[index] = h_acceptance_BGO_fiducial_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_nBarLayer13_cut[index] = h_acceptance_BGO_fiducial_nBarLayer13_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_maxRms_cut[index] = h_acceptance_BGO_fiducial_maxRms_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_track_selection_cut[index] = h_acceptance_BGO_fiducial_track_selection_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_psd_stk_match_cut[index] = h_acceptance_BGO_fiducial_psd_stk_match_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_psd_charge_cut[index] = h_acceptance_BGO_fiducial_psd_charge_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_stk_charge_cut[index] = h_acceptance_BGO_fiducial_stk_charge_cut->GetBinError(index + 1);
		acceptanceError_BGO_fiducial_xtrl_cut[index] = h_acceptance_BGO_fiducial_xtrl_cut->GetBinError(index + 1);
		acceptanceError_all_cut[index] = h_acceptance_all_cut->GetBinError(index + 1);
	}

	// Building graphs
	TGraphErrors gr_acceptance_geometric_factor(energyValues.size(), &(energyValues[0]), &(acceptanceValues_geometric_factor[0]), &(energyError[0]), &(acceptanceError_geometric_factor[0]));
	TGraphErrors gr_acceptance_gometric_cut(energyValues.size(), &energyValues[0], &acceptanceValues_gometric_cut[0], &(energyError[0]), &(acceptanceError_gometric_cut[0]));
	TGraphErrors gr_acceptance_maxElayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxElayer_cut[0], &(energyError[0]), &(acceptanceError_maxElayer_cut[0]));
	TGraphErrors gr_acceptance_maxBarLayer_cut(energyValues.size(), &energyValues[0], &acceptanceValues_maxBarLayer_cut[0], &(energyError[0]), &(acceptanceError_maxBarLayer_cut[0]));
	TGraphErrors gr_acceptance_BGOTrackContainment_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGOTrackContainment_cut[0], &(energyError[0]), &(acceptanceError_BGOTrackContainment_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_nBarLayer13_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_nBarLayer13_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_nBarLayer13_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_maxRms_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_maxRms_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_maxRms_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_track_selection_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_track_selection_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_track_selection_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_psd_stk_match_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_psd_stk_match_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_psd_stk_match_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_psd_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_psd_charge_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_psd_charge_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_stk_charge_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_stk_charge_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_stk_charge_cut[0]));
	TGraphErrors gr_acceptance_BGO_fiducial_xtrl_cut(energyValues.size(), &energyValues[0], &acceptanceValues_BGO_fiducial_xtrl_cut[0], &(energyError[0]), &(acceptanceError_BGO_fiducial_xtrl_cut[0]));
	TGraphErrors gr_acceptance_all_cut(energyValues.size(), &energyValues[0], &acceptanceValues_all_cut[0], &(energyError[0]), &(acceptanceError_all_cut[0]));

	gr_acceptance_geometric_factor.SetName("gr_acceptance_geometric_factor");
	gr_acceptance_gometric_cut.SetName("gr_acceptance_gometric_cut");
	gr_acceptance_maxElayer_cut.SetName("gr_acceptance_maxElayer_cut");
	gr_acceptance_maxBarLayer_cut.SetName("gr_acceptance_maxBarLayer_cut");
	gr_acceptance_BGOTrackContainment_cut.SetName("gr_acceptance_BGOTrackContainment_cut");
	gr_acceptance_BGO_fiducial_cut.SetName("gr_acceptance_BGO_fiducial_cut");
	gr_acceptance_BGO_fiducial_nBarLayer13_cut.SetName("gr_acceptance_BGO_fiducial_nBarLayer13_cut");
	gr_acceptance_BGO_fiducial_maxRms_cut.SetName("gr_acceptance_BGO_fiducial_maxRms_cut");
	gr_acceptance_BGO_fiducial_track_selection_cut.SetName("gr_acceptance_BGO_fiducial_track_selection_cut");
	gr_acceptance_BGO_fiducial_psd_stk_match_cut.SetName("gr_acceptance_BGO_fiducial_psd_stk_match_cut");
	gr_acceptance_BGO_fiducial_psd_charge_cut.SetName("gr_acceptance_BGO_fiducial_psd_charge_cut");
	gr_acceptance_BGO_fiducial_stk_charge_cut.SetName("gr_acceptance_BGO_fiducial_stk_charge_cut");
	gr_acceptance_BGO_fiducial_xtrl_cut.SetName("gr_acceptance_BGO_fiducial_xtrl_cut");
	gr_acceptance_all_cut.SetName("gr_acceptance_all_cut");

	gr_acceptance_geometric_factor.SetTitle("Geometric Factor");
	gr_acceptance_gometric_cut.SetTitle("Acceptance - geometric cut");
	gr_acceptance_maxElayer_cut.SetTitle("Acceptance - maxElateral cut");
	gr_acceptance_maxBarLayer_cut.SetTitle("Acceptance - maxBarLayer cut");
	gr_acceptance_BGOTrackContainment_cut.SetTitle("Acceptance - BGOTrackContainment cut");
	gr_acceptance_BGO_fiducial_cut.SetTitle("Acceptance - BGO fiducial volume cut");
	gr_acceptance_BGO_fiducial_nBarLayer13_cut.SetTitle("Acceptance - BGO fiducial volume + nBarLayer13 cut");
	gr_acceptance_BGO_fiducial_maxRms_cut.SetTitle("Acceptance - BGO fiducial volume + maxRms cut");
	gr_acceptance_BGO_fiducial_track_selection_cut.SetTitle("Acceptance - BGO fiducial volume + track selection cut");
	gr_acceptance_BGO_fiducial_psd_stk_match_cut.SetTitle("Acceptance - BGO fiducial volume + PSD-STK match cut");
	gr_acceptance_BGO_fiducial_psd_charge_cut.SetTitle("Acceptance - BGO fiducial volume + STK charge cut");
	gr_acceptance_BGO_fiducial_stk_charge_cut.SetTitle("Acceptance - BGO fiducial volume + STK charge cut");
	gr_acceptance_BGO_fiducial_xtrl_cut.SetTitle("Acceptance - BGO fiducial volume + XTRL cut");
	gr_acceptance_all_cut.SetTitle("Acceptance - all cut");

	// Write histos to file
	// Cut histos
	h_geo_factor.Write();
	h_incoming.Write();
	h_trigger.Write();
	h_geometric_cut.Write();
	h_maxElayer_cut.Write();
	h_maxBarLayer_cut.Write();
	h_BGOTrackContainment_cut.Write();
	h_BGO_fiducial_cut.Write();
	h_all_cut.Write();
	// Cuts && Geometric Cut
	h_geometric_maxElayer_cut.Write();
	h_geometric_maxBarLayer_cut.Write();
	h_geometric_BGOTrackContainment_cut.Write();
	h_geometric_BGO_fiducial_cut.Write();
	h_geometric_all_cut.Write();
	// && BGO fiducial volume cut
	h_BGOfiducial_nBarLayer13_cut.Write();
	h_BGOfiducial_maxRms_cut.Write();
	h_BGOfiducial_track_selection_cut.Write();
	h_BGOfiducial_psd_stk_match_cut.Write();
	h_BGOfiducial_psd_charge_cut.Write();
	h_BGOfiducial_stk_charge_cut.Write();
	h_BGOfiducial_xtrl_cut.Write();
	h_BGOfiducial_all_cut.Write();

	auto nw_histos_dir = outFile.mkdir("nw_histos");
	nw_histos_dir->cd();

	// Cut histos
	h_geo_factor_nw.Write();
	h_incoming_nw.Write();
	h_trigger_nw.Write();
	h_geometric_cut_nw.Write();
	h_maxElayer_cut_nw.Write();
	h_maxBarLayer_cut_nw.Write();
	h_BGOTrackContainment_cut_nw.Write();
	h_BGO_fiducial_cut_nw.Write();
	h_all_cut_nw.Write();
	// Cuts && Geometric Cut
	h_geometric_maxElayer_cut_nw.Write();
	h_geometric_maxBarLayer_cut_nw.Write();
	h_geometric_BGOTrackContainment_cut_nw.Write();
	h_geometric_BGO_fiducial_cut_nw.Write();
	h_geometric_all_cut_nw.Write();
	// && BGO fiducial volume cut
	h_BGOfiducial_nBarLayer13_cut_nw.Write();
	h_BGOfiducial_maxRms_cut_nw.Write();
	h_BGOfiducial_track_selection_cut_nw.Write();
	h_BGOfiducial_psd_stk_match_cut_nw.Write();
	h_BGOfiducial_psd_charge_cut_nw.Write();
	h_BGOfiducial_stk_charge_cut_nw.Write();
	h_BGOfiducial_xtrl_cut_nw.Write();
	h_BGOfiducial_all_cut_nw.Write();

	// Create output acceptance_histo dir in the output TFile
	auto acceptanceHistoDir = outFile.mkdir("Acceptance_histos");
	acceptanceHistoDir->cd();

	h_acceptance_geometric_factor->Write();
	h_acceptance_gometric_cut->Write();
	h_acceptance_maxElayer_cut->Write();
	h_acceptance_maxBarLayer_cut->Write();
	h_acceptance_BGOTrackContainment_cut->Write();
	h_acceptance_BGO_fiducial_cut->Write();
	h_acceptance_BGO_fiducial_nBarLayer13_cut->Write();
	h_acceptance_BGO_fiducial_maxRms_cut->Write();
	h_acceptance_BGO_fiducial_track_selection_cut->Write();
	h_acceptance_BGO_fiducial_psd_stk_match_cut->Write();
	h_acceptance_BGO_fiducial_psd_charge_cut->Write();
	h_acceptance_BGO_fiducial_stk_charge_cut->Write();
	h_acceptance_BGO_fiducial_xtrl_cut->Write();
	h_acceptance_all_cut->Write();

	// Create output acceptance dir in the output TFile
	auto acceptanceDir = outFile.mkdir("Acceptance");
	acceptanceDir->cd();

	// Write final TGraphs
	gr_acceptance_gometric_cut.Write();
	gr_acceptance_maxElayer_cut.Write();
	gr_acceptance_maxBarLayer_cut.Write();
	gr_acceptance_BGOTrackContainment_cut.Write();
	gr_acceptance_BGO_fiducial_cut.Write();
	gr_acceptance_BGO_fiducial_nBarLayer13_cut.Write();
	gr_acceptance_BGO_fiducial_maxRms_cut.Write();
	gr_acceptance_BGO_fiducial_track_selection_cut.Write();
	gr_acceptance_BGO_fiducial_psd_stk_match_cut.Write();
	gr_acceptance_BGO_fiducial_psd_charge_cut.Write();
	gr_acceptance_BGO_fiducial_stk_charge_cut.Write();
	gr_acceptance_BGO_fiducial_xtrl_cut.Write();
	gr_acceptance_all_cut.Write();

	auto geoFactor = outFile.mkdir("GeometricFactor");
	geoFactor->cd();

	gr_acceptance_geometric_factor.Write();
	
	// Create output ratio dir in the output TFile
	auto ratioDir = outFile.mkdir("Efficiency");

	// Create trigger folder
	auto trigger_dir = ratioDir->mkdir("Trigger");
	trigger_dir->cd();

	// Define TEfficiency pointers
	std::shared_ptr<TEfficiency> trigger_efficiency;
	std::shared_ptr<TEfficiency> tr_eff_gometric_cut;
	std::shared_ptr<TEfficiency> tr_eff_maxElayer_cut;
	std::shared_ptr<TEfficiency> tr_eff_maxBarLayer_cut;
	std::shared_ptr<TEfficiency> tr_eff_BGOTrackContainment_cut;
	std::shared_ptr<TEfficiency> tr_eff_BGO_fiducial_cut;
	std::shared_ptr<TEfficiency> tr_eff_all_cut;

	if (TEfficiency::CheckConsistency(h_geometric_cut_nw, h_geo_factor_nw))
		trigger_efficiency = std::make_shared<TEfficiency>(h_geometric_cut_nw, h_geo_factor_nw);

	if (TEfficiency::CheckConsistency(h_geometric_cut_nw, h_trigger_nw))
		tr_eff_gometric_cut = std::make_shared<TEfficiency>(h_geometric_cut_nw, h_trigger_nw);

	if (TEfficiency::CheckConsistency(h_maxElayer_cut_nw, h_trigger_nw))
		tr_eff_maxElayer_cut = std::make_shared<TEfficiency>(h_maxElayer_cut_nw, h_trigger_nw);

	if (TEfficiency::CheckConsistency(h_maxBarLayer_cut_nw, h_trigger_nw))
		tr_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(h_maxBarLayer_cut_nw, h_trigger_nw);

	if (TEfficiency::CheckConsistency(h_BGOTrackContainment_cut_nw, h_trigger_nw))
		tr_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(h_BGOTrackContainment_cut_nw, h_trigger_nw);

	if (TEfficiency::CheckConsistency(h_BGO_fiducial_cut_nw, h_trigger_nw))
		tr_eff_BGO_fiducial_cut = std::make_shared<TEfficiency>(h_BGO_fiducial_cut_nw, h_trigger_nw);

	if (TEfficiency::CheckConsistency(h_all_cut_nw, h_trigger_nw))
		tr_eff_all_cut = std::make_shared<TEfficiency>(h_all_cut_nw, h_trigger_nw);

	// Set uniform statistic option
	trigger_efficiency->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_gometric_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_maxElayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_maxBarLayer_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_BGOTrackContainment_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_BGO_fiducial_cut->SetStatisticOption(TEfficiency::kBUniform);
	tr_eff_all_cut->SetStatisticOption(TEfficiency::kBUniform);

	trigger_efficiency->SetName("trigger_efficiency");
	tr_eff_gometric_cut->SetName("tr_eff_gometric_cut");
	tr_eff_maxElayer_cut->SetName("tr_eff_maxElayer_cut");
	tr_eff_maxBarLayer_cut->SetName("tr_eff_maxBarLayer_cut");
	tr_eff_BGOTrackContainment_cut->SetName("tr_eff_BGOTrackContainment_cut");
	tr_eff_BGO_fiducial_cut->SetName("tr_eff_BGO_fiducial_cut");
	tr_eff_all_cut->SetName("tr_eff_all_cut");

	trigger_efficiency->SetTitle("Trigger efficiency");
	tr_eff_gometric_cut->SetTitle("Gometric cut efficiency");
	tr_eff_maxElayer_cut->SetTitle("maxElayer cut efficiency");
	tr_eff_maxBarLayer_cut->SetTitle("maxBarLayer cut efficiency");
	tr_eff_BGOTrackContainment_cut->SetTitle("BGOTrackContainment cut efficiency");
	tr_eff_all_cut->SetTitle("all cut efficiency");

	// Write histos to disk
	trigger_efficiency->Write();
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

	if (TEfficiency::CheckConsistency(h_geometric_maxElayer_cut_nw, h_geometric_cut_nw))
		geo_eff_maxElayer_cut = std::make_shared<TEfficiency>(h_geometric_maxElayer_cut_nw, h_geometric_cut_nw);

	if (TEfficiency::CheckConsistency(h_geometric_maxBarLayer_cut_nw, h_geometric_cut_nw))
		geo_eff_maxBarLayer_cut = std::make_shared<TEfficiency>(h_geometric_maxBarLayer_cut_nw, h_geometric_cut_nw);

	if (TEfficiency::CheckConsistency(h_geometric_BGOTrackContainment_cut_nw, h_geometric_cut_nw))
		geo_eff_BGOTrackContainment_cut = std::make_shared<TEfficiency>(h_geometric_BGOTrackContainment_cut_nw, h_geometric_cut_nw);

	if (TEfficiency::CheckConsistency(h_geometric_BGO_fiducial_cut_nw, h_geometric_cut_nw))
		geo_eff_BGO_fiducial = std::make_shared<TEfficiency>(h_geometric_BGO_fiducial_cut_nw, h_geometric_cut_nw);

	if (TEfficiency::CheckConsistency(h_geometric_all_cut_nw, h_geometric_cut_nw))
		geo_eff_all_cut = std::make_shared<TEfficiency>(h_geometric_all_cut_nw, h_geometric_cut_nw);

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


	if (TEfficiency::CheckConsistency(h_BGOfiducial_nBarLayer13_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_nBarLayer13_cut = std::make_shared<TEfficiency>(h_BGOfiducial_nBarLayer13_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_maxRms_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_l13_maxRms_cut = std::make_shared<TEfficiency>(h_BGOfiducial_maxRms_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_track_selection_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_l13_rms_track_selection_cut = std::make_shared<TEfficiency>(h_BGOfiducial_track_selection_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_stk_match_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_l13_rms_ts_psd_stk_match_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_stk_match_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_charge_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_l13_rms_ts_psdstk_psd_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_charge_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_stk_charge_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_l13_rms_ts_psdstk_pc_stk_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_stk_charge_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_xtrl_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_l13_rms_ts_psdstk_pc_sc_xtrl_cut = std::make_shared<TEfficiency>(h_BGOfiducial_xtrl_cut_nw, h_BGO_fiducial_cut_nw);

	if (TEfficiency::CheckConsistency(h_BGOfiducial_all_cut_nw, h_BGO_fiducial_cut_nw))
		BGOfiducial_eff_all_cut = std::make_shared<TEfficiency>(h_BGOfiducial_all_cut_nw, h_BGO_fiducial_cut_nw);

	if (active_cuts.maxRms && active_cuts.nBarLayer13)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_maxRms_cut_nw, h_BGOfiducial_nBarLayer13_cut_nw))
			BGOfiducial_eff_maxRms_cut = std::make_shared<TEfficiency>(h_BGOfiducial_maxRms_cut_nw, h_BGOfiducial_nBarLayer13_cut_nw);

	if (active_cuts.track_selection && active_cuts.maxRms)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_track_selection_cut_nw, h_BGOfiducial_maxRms_cut_nw))
			BGOfiducial_eff_track_selection_cut = std::make_shared<TEfficiency>(h_BGOfiducial_track_selection_cut_nw, h_BGOfiducial_maxRms_cut_nw);

	if (active_cuts.psd_stk_match && active_cuts.track_selection)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_stk_match_cut_nw, h_BGOfiducial_track_selection_cut_nw))
			BGOfiducial_eff_psd_stk_match_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_stk_match_cut_nw, h_BGOfiducial_track_selection_cut_nw);

	if (active_cuts.psd_charge && active_cuts.psd_stk_match)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_psd_charge_cut_nw, h_BGOfiducial_psd_stk_match_cut_nw))
			BGOfiducial_eff_psd_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_psd_charge_cut_nw, h_BGOfiducial_psd_stk_match_cut_nw);

	if (active_cuts.stk_charge && active_cuts.psd_charge)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_stk_charge_cut_nw, h_BGOfiducial_psd_stk_match_cut_nw))
			BGOfiducial_eff_stk_charge_cut = std::make_shared<TEfficiency>(h_BGOfiducial_stk_charge_cut_nw, h_BGOfiducial_psd_stk_match_cut_nw);

	if (active_cuts.xtrl && active_cuts.stk_charge)
		if (TEfficiency::CheckConsistency(h_BGOfiducial_xtrl_cut_nw, h_BGOfiducial_stk_charge_cut_nw))
			BGOfiducial_eff_xtrl_cut = std::make_shared<TEfficiency>(h_BGOfiducial_xtrl_cut_nw, h_BGOfiducial_stk_charge_cut_nw);

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
	
	// Create output analysis dir in the output TFile
	auto preGeo_analysisDir = outFile.mkdir("Analysis_preGeoCut");
	preGeo_analysisDir->cd();

	h_preGeo_BGOrec_topX_vs_realX.Write();
	h_preGeo_BGOrec_topY_vs_realY.Write();
	h_preGeo_real_slopeX.Write();
	h_preGeo_real_slopeY.Write();
	h_preGeo_BGOrec_slopeX.Write();
	h_preGeo_BGOrec_slopeY.Write();
	h_preGeo_real_interceptX.Write();
	h_preGeo_real_interceptY.Write();
	h_preGeo_BGOrec_interceptX.Write();
	h_preGeo_BGOrec_interceptY.Write();
	h_preGeo_real_topMap.Write();
	h_preGeo_BGOreco_topMap.Write();
	h_preGeo_real_bottomMap.Write();
	h_preGeo_BGOreco_bottomMap.Write();

	h_noBGOenergy_real_topMap.Write();

	auto geo_analysisDir = outFile.mkdir("Analysis_GeoCut");
	geo_analysisDir->cd();

	h_geo_BGOrec_topX_vs_realX.Write();
	h_geo_BGOrec_topY_vs_realY.Write();
	h_geo_real_slopeX.Write();
	h_geo_real_slopeY.Write();
	h_geo_BGOrec_slopeX.Write();
	h_geo_BGOrec_slopeY.Write();
	h_geo_real_interceptX.Write();
	h_geo_real_interceptY.Write();
	h_geo_BGOrec_interceptX.Write();
	h_geo_BGOrec_interceptY.Write();
	h_geo_real_topMap.Write();
	h_geo_BGOreco_topMap.Write();
	h_geo_real_bottomMap.Write();
	h_geo_BGOreco_bottomMap.Write();

	auto BGOdir = outFile.mkdir("BGO_Energy");
	BGOdir->cd();

	h_BGOrec_energy.Write();
	h_simu_energy.Write();
	h_energy_diff.Write();
	h_energy_diff2D.Write();
	h_energy_unfold.Write();
	h_layer_max_energy_ratio.Write();

	sumRms_cosine_20_100.Write();
	sumRms_cosine_100_250.Write();
	sumRms_cosine_250_500.Write();
	sumRms_cosine_500_1000.Write();
	sumRms_cosine_1000_3000.Write();
	sumRms_cosine_3000_10000.Write();

	for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
		h_layer_energy_ratio[lIdx].Write();

	for (auto it = sumRms_cosine.begin(); it != sumRms_cosine.end(); ++it)
		(*it)->Write();

	auto XTRLdir = outFile.mkdir("xtrl");
	XTRLdir->cd();

	h_xtrl_energy_int.Write();
	h_xtrl.Write();
	e_discrimination.Write();
	e_discrimination_20_100.Sumw2();
	e_discrimination_100_250.Sumw2();
	e_discrimination_250_500.Sumw2();
	e_discrimination_500_1000.Sumw2();
	e_discrimination_1000_3000.Sumw2();
	e_discrimination_3000_10000.Sumw2();

	e_discrimination_last.Write();
	e_discrimination_last_20_100.Sumw2();
	e_discrimination_last_100_250.Sumw2();
	e_discrimination_last_250_500.Sumw2();
	e_discrimination_last_500_1000.Sumw2();
	e_discrimination_last_1000_3000.Sumw2();
	e_discrimination_last_3000_10000.Sumw2();

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
	if (ancillary_cuts.compute_proton_background)
	{
		auto ancillaryDir = outFile.mkdir("mc_ancillary");
		ancillaryDir->cd();

		h_background_under_xtrl_cut.Write();
		h_background_over_xtrl_cut.Write();

		// Create proton background ratio
		auto proton_background_ratio = static_cast<TH1D *>(h_background_under_xtrl_cut.Clone("proton_background_ratio"));
		proton_background_ratio->SetTitle("Proton background ratio");
		proton_background_ratio->Divide(&h_background_over_xtrl_cut);

		proton_background_ratio->Write();
	}

	outFile.cd();
}