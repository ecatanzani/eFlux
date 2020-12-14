#include "mc.h"
#include "utils.h"
#include "energy.h"
#include "config.h"
#include "particle.h"
#include "mc_tuple.h"
#include "DmpBgoContainer.h"
#include "DmpPsdContainer.h"
#include "DmpNudContainer.h"
#include "aggregate_events.h"
#include "DmpFilterContainer.h"
#include "DAMPE_geo_structure.h"

#include "TClonesArray.h"

#include "DmpStkTrack.h"
#include "DmpEvtNudRaw.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "DmpStkSiCluster.h"
#include "DmpSimuTrajectory.h"
#include "DmpEvtSimuPrimaries.h"

#include <vector>
#include <memory>

void mcLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool _VERBOSE,
	std::shared_ptr<ofstream> evlogger,
	const std::string wd)
{
	bool _MC = true;
	double _GeV = 0.001;

	std::shared_ptr<TChain> dmpch;
	event_collector collector(
		inputPath,
		_VERBOSE,
		_MC);
	if (collector.GetChainStatus())
		dmpch = collector.GetChain();
	else
	{
		std::cerr << "\n\nERROR: Corrupted TChain object...";
		exit(100);
	}

	// Register Header container
	std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
	dmpch->SetBranchAddress("EventHeader", &evt_header);

	// Register SimuPrimaries container
	std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::make_shared<DmpEvtSimuPrimaries>();
	dmpch->SetBranchAddress("DmpEvtSimuPrimaries", &simu_primaries);

	// Register SimuTrajectory container
	std::shared_ptr<TClonesArray> simu_trajectories = std::make_shared<TClonesArray>("DmpSimuTrajectory");
	dmpch->SetBranchAddress("DmpTruthTrajectoriesCollection", &simu_trajectories);

	// Register BGO container
	std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
	dmpch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

	// Register BGO REC container
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

	// Register NUD container
	std::shared_ptr<DmpEvtNudRaw> nudraw = std::make_shared<DmpEvtNudRaw>();
	dmpch->SetBranchAddress("DmpEvtNudRaw", &nudraw);

	// Load particle ID class
	mc_particle simu_particle;
	// Load config class
	config mc_config(wd, _MC);
	// Compute energy binning
	auto logEBins = mc_config.GetEnergyBinning();
	// Create MC tuple objects
	std::unique_ptr<mc_tuple> simu_tuple = std::make_unique<mc_tuple>(mc_config.GetActiveCuts());
	// Load energy class
	energy evt_energy;
	// Load filter class
	DmpFilterContainer filter;
	// Load output file
	outFile.cd();

	auto nevents = dmpch->GetEntries();
	int kStep = 10;
	if (_VERBOSE)
	{
		mc_config.PrintActiveFilters();
		filter.PrintDataInfo(dmpch, _MC);
		std::cout << "Analysing MC...\n\n";
	}
	
	for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
	{
		// Read tree event
		dmpch->GetEvent(evIdx);
		// Reset filter event flags
		filter.Reset();
		// Reset tuple
		simu_tuple->Reset();
		// Build BGO, PSD and NUD vault objects
		DmpBgoContainer bgoVault;
		DmpPsdContainer psdVault;
		DmpNudContainer nudVault;
		// Load particle ID
		simu_particle.Load(simu_primaries->pvpart_pdg);
		// Update event counter
		filter.UpdateEvtCounter();
		// Status printout
		if (_VERBOSE)
			UpdateProcessStatus(evIdx, kStep, nevents);
		// Reset energy class
		if (evIdx)
			evt_energy.Reset();
		// Load energy class
		evt_energy.SetEnergies(
			simu_primaries->pvpart_ekin,
			bgorec->GetTotalEnergy(),
			bgorec->GetElectronEcor());
		// Compute energy weights for MC event
		if (simu_particle.IsElectron())
		{
			evt_energy.SetSimuEnergyWeight(pow(evt_energy.GetSimuEnergy() * _GeV, -2)*pow(logEBins[0], 2));
			evt_energy.SetCorrEnergyWeight(pow(evt_energy.GetCorrEnergy() * _GeV, -2)*pow(logEBins[0], 2));
		}
		else if (simu_particle.IsProton())
		{
			evt_energy.SetSimuEnergyWeight(pow(evt_energy.GetSimuEnergy() * _GeV, -1.7)*pow(logEBins[0], 2));
			evt_energy.SetCorrEnergyWeight(pow(evt_energy.GetCorrEnergy() * _GeV, -1.7)*pow(logEBins[0], 2));
		}
		// Check particle energy
		filter.EnergyCheck(
			mc_config.GetCutsConfigValues(),
			evt_energy.GetCorrEnergy(),
			mc_config.GetMinEnergyRange(),
			mc_config.GetMaxEnergyRange());
		// Check BGO geometry before trigger
		filter.CheckGeometry(simu_primaries);
		// Check current event (Trigger and BGO reco)
		auto _good_evt =
			filter.CheckIncomingEvent(
				evt_header,
				bgoVault.FastBGOslope(bgorec),
				bgoVault.FastBGOintercept(bgorec),
				simu_primaries);
		if (_good_evt)
		{
			// Check BGO geometry after trigger
			filter.CheckGeometry(simu_primaries);
			// Load BGO class
			bgoVault.scanBGOHits(
				bgohits,
				bgorec,
				evt_energy.GetRawEnergy(),
				mc_config.GetBGOLayerMinEnergy());
			// Load PSD class
			psdVault.scanPSDHits(
				psdhits,
				mc_config.GetPSDBarMinEnergy());
			// Load NUD class
			nudVault.scanNudHits(nudraw);
			// Filter event
			filter.Pipeline(
				evIdx,
				bgorec,
				bgohits,
				mc_config.GetCutsConfigValues(),
				mc_config.GetLoggerCutsConfigValues(),
				evt_energy.GetRawEnergy(),
				evt_energy.GetCorrEnergy(),
				bgoVault,
				psdVault,
				stkclusters,
				stktracks,
				mc_config.GetActiveCuts(),
				mc_config.GetLoggerActiveCuts(),
				evlogger);
		}
		
		// Fill output structures
		simu_tuple->Fill(
			fillFilterTmpStruct(filter),
			fillBGOTmpStruct(bgoVault),
			fillSimuTmpStruct(simu_primaries, simu_trajectories),
			fillEnergyTmpStruct(evt_energy),
			fillNUDTmpStruct(nudVault));
	}

	if (_VERBOSE)
	{
		std::cout << "\n\n ****** \n\n";
		std::cout << "Generated events: " << filter.GetStatEvtCounter() << std::endl;
		std::cout << "Generated events in energy range: " << filter.GetStatEvtInRange() << std::endl;
		std::cout << "Generated events out of energy range: " << filter.GetStatEvtOutRange() << std::endl;
		std::cout << "Triggered events: " << filter.GetStatiEvtTrigger() << std::endl;
		std::cout << "\n\n**** Filter result ****\n";
		std::cout << "***********************\n\n";
		std::cout << "Particles surviving the selection cuts: " << filter.GetStatEvtSelection() << "\t | " << ((double)filter.GetStatEvtSelection() / filter.GetStatiEvtTrigger()) * 100 << "%";
		std::cout << "\n\n ***********************\n\n";
	}
	
	simu_tuple->Write(outFile);
}

std::shared_ptr<_tmp_simu> fillSimuTmpStruct(
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
	const std::shared_ptr<TClonesArray> simu_trajectories)
{
	std::shared_ptr<_tmp_simu> _simu_res = std::make_shared<_tmp_simu>();
	_simu_res->position = TVector3(simu_primaries->pv_x, simu_primaries->pv_y, simu_primaries->pv_z);
	_simu_res->momentum = TVector3(simu_primaries->pvpart_px, simu_primaries->pvpart_py, simu_primaries->pvpart_pz);
	_simu_res->radius = simu_primaries->pv_r;
	_simu_res->theta = simu_primaries->pv_theta;
	_simu_res->phi = simu_primaries->pv_phi;
	_simu_res->flux_w = simu_primaries->pv_fluxw;
	_simu_res->n_particle = simu_primaries->npvpart;
	_simu_res->cos_x = simu_primaries->pvpart_cosx;
	_simu_res->cos_y = simu_primaries->pvpart_cosy;
	_simu_res->cos_z = simu_primaries->pvpart_cosz;
	_simu_res->charge = simu_primaries->pvpart_q;
	_simu_res->zenith = simu_primaries->pvpart_zenith;
	_simu_res->azimuth = simu_primaries->pvpart_azimuth;
	_simu_res->w = simu_primaries->pvpart_weight;
	_simu_res->PDG = simu_primaries->pvpart_pdg;
	_simu_res->geocut = simu_primaries->pvpart_geocut;

	auto n_simu_trajectories = simu_trajectories->GetEntries();

	_simu_res->thruthtrajectory_x.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_y.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_z.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_energy.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_start_x.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_start_y.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_start_z.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_x.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_y.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_z.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_trackID.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_parentID.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_charge.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_PDG.resize(n_simu_trajectories);
	_simu_res->thruthtrajectory_stop_index.resize(n_simu_trajectories);

	for (int trIdx = 0; trIdx<n_simu_trajectories; ++trIdx)
	{
		auto simu_track = static_cast<DmpSimuTrajectory*>(simu_trajectories->ConstructedAt(trIdx));

		_simu_res->thruthtrajectory_x[trIdx] = simu_track->px;
		_simu_res->thruthtrajectory_y[trIdx] = simu_track->py;
		_simu_res->thruthtrajectory_z[trIdx] = simu_track->pz;
		_simu_res->thruthtrajectory_energy[trIdx] = simu_track->ekin;
		_simu_res->thruthtrajectory_start_x[trIdx] = simu_track->start_x;
		_simu_res->thruthtrajectory_start_y[trIdx] = simu_track->start_y;
		_simu_res->thruthtrajectory_start_z[trIdx] = simu_track->start_z;	
		_simu_res->thruthtrajectory_stop_x[trIdx] = simu_track->stop_x;
		_simu_res->thruthtrajectory_stop_y[trIdx] = simu_track->stop_y;
		_simu_res->thruthtrajectory_stop_z[trIdx] = simu_track->stop_z;
		_simu_res->thruthtrajectory_trackID[trIdx] = simu_track->trackID;
		_simu_res->thruthtrajectory_parentID[trIdx] = simu_track->parentID;
		_simu_res->thruthtrajectory_charge[trIdx] = simu_track->charge;
		_simu_res->thruthtrajectory_PDG[trIdx] = simu_track->pdg_id;
		_simu_res->thruthtrajectory_stop_index[trIdx] = simu_track->stop_index;
	}

	return _simu_res;
}

std::shared_ptr<_tmp_filter> fillFilterTmpStruct(DmpFilterContainer &filter)
{
	std::shared_ptr<_tmp_filter> _filter_res = std::make_shared<_tmp_filter>();	
	_filter_res->output = filter.GetFilterOutput();
	_filter_res->evt_best_track = filter.GetBestTrack();
	_filter_res->evt_psd_charge = filter.GetPSDCharge();
	_filter_res->evt_stk_charge = filter.GetSTKCharge();
	_filter_res->evt_bgo_classifier = filter.GetClassifiers();
	_filter_res->evt_trigger_info = filter.GetTrigger();

	return _filter_res;
}

std::shared_ptr<_tmp_bgo> fillBGOTmpStruct(DmpBgoContainer &bgoVault)
{
	std::shared_ptr<_tmp_bgo> _bgo_res = std::make_shared<_tmp_bgo>();

	_bgo_res->layer_energies = bgoVault.GetLayerEnergies();
	_bgo_res->layer_bar_energies = bgoVault.GetLayerBarEnergies();
	_bgo_res->slope = bgoVault.GetBGOslope();
	_bgo_res->intercept = bgoVault.GetBGOintercept();
	_bgo_res->trajectory2D = bgoVault.GetBGOTrajectory2D();
	_bgo_res->sumrms = bgoVault.GetSumRMS();
	_bgo_res->sumrms_layer = bgoVault.GetRmsLayer();
	_bgo_res->energy_fraction_layer = bgoVault.GetFracLayer();
	_bgo_res->energy_fraction_last_layer = bgoVault.GetSingleFracLayer(bgoVault.GetLastEnergyLayer());
	_bgo_res->energy_fraction_13th_layer = bgoVault.GetSingleFracLayer(13);
	_bgo_res->last_energy_layer = bgoVault.GetLastEnergyLayer();
	_bgo_res->hits = bgoVault.GetNhits();
	_bgo_res->energy_1mr = bgoVault.GetEnergy1MR();
	_bgo_res->energy_2mr = bgoVault.GetEnergy2MR();
	_bgo_res->energy_3mr = bgoVault.GetEnergy3MR();
	_bgo_res->energy_5mr = bgoVault.GetEnergy5MR();

	return _bgo_res;
}

std::shared_ptr<_tmp_energy> fillEnergyTmpStruct(energy &evt_energy)
{
	std::shared_ptr<_tmp_energy> _energy_res = std::make_shared<_tmp_energy>();

	_energy_res->simu = evt_energy.GetSimuEnergy();
	_energy_res->simu_w = evt_energy.GetSimuWeight();
	_energy_res->raw = evt_energy.GetRawEnergy();
	_energy_res->correct = evt_energy.GetCorrEnergy();
	_energy_res->correct_w = evt_energy.GetCorrWeight();

	return _energy_res;
}

std::shared_ptr<_tmp_nud> fillNUDTmpStruct(DmpNudContainer &nudVault)
{
	std::shared_ptr<_tmp_nud> _nud_res = std::make_shared<_tmp_nud>();
	
	_nud_res->adc = nudVault.GetADC();
	_nud_res->total_adc = nudVault.GetTotalADC();
	_nud_res->max_adc = nudVault.GetMaxADC();
	_nud_res->max_channel_ID = nudVault.GetMaxChannelID();

	return _nud_res;
}