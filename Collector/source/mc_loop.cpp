#include "mc.h"
#include "utils.h"
#include "energy.h"
#include "tmpstruct.h"
#include "mc_tmpstruct.h"
#include "config.h"
#include "particle.h"
#include "mc_tuple.h"
#include "DmpStkContainer.h"
#include "DmpBgoContainer.h"
#include "DmpPsdContainer.h"
#include "DmpNudContainer.h"
#include "aggregate_events.h"
#include "DmpFilterContainer.h"
#include "DAMPE_geo_structure.h"

#include "TClonesArray.h"

#include "DmpStkHit.hh"
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
	// Load energy class
	energy evt_energy;
	// Load filter class
	DmpFilterContainer filter;
	// Load output file
	outFile.cd();
	// Create MC tuple objects
	std::unique_ptr<mc_tuple> simu_tuple = std::make_unique<mc_tuple>(mc_config.GetActiveCuts());
	
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
		DmpStkContainer stkVault;
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
			evt_energy.SetSimuEnergyWeight(pow(evt_energy.GetSimuEnergy() * _GeV, -2));
			evt_energy.SetCorrEnergyWeight(pow(evt_energy.GetCorrEnergy() * _GeV, -2));
		}
		else if (simu_particle.IsProton())
		{
			evt_energy.SetSimuEnergyWeight(pow(evt_energy.GetSimuEnergy() * _GeV, -1.7));
			evt_energy.SetCorrEnergyWeight(pow(evt_energy.GetCorrEnergy() * _GeV, -1.7));
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
			// Load STK class
			stkVault.scanSTKHits(stkclusters);
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
				bgorec,
				bgohits,
				mc_config.GetCutsConfigValues(),
				evt_energy.GetRawEnergy(),
				evt_energy.GetCorrEnergy(),
				bgoVault,
				psdVault,
				stkclusters,
				stktracks,
				mc_config.GetActiveCuts());
		}
		
		// Fill output structures
		simu_tuple->Fill(
			fillFilterTmpStruct(filter),
			stkVault.GetNPlaneClusters(),
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