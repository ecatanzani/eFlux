#include "utils.h"
#include "config.h"
#include "energy.h"
#include "Data/data.h"
#include "tmpstruct.h"
#include "Tuple/data_tuple.h"
#include "aggregate_events.h"
#include "Data/data_tmpstruct.h"
#include "Efficiency/efficiency.h"

#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpPsdContainer.h"
#include "Dmp/DmpNudContainer.h"
#include "Dmp/DmpFilterContainer.h"

#include "TClonesArray.h"

#include "DmpCore.h"
#include "DmpIOSvc.h"
#include "DmpStkTrack.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "DmpFilterOrbit.h"
#include "DmpEvtAttitude.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtSimuPrimaries.h"

#include <vector>
#include <memory>

void rawDataLoop(
	const std::string inputPath,
	TFile &outFile,
	const bool verbose,
	const std::string wd)
	{
		bool _MC = false;

		std::shared_ptr<TChain> dmpch;
		event_collector collector(
			inputPath,
			verbose,
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

		// Register attitude container
		std::shared_ptr<DmpEvtAttitude> attitude = std::make_shared<DmpEvtAttitude>();
		dmpch->SetBranchAddress("EvtAttitudeContainer", &attitude);

		// Orbit filter
		// Set gIOSvc
		gIOSvc->Set("OutData/NoOutput", "True");
		gIOSvc->Initialize();
		// Create orbit filter
		std::shared_ptr<DmpFilterOrbit> pFilter = std::make_shared<DmpFilterOrbit>("EventHeader");
		// Activate orbit filter
		pFilter->ActiveMe(); // Call this function to calculate SAA through House Keeping Data

		// Load config class
		config data_config(wd, _MC);
		// Load energy class
		energy evt_energy;
		// Load filter class
		DmpFilterContainer filter;
		// Load BGO, PSD and NUD high level classes
		DmpBgoContainer bgoVault;
		DmpStkContainer stkVault;
		DmpPsdContainer psdVault;
		DmpNudContainer nudVault;
		// Load efficiency class
		efficiency cuts_efficiency;
		// Load output file
		outFile.cd();
		// Create DATA tuple objects
		std::unique_ptr<data_tuple> tuple = std::make_unique<data_tuple>(data_config.GetActiveCuts());
		
		auto nevents = dmpch->GetEntries();
		int kStep = 10;
		// Update event time
		filter.UpdateEvtTime(dmpch, evt_header);
		if (verbose)
		{
			data_config.PrintActiveFilters();
			filter.PrintDataInfo(dmpch, _MC);
			std::cout << "Analysing RawData...\n\n";
		}

		for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
		{
			// Read tree event
			dmpch->GetEvent(evIdx);

			// Reset classes
			if (evIdx) {
				// Reset filter event flags
				filter.Reset();
				// Reset efficiency filter event flags
				cuts_efficiency.Reset();
				// Reset tuple
				tuple->Reset();
				// Reset energy class
				evt_energy.Reset();
				// Reset BGO class
				bgoVault.Reset();
				// Reset STK class
				stkVault.Reset();
				// Reset PSD class
				psdVault.Reset();
				// Reset NUD class
				nudVault.Reset();
			}
			
			// Update event counter
			filter.UpdateEvtCounter();
			// Status printout
			if (verbose)
				UpdateProcessStatus(evIdx, kStep, nevents);

			// Load energy class
			evt_energy.SetRawEnergy(bgorec->GetTotalEnergy());
			evt_energy.SetCorrEnergy(bgorec->GetElectronEcor());
			// Check SAA
			filter.SAACheck(evt_header, pFilter);
			// Check particle energy
			filter.EnergyCheck(
				data_config.GetCutsConfigValues(),
				evt_energy.GetCorrEnergy(),
				data_config.GetMinEnergyRange(),
				data_config.GetMaxEnergyRange());

			// Load sub-detectors high level class
			// Load STK class
			stkVault.scanSTKHits(stkclusters);
			// Load BGO class
			bgoVault.scanBGOHits(bgohits, bgorec, evt_energy.GetRawEnergy(), data_config.GetBGOLayerMinEnergy());
			// Load PSD class
			psdVault.scanPSDHits(psdhits, data_config.GetPSDBarMinEnergy());
			// Load NUD class
			nudVault.scanNudHits(nudraw);


			// Check current event (Trigger and BGO reco)
			if (filter.CheckIncomingEvent(evt_header, bgoVault.GetBGOslope(), bgoVault.GetBGOintercept()))
			{
				// Check BGO geometry after trigger
				filter.CheckGeometry(
					std::shared_ptr<DmpEvtSimuPrimaries>(nullptr),
					bgoVault.FastBGOslope(bgorec),
					bgoVault.FastBGOintercept(bgorec));
				// Filter event
				filter.Pipeline(
					bgorec,
					bgohits,
					data_config.GetCutsConfigValues(),
					evt_energy.GetRawEnergy(),
					evt_energy.GetCorrEnergy(),
					bgoVault,
					psdVault,
					stkclusters,
					stktracks,
					data_config.GetActiveCuts());
			}
			
			// Efficiency pipeline
			if (filter.CheckIncomingEventNoTrigger())
				cuts_efficiency.Pipeline(
					bgorec,
					bgohits,
					data_config.GetCutsConfigValues(),
					evt_energy.GetRawEnergy(),
					evt_energy.GetCorrEnergy(),
					bgoVault,
					psdVault,
					stkclusters,
					stktracks,
					data_config.GetActiveCuts(),
					filter.GetTrigger());

			// Fill output structures
			tuple->Fill(
				fillFilterTmpStruct(filter, cuts_efficiency),
				stkVault.GetNPlaneClusters(),
				fillBGOTmpStruct(bgoVault),
				fillDataEnergyTmpStruct(evt_energy),
				attitude,
				evt_header,
				fillNUDTmpStruct(nudVault));
		}

		if (verbose)
		{
			std::cout << "\n\n ****** \n\n";
			std::cout << "Triggered events: " << filter.GetStatiEvtTrigger() << std::endl;
			std::cout << "Triggered events in energy range: " << filter.GetStatEvtInRange() << std::endl;
			std::cout << "Triggered events out of energy range: " << filter.GetStatEvtOutRange() << std::endl;
			std::cout << "Events in SAA: " << filter.GetStatEvtSAA() << std::endl;
			std::cout << "\n\n**** Filter result ****\n";
			std::cout << "***********************\n\n";
			std::cout << "Particles surviving the selection cuts: " << filter.GetStatEvtSelection() << "\t | " << ((double)filter.GetStatEvtSelection() / filter.GetStatiEvtTrigger()) * 100 << "%";
			std::cout << "\n\n ***********************\n\n";
		}
		
		tuple->Write(outFile);
	}