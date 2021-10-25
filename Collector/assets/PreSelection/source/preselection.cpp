#include "bgo.h"
#include "psd.h"
#include "chain.h"
#include "histos.h"
#include "preselection.h"
#include "energy_config.h"

#include <iostream>
#include <memory>

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtPsdHits.h"
#include "DmpSimuTrajectory.h"
#include "DmpEvtSimuPrimaries.h"

#include "TChain.h"
#include "TClonesArray.h"

void preselection(const in_pars &input_pars) {
    
    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose);
    std::shared_ptr<energy_config> econfig = std::make_shared<energy_config>(input_pars.config_wd);

    // Register Header container
    std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
    evtch->SetBranchAddress("EventHeader", &evt_header);

    // Register BGO container
    std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
    evtch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC container
    std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
    evtch->SetBranchAddress("DmpEvtBgoRec", &bgorec);
    
    // Register STK container
	std::shared_ptr<TClonesArray> stkclusters = std::make_shared<TClonesArray>("DmpStkSiCluster");
	evtch->SetBranchAddress("StkClusterCollection", &stkclusters);

    // Register SimuPrimaries container
	std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries;
    // Register SimuTrajectory container
	std::shared_ptr<TClonesArray> simu_trajectories;
    if (input_pars.mc) {
        simu_primaries = std::make_shared<DmpEvtSimuPrimaries>();
        simu_trajectories = std::make_shared<TClonesArray>("DmpSimuTrajectory");    
        evtch->SetBranchAddress("DmpEvtSimuPrimaries", &simu_primaries);
        evtch->SetBranchAddress("DmpTruthTrajectoriesCollection", &simu_trajectories);
    }

	// Check if STK tracks collection exists
	bool fStkKalmanTracksFound = false;
	for (int brIdx = 0; brIdx < evtch->GetListOfBranches()->GetEntries(); ++brIdx)
		if (strcmp(evtch->GetListOfBranches()->At(brIdx)->GetName(), "StkKalmanTracks")) {
			fStkKalmanTracksFound = true;
			break;
		}

	// Register STK tracks collection
	std::shared_ptr<TClonesArray> stktracks = std::make_shared<TClonesArray>("DmpStkTrack", 200);
	if (fStkKalmanTracksFound)
		evtch->SetBranchAddress("StkKalmanTracks", &stktracks);

	// Register PSD container
	std::shared_ptr<DmpEvtPsdHits> psdhits = std::make_shared<DmpEvtPsdHits>();
	evtch->SetBranchAddress("DmpPsdHits", &psdhits);

    auto min_evt_energy = econfig->GetMinEvtEnergy();
    auto max_evt_energy = econfig->GetMaxEvtEnergy();

    double evt_corr_energy {0};
    double evt_energy {0};
    double evt_corr_energy_gev {0};
    double evt_energy_gev {0};
    double gev {0.001};
    unsigned int perc {10};

    auto nevents {evtch->GetEntries()};

    std::shared_ptr<histos> ps_histos = std::make_shared<histos>(econfig, input_pars.mc);

    if (input_pars.verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx) {

        evtch->GetEvent(evIdx);
        if (input_pars.verbose && !((evIdx+1)%(nevents*perc/100))) {
            std::cout << "\nProcessing ... [ " << perc << "% ]";
            perc += 10;
        }
    
        evt_energy = bgorec->GetTotalEnergy();
        evt_corr_energy = bgorec->GetElectronEcor();
        evt_energy_gev = evt_energy*gev;
        evt_corr_energy_gev = evt_corr_energy*gev;

        if (evt_corr_energy_gev>=min_evt_energy && evt_corr_energy_gev<=max_evt_energy) {
            bgofiducial_distributions(
                bgohits, 
                bgorec, 
                evt_header, 
                simu_primaries,
                evt_corr_energy_gev, 
                ps_histos);
            psd_stk_distributions(
                bgohits, 
                bgorec, 
                evt_header, 
                stkclusters, 
                stktracks, 
                psdhits, 
                evt_energy, 
                evt_corr_energy, 
                evt_energy_gev, 
                evt_corr_energy_gev, 
                ps_histos);

        }

    }

    ps_histos->Write(input_pars.output_wd, input_pars.verbose);
}