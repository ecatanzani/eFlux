#include "focus.h"
#include "bgo.h"
#include "stk.h"
#include "psd.h"
#include "chain.h"
#include "histos.h"
#include "charge.h"
#include "config.h"
#include "preselection.h"
#include "energy_config.h"
#include "lateral_showering.h"

#include <iostream>
#include <memory>
#include <vector>
#include <tuple>

#include "Dmp/DmpGeoStruct.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtPsdHits.h"
#include "DmpSimuTrajectory.h"
#include "DmpEvtSimuPrimaries.h"

#include "DmpStkTrack.h"

#include "DmpIOSvc.h"
#include "DmpCore.h"
#include "DmpFilterOrbit.h"

#include "TChain.h"
#include "TClonesArray.h"

#include "TFile.h"
#include "TTree.h"

// Event display parameters
#define _BEST_TRACK         false
#define _BGO_EDGE           true
#define _OUTSIDE_BGO        false
#define _N_BAR_LAST_LAYER   false
#define _PERC_35_LAYER      false

inline bool SAACheck(const std::shared_ptr<DmpEvtHeader> evt_header, const std::shared_ptr<DmpFilterOrbit> pFilter) {
    return pFilter->IsInSAA(evt_header->GetSecond()) ? false : true;
}

inline bool checkBGOreco(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries) {
        return ((bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0) || (bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0)) ? false : true;
    }

inline bool check_event(
    const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries, 
    const std::shared_ptr<DmpEvtHeader> evt_header, 
    const std::shared_ptr<DmpFilterOrbit> pFilter, 
    const bool mc) {
        if (mc) {
            return checkBGOreco(bgoRec_slope, bgoRec_intercept, simu_primaries);
        }
        else {
            if (SAACheck(evt_header, pFilter))
                return checkBGOreco(bgoRec_slope, bgoRec_intercept, simu_primaries);
            else return false;
        }
    }

inline const std::vector<double> get_bgo_slope(std::shared_ptr<DmpEvtBgoRec> bgorec) {
    return std::vector<double> {bgorec->GetSlopeXZ(), bgorec->GetSlopeYZ()};
}

inline const std::vector<double> get_bgo_intercept(std::shared_ptr<DmpEvtBgoRec> bgorec) {
    return std::vector<double> {bgorec->GetInterceptXZ(), bgorec->GetInterceptYZ()};
}

void focus(const in_pars &input_pars) {

    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose, input_pars.mc_flag);
    std::shared_ptr<energy_config> econfig = std::make_shared<energy_config>(input_pars.config_wd);
    std::shared_ptr<config> _config = std::make_shared<config>(input_pars.config_wd);

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
    if (input_pars.mc_flag) {
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
	std::shared_ptr<TClonesArray> stktracks = std::make_shared<TClonesArray>("DmpStkTrack");
	if (fStkKalmanTracksFound)
		evtch->SetBranchAddress("StkKalmanTracks", &stktracks);

    // Register PSD container
	std::shared_ptr<DmpEvtPsdHits> psdhits = std::make_shared<DmpEvtPsdHits>();
	evtch->SetBranchAddress("DmpPsdHits", &psdhits);

    // Register the orbit filter for DATA SAA check
    std::shared_ptr<DmpFilterOrbit> pFilter;
    if (input_pars.rawdata_flag) {
        gIOSvc->Set("OutData/NoOutput", "True");
	    gIOSvc->Initialize();
        pFilter = std::make_shared<DmpFilterOrbit>("EventHeader");
        pFilter->ActiveMe(); // Call this function to calculate SAA through House Keeping Data
    }

    auto min_evt_energy = econfig->GetMinEvtEnergy();
    auto max_evt_energy = econfig->GetMaxEvtEnergy();

    double evt_corr_energy {0};
    double evt_energy {0};
    double evt_corr_energy_gev {0};
    double evt_energy_gev {0};
    double gev {0.001};
    unsigned int perc {10};

    auto nevents {evtch->GetEntries()};

    if (input_pars.verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;

    std::vector<unsigned int> focus_events;
    std::vector<DmpStkTrack> best_track;

    const int number_of_track_clusters {3};
    
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

        bool isTrackInteresting {false};
        DmpStkTrack tmp_best_track;

        if (evt_corr_energy_gev>=(_config->GetEventDisplayConfig()).evt_min_energy && evt_corr_energy_gev<=(_config->GetEventDisplayConfig()).evt_max_energy) {
            
            if (check_event(get_bgo_slope(bgorec), get_bgo_intercept(bgorec), simu_primaries, evt_header, pFilter, input_pars.mc_flag)) {
                
                if (_BEST_TRACK) {
                    std::tie(isTrackInteresting, tmp_best_track) = 
                        get_best_track(
                            bgorec,
                            evt_header,
                            get_bgo_slope(bgorec),
                            get_bgo_intercept(bgorec),
                            bgohits,
                            stkclusters,
                            stktracks,
                            evt_energy,
                            _config);

                    if (isTrackInteresting) {
                        focus_events.push_back(evIdx);
                        best_track.push_back(tmp_best_track);
                    }
                }

                if (_BGO_EDGE) {
                    if (isEventHittingEdgeBars(bgohits, bgorec, evt_header, stkclusters, stktracks, psdhits, evt_energy, _config))
                        focus_events.push_back(evIdx);
                }

                if (_OUTSIDE_BGO) {
                    if (isEventOutsideBGOFiducial(bgohits, bgorec, evt_header, evt_energy, _config))
                        focus_events.push_back(evIdx);
                }

                if (_N_BAR_LAST_LAYER) {
                    if (isNumberOfLastLayerBarsTooHigh(bgohits, bgorec, evt_header, stkclusters, stktracks, psdhits, evt_energy, _config))
                        focus_events.push_back(evIdx);
                }

                if (_PERC_35_LAYER) {
                    if (isEnergyPercentageGT35(bgohits, bgorec, evt_header, stkclusters, stktracks, psdhits, evt_energy, _config))
                        focus_events.push_back(evIdx);
                }
            }
        }
    }

    if (input_pars.verbose) std::cout << "\n\nNumber of selected events: " << focus_events.size() << std::endl;

    std::string output_file = input_pars.output_wd + "/focus.root";
    TFile* focus_tuple_output = TFile::Open(output_file.c_str(), "RECREATE");
    auto focus_tree = static_cast<TTree*>(evtch->CloneTree(0));
    DmpStkTrack track;
    //focus_tree->Branch("STKBestTrack", &track);

    for (unsigned int ev {0}; ev<focus_events.size(); ++ev) {
        evtch->GetEvent(focus_events[ev]);
        //track = best_track[ev];
        focus_tree->Fill();
    }

    focus_tree->Write();
    focus_tuple_output->Close();

    if (input_pars.verbose) std::cout << "\n\nOutput file has been written [" << output_file << "]\n\n"; 
}