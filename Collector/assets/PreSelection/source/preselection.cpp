#include "bgo.h"
#include "stk.h"
#include "psd.h"
#include "chain.h"
#include "histos.h"
#include "charge.h"
#include "config.h"
#include "efficiency.h"
#include "preselection.h"
#include "energy_config.h"
#include "lateral_showering.h"

#include <iostream>
#include <memory>
#include <vector>

#include "Dmp/DmpGeoStruct.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtPsdHits.h"
#include "DmpSimuTrajectory.h"
#include "DmpEvtSimuPrimaries.h"

#include "DmpIOSvc.h"
#include "DmpCore.h"
#include "DmpFilterOrbit.h"

#include "TChain.h"
#include "TClonesArray.h"

#include "TFile.h"
#include "TTree.h"

inline bool SAACheck(const std::shared_ptr<DmpEvtHeader> evt_header, const std::shared_ptr<DmpFilterOrbit> pFilter) {
    return pFilter->IsInSAA(evt_header->GetSecond()) ? false : true;
}

inline bool checkBGOreco(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries) {
        
#if 0        
        bool bgoreco_status {false};
        int position_sensitivity {20};

        if ((bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0) || (bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0)) {
            if (simu_primaries!=nullptr) {
                TVector3 orgPosition;
                orgPosition.SetX(simu_primaries->pv_x);
                orgPosition.SetY(simu_primaries->pv_y);
                orgPosition.SetZ(simu_primaries->pv_z);

                TVector3 orgMomentum;
                orgMomentum.SetX(simu_primaries->pvpart_px);
                orgMomentum.SetY(simu_primaries->pvpart_py);
                orgMomentum.SetZ(simu_primaries->pvpart_pz);

                std::vector<double> slope(2, 0);
                std::vector<double> intercept(2, 0);

                slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
                slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
                intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
                intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

                double actual_X = slope[0] * BGO_TopZ + intercept[0];
                double actual_Y = slope[1] * BGO_TopZ + intercept[1];
                double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
                double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

                if (fabs(actual_X - topX) < position_sensitivity && fabs(actual_Y - topY) < position_sensitivity)
                    bgoreco_status = true;
            }
        }
        else
            bgoreco_status = true;

        return bgoreco_status;
#endif
        
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

void preselection(const in_pars &input_pars) {
    
    auto evtch = GetFileChain(input_pars.input_path, input_pars.verbose, input_pars.mc_flag);
    std::shared_ptr<energy_config> econfig = std::make_shared<energy_config>(input_pars.config_wd);
    std::shared_ptr<config> cuts_config = std::make_shared<config>(input_pars.config_wd);

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
    double simu_energy {0};
    double simu_energy_gev {0};
    double gev {0.001};
    unsigned int perc {10};

    auto nevents {evtch->GetEntries()};

    std::shared_ptr<histos> ps_histos = std::make_shared<histos>(econfig, input_pars.mc_flag);

    if (input_pars.verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;
    
    //std::vector<DmpStkTrack> myBestTrack (nevents);

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
        ps_histos->SetWeight(simu_primaries, evt_corr_energy_gev);

        if (input_pars.mc_flag) {
            simu_energy = simu_primaries->pvpart_ekin;
            simu_energy_gev = simu_energy*gev;
        }

        if (evt_corr_energy_gev>=min_evt_energy && evt_corr_energy_gev<=max_evt_energy) {
            
            if (check_event(get_bgo_slope(bgorec), get_bgo_intercept(bgorec), simu_primaries, evt_header, pFilter, input_pars.mc_flag)) {
                
                // BGO distributions without any prior cut
                bgo_distributions(
                    bgohits, 
                    bgorec, 
                    evt_header, 
                    simu_primaries,
                    evt_energy,
                    evt_corr_energy,
                    evt_energy_gev,
                    evt_corr_energy_gev, 
                    ps_histos,
                    cuts_config);

                // BGO distributions after the BGO fiducial volume cut
                bgofiducial_distributions(
                    bgohits, 
                    bgorec, 
                    evt_header, 
                    simu_primaries,
                    evt_energy,
                    evt_corr_energy,
                    evt_energy_gev,
                    evt_corr_energy_gev, 
                    ps_histos,
                    cuts_config);

                // BGO distributions as last cut
                bgofiducial_distributions_lastcut(
                    bgohits, 
                    bgorec, 
                    evt_header,
                    stkclusters,
                    stktracks,
                    psdhits,
                    simu_primaries,
                    evt_energy,
                    evt_corr_energy,
                    evt_energy_gev,
                    evt_corr_energy_gev, 
                    ps_histos,
                    cuts_config);

                // Remove lateral and large showering events
                lateral_showering_distributions(
                    bgohits, 
                    bgorec, 
                    evt_header,
                    evt_energy, 
                    evt_corr_energy, 
                    evt_energy_gev, 
                    evt_corr_energy_gev,
                    ps_histos,
                    cuts_config);

                /*
                // Remove not showering events using the sumrms low energy cut
                lateral_showering_distributions_after_sumrms_clean_cut(
                    bgohits, 
                    bgorec, 
                    evt_header,
                    evt_energy, 
                    evt_corr_energy, 
                    evt_energy_gev, 
                    evt_corr_energy_gev,
                    ps_histos,
                    cuts_config);
                */

                // Remove lateral and large showering events as last cut
                lateral_showering_distributions_lastcut(
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
                    ps_histos,
                    cuts_config);
                
                // STK distributions and charges
                stk_distributions(
                    bgohits, 
                    bgorec, 
                    evt_header, 
                    stkclusters, 
                    stktracks, 
                    evt_energy, 
                    evt_corr_energy, 
                    evt_energy_gev, 
                    evt_corr_energy_gev, 
                    ps_histos,
                    cuts_config);

                // PSD-STK match distributions and charges
                psd_stk_match_distributions(
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
                    ps_histos,
                    cuts_config);
                
                // PSD-STK charge distributions
                charge_distributions(
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
                    ps_histos,
                    cuts_config);

                buildEfficiencies(
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
                    simu_energy_gev,
                    ps_histos,
                    cuts_config);
            }
        }

    }

    ps_histos->Write(input_pars.output_wd, input_pars.verbose);
}