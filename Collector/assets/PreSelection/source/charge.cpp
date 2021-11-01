#include "cuts.h"
#include "charge.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"

inline bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
    return trigger_general;
}

inline bool stk_charge(
    const std::shared_ptr<TClonesArray> stkclusters, 
    best_track &event_best_track, 
    std::shared_ptr<histos> ps_histos,
    const bool lastcut = false) {
    
    auto weight {ps_histos->GetWeight()};

    double cluster_chargeX = -999;
	double cluster_chargeY = -999;

	// Charge correction
	auto track_correction = (event_best_track.myBestTrack).getDirection().CosTheta();

	// Compute charges
	for (auto clIdx = 0; clIdx < event_best_track.n_points; ++clIdx) {
		auto cluster_x = (event_best_track.myBestTrack).GetClusterX(clIdx, stkclusters.get());
		auto cluster_y = (event_best_track.myBestTrack).GetClusterY(clIdx, stkclusters.get());
		if (cluster_x && !cluster_x->getPlane())
			cluster_chargeX = sqrt(cluster_x->getEnergy() * track_correction);
		if (cluster_y && !cluster_y->getPlane())
			cluster_chargeY = sqrt(cluster_y->getEnergy() * track_correction);
	}

	// Check charges
	if (cluster_chargeX != -999 && cluster_chargeY != -999) {
        if (!lastcut) {
            ps_histos->h_STK_charge_X_nocut->Fill(cluster_chargeX, weight);
            ps_histos->h_STK_charge_Y_nocut->Fill(cluster_chargeY, weight);
            ps_histos->h_STK_charge_nocut->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
            ps_histos->h_STK_charge_2D_nocut->Fill(cluster_chargeX, cluster_chargeY, weight);
        }
        else {
            ps_histos->h_STK_charge_X_PSD_charge_cut->Fill(cluster_chargeX, weight);
            ps_histos->h_STK_charge_Y_PSD_charge_cut->Fill(cluster_chargeY, weight);
            ps_histos->h_STK_charge_PSD_charge_cut->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
            ps_histos->h_STK_charge_2D_PSD_charge_cut->Fill(cluster_chargeX, cluster_chargeY, weight);
        }

    }
}

inline bool psd_charge(
	const std::vector<std::vector<double>> psdCluster_maxE,
    best_track &event_best_track,
    psd_cluster_match &clu_matching,
	std::shared_ptr<histos> ps_histos,
    const bool lastcut = false) {

        double psd_chargeX, psd_chargeY;

        // Charge correction
        auto track_correction = (event_best_track.myBestTrack).getDirection().CosTheta();

        // Get Y charge
        if (clu_matching.icloPsdClu_track[0] > -1) {
            auto energy_ClusterYTrack = psdCluster_maxE[0][clu_matching.icloPsdClu_track[0]];
            auto energy_ClusterYTrack_corr = track_correction * energy_ClusterYTrack;
            psd_chargeY = sqrt(energy_ClusterYTrack_corr);
        }

        // Get X charge
        if (clu_matching.icloPsdClu_track[1] > -1)
        {
            auto energy_ClusterXTrack = psdCluster_maxE[1][clu_matching.icloPsdClu_track[1]];
            auto energy_ClusterXTrack_corr = track_correction * energy_ClusterXTrack;
            psd_chargeX = sqrt(energy_ClusterXTrack_corr);
        }

        auto weight {ps_histos->GetWeight()};

        if (!lastcut) {
            ps_histos->h_PSD_charge_X_nocut->Fill(psd_chargeX, weight);
            ps_histos->h_PSD_charge_Y_nocut->Fill(psd_chargeY, weight);
            ps_histos->h_PSD_charge_nocut->Fill(0.5*(psd_chargeX+psd_chargeY), weight);
            ps_histos->h_PSD_charge_2D_nocut->Fill(psd_chargeX, psd_chargeY, weight);
            ps_histos->h_PSD_sum_of_XY_charges_nocut->Fill(psd_chargeX+psd_chargeY, weight);
        }
        else {
            ps_histos->h_PSD_charge_X_STK_charge_cut->Fill(psd_chargeX, weight);
            ps_histos->h_PSD_charge_Y_STK_charge_cut->Fill(psd_chargeY, weight);
            ps_histos->h_PSD_charge_STK_charge_cut->Fill(0.5*(psd_chargeX+psd_chargeY), weight);
            ps_histos->h_PSD_charge_2D_STK_charge_cut->Fill(psd_chargeX, psd_chargeY, weight);
            ps_histos->h_PSD_sum_of_XY_charges_STK_charge_cut->Fill(psd_chargeX+psd_chargeY, weight);
        }
       
    }

void charge_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        std::unique_ptr<DmpPsdContainer> psdVault = std::make_unique<DmpPsdContainer>();

        best_track event_best_track;
        psd_cluster_match clu_matching;

        double bgo_layer_min_energy     {0};    // Minimum energy per BGO layer
        double bgo_max_energy_ratio     {0.35}; // Maximum energy ratio per layer
        double bgo_shower_axis_delta    {280};  // BGO maximum shower axis delta (mm)
        double bgo_shower_width         {100};  // BGO maximum shower width (mm)
        double STK_BGO_delta_position   {40};   // Linear distance between STK and BGO projections
        double STK_BGO_delta_track      {10};   // Angular distance between BGO/STK tracks (deg)
        int track_X_clusters            {4};    // Number of requested clusters on X tracks
        int track_Y_clusters            {4};    // Number of requested clusters on Y tracks
        double STK_PSD_delta_position   {40};   // Linear distance between STK progection and PSD seed strip
        double psd_min_energy           {0};    // Minimum energy per PSD bar
        double PSD_sharge_sum           {10};   // Sum of PSD charges on X and Y views
        double PSD_single_charge        {2.6};  // PSD charge cut on single view
        double STK_single_charge        {40};   // STK charge cut on single view

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);
        psdVault->scanPSDHits(psdhits, psd_min_energy);

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);
                if (nbarlayer13_cut) {
                    auto maxrms_cut = maxRms_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetRmsLayer(), evt_energy, bgo_shower_width);

                    if (maxrms_cut) {
                        auto trackselection_cut = track_selection_cut(
                        bgorec, 
                        bgoVault->GetBGOslope(), 
                        bgoVault->GetBGOintercept(), 
                        bgohits, 
                        stkclusters, 
                        stktracks, 
                        event_best_track,
                        STK_BGO_delta_position,
                        STK_BGO_delta_track,
                        track_X_clusters,
                        track_Y_clusters);
                
                        if (trackselection_cut) {
                            auto psdstkmatch_cut = psd_stk_match_cut(
                            bgoVault->GetBGOslope(), 
                            bgoVault->GetBGOintercept(),
                            psdVault->getPsdClusterIdxBegin(),
                            psdVault->getPsdClusterZ(),
                            psdVault->getPsdClusterMaxECoo(),
                            event_best_track,
                            clu_matching,
                            STK_PSD_delta_position);

                            if (psdstkmatch_cut) {
                                
                                auto psdcharge_cut = psd_charge_cut(
                                    psdVault->getPsdClusterMaxE(),
                                    event_best_track,
                                    clu_matching,
                                    PSD_sharge_sum,
                                    PSD_single_charge);

                                auto stkcharge_cut = stk_charge_cut(stkclusters, event_best_track, STK_single_charge);

                                stk_charge(stkclusters, event_best_track, ps_histos);
                                psd_charge(psdVault->getPsdClusterMaxE(), event_best_track, clu_matching, ps_histos);

                                if (psdcharge_cut)
                                    stk_charge(stkclusters, event_best_track, ps_histos, true);

                                if (stkcharge_cut)
                                    psd_charge(psdVault->getPsdClusterMaxE(), event_best_track, clu_matching, ps_histos, true);
                            }
                        }
                    }
                }
            }
        }
    }