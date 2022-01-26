#include "cuts.h"
#include "charge.h"
#include "trigger.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"

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
    std::shared_ptr<histos> ps_histos,
    std::shared_ptr<config> cuts_config) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        std::unique_ptr<DmpPsdContainer> psdVault = std::make_unique<DmpPsdContainer>();

        best_track event_best_track;
        psd_cluster_match clu_matching;
        
        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);
        psdVault->scanPSDHits(psdhits, (cuts_config->GetCutsConfig()).psd_min_energy);

        auto weight {ps_histos->GetWeight()};
        
        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);
                if (nbarlayer13_cut) {
                    auto maxrms_cut = maxRms_cut(bgoVault->GetELayer(), bgoVault->GetRmsLayer(), evt_energy, (cuts_config->GetCutsConfig()).bgo_shower_width);

                    if (maxrms_cut) {
                        auto trackselection_cut = track_selection_cut(
                        bgorec, 
                        bgoVault->GetBGOslope(), 
                        bgoVault->GetBGOintercept(), 
                        bgohits, 
                        stkclusters, 
                        stktracks, 
                        event_best_track,
                        (cuts_config->GetCutsConfig()).STK_BGO_delta_position,
                        (cuts_config->GetCutsConfig()).STK_BGO_delta_track,
                        (cuts_config->GetCutsConfig()).track_X_clusters,
                        (cuts_config->GetCutsConfig()).track_Y_clusters,
                        (cuts_config->GetCutsConfig()).track_X_holes,
                        (cuts_config->GetCutsConfig()).track_Y_holes);
                
                        if (trackselection_cut) {
                            auto psdstkmatch_cut = psd_stk_match_cut(
                            bgoVault->GetBGOslope(), 
                            bgoVault->GetBGOintercept(),
                            psdVault->getPsdClusterIdxBegin(),
                            psdVault->getPsdClusterZ(),
                            psdVault->getPsdClusterMaxECoo(),
                            event_best_track,
                            clu_matching,
                            (cuts_config->GetCutsConfig()).STK_PSD_delta_position);

                            if (psdstkmatch_cut) {
                                
                                auto psdcharge_cut = psd_charge_cut(
                                    psdVault->getPsdClusterMaxE(),
                                    event_best_track,
                                    clu_matching,
                                    (cuts_config->GetCutsConfig()).PSD_sharge_sum,
                                    (cuts_config->GetCutsConfig()).PSD_single_charge);

                                auto stkcharge_cut = stk_charge_cut(stkclusters, event_best_track, (cuts_config->GetCutsConfig()).STK_single_charge);

                                auto stk_charges = stk_charge(stkclusters, event_best_track, ps_histos);
                                auto psd_charges = psd_charge(psdVault->getPsdClusterMaxE(), event_best_track, clu_matching, ps_histos);

                                ps_histos->h_PSD_STK_X_charge->Fill(stk_charges[0], psd_charges[0], weight);
                                ps_histos->h_PSD_STK_Y_charge->Fill(stk_charges[1], psd_charges[1], weight);
                                ps_histos->h_PSD_STK_charge->Fill(stk_charges[2], psd_charges[2], weight);

                                if (abs(psd_charges[0]-psd_charges[1])<=2) {
                                    ps_histos->h_PSD_charge_X_cleaned->Fill(psd_charges[0], weight);
                                    ps_histos->h_PSD_charge_Y_cleaned->Fill(psd_charges[1], weight);
                                    ps_histos->h_PSD_charge_cleaned->Fill(psd_charges[2], weight);
                                    ps_histos->h_PSD_charge_2D_cleaned->Fill(psd_charges[0], psd_charges[1], weight);
                                }
                                
                                if (psdcharge_cut) {
                                    stk_charge(stkclusters, event_best_track, ps_histos, true);

                                    ps_histos->h_PSD_STK_X_charge_after_PSD_charge_cut->Fill(stk_charges[0], psd_charges[0], weight);
                                    ps_histos->h_PSD_STK_Y_charge_after_PSD_charge_cut->Fill(stk_charges[1], psd_charges[1], weight);
                                    ps_histos->h_PSD_STK_charge_after_PSD_charge_cut->Fill(stk_charges[2], psd_charges[2], weight);
                                }

                                if (stkcharge_cut) {
                                    psd_charge(psdVault->getPsdClusterMaxE(), event_best_track, clu_matching, ps_histos, true);

                                    ps_histos->h_PSD_STK_X_charge_after_STK_charge_cut->Fill(stk_charges[0], psd_charges[0], weight);
                                    ps_histos->h_PSD_STK_Y_charge_after_STK_charge_cut->Fill(stk_charges[1], psd_charges[1], weight);
                                    ps_histos->h_PSD_STK_charge_after_STK_charge_cut->Fill(stk_charges[2], psd_charges[2], weight);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

void stk_charge_explore(
    const std::shared_ptr<TClonesArray> stkclusters,
    DmpStkTrack* track,
    std::shared_ptr<histos> ps_histos) {
        
        auto weight {ps_histos->GetWeight()};

        double cluster_chargeX = -999;
        double cluster_chargeY = -999;

        // Charge correction
        auto track_correction = track->getDirection().CosTheta();

        // Compute charges
        for (auto clIdx = 0; clIdx < track->GetNPoints(); ++clIdx) {
            auto cluster_x = track->GetClusterX(clIdx, stkclusters.get());
            auto cluster_y = track->GetClusterY(clIdx, stkclusters.get());
            if (cluster_x && !cluster_x->getPlane())
                cluster_chargeX = sqrt(cluster_x->getEnergy() * track_correction);
            if (cluster_y && !cluster_y->getPlane())
                cluster_chargeY = sqrt(cluster_y->getEnergy() * track_correction);
        }

        // Check charges
        if (cluster_chargeX != -999 && cluster_chargeY != -999) {
            ps_histos->h_STK_charge_X->Fill(cluster_chargeX, weight);
            ps_histos->h_STK_charge_Y->Fill(cluster_chargeY, weight);
            ps_histos->h_STK_charge->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
            ps_histos->h_STK_charge_2D->Fill(cluster_chargeX, cluster_chargeY, weight);

            if (track->getNhitX() == 3 && track->getNhitY() == 3) {
                ps_histos->h_STK_charge_X_3_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_3_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_3_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_3_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (track->getNhitX() == 4 && track->getNhitY() == 4) {
                ps_histos->h_STK_charge_X_4_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_4_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_4_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_4_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (track->getNhitX() == 5 && track->getNhitY() == 5) {
                ps_histos->h_STK_charge_X_5_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_5_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_5_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_5_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
            else if (track->getNhitX() == 6 && track->getNhitY() == 6) {
                ps_histos->h_STK_charge_X_6_clusters->Fill(cluster_chargeX, weight);
                ps_histos->h_STK_charge_Y_6_clusters->Fill(cluster_chargeY, weight);
                ps_histos->h_STK_charge_6_clusters->Fill(0.5 * (cluster_chargeX + cluster_chargeY), weight);
                ps_histos->h_STK_charge_2D_6_clusters->Fill(cluster_chargeX, cluster_chargeY, weight);
            }
        }
    }

std::vector<double> stk_charge(
    const std::shared_ptr<TClonesArray> stkclusters, 
    best_track &event_best_track, 
    std::shared_ptr<histos> ps_histos,
    const bool lastcut) {
    
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

        return std::vector<double> {cluster_chargeX, cluster_chargeY, 0.5 * (cluster_chargeX + cluster_chargeY)};
    }

std::vector<double> psd_charge(
	const std::vector<std::vector<double>> psdCluster_maxE,
    best_track &event_best_track,
    psd_cluster_match &clu_matching,
	std::shared_ptr<histos> ps_histos,
    const bool lastcut) {
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

        return std::vector<double> {psd_chargeX, psd_chargeY, 0.5*(psd_chargeX+psd_chargeY)};
    }

void psd_charge_explore(
	const std::vector<std::vector<double>> psdCluster_maxE,
    best_track &event_best_track,
    psd_cluster_match &clu_matching,
	std::shared_ptr<histos> ps_histos) {
        
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
        
        ps_histos->h_PSD_charge_X->Fill(psd_chargeX, weight);
        ps_histos->h_PSD_charge_Y->Fill(psd_chargeY, weight);
        ps_histos->h_PSD_charge->Fill(0.5*(psd_chargeX+psd_chargeY), weight);
        ps_histos->h_PSD_charge_2D->Fill(psd_chargeX, psd_chargeY, weight);
        ps_histos->h_PSD_sum_of_XY_charges->Fill(psd_chargeX+psd_chargeY, weight);
        
    }