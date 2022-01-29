#include "cuts.h"
#include "trigger.h"
#include "efficiency.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"

void buildEfficiencies(
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
    const double simu_energy_gev,
    std::shared_ptr<histos> ps_histos,
    std::shared_ptr<config> cuts_config) {

        // PSD charge efficiency
        psd_charge_efficiency(
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
        
        // PSD-STK matching efficiency
        psd_stk_matching_efficiency(
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

        // Track selection efficiency
        stk_track_efficiency(
            bgohits, 
            bgorec, 
            evt_header, 
            stkclusters, 
            stktracks,  
            evt_energy, 
            evt_corr_energy, 
            evt_energy_gev, 
            evt_corr_energy_gev,
            simu_energy_gev,
            ps_histos, 
            cuts_config);

        // nbarlayer13 selection efficiency
        nbarlayer13_efficiency(
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

        // maxrms selection efficiency
        maxrms_efficiency(
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

void psd_charge_efficiency(
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
    const double simu_energy_gev, 
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

                                if (psdcharge_cut)
                                    ps_histos->h_psdcharge_lastcut_pass->Fill(simu_energy_gev);
                                ps_histos->h_psdcharge_lastcut->Fill(simu_energy_gev);
                            }
                        }
                    }
                }
            }
        }
    }

void psd_stk_matching_efficiency(
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
    const double simu_energy_gev,
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

                            if (psdstkmatch_cut)
                                ps_histos->h_psdstkmatch_lastcut_pass->Fill(simu_energy_gev);
                            ps_histos->h_psdstkmatch_lastcut->Fill(simu_energy_gev);
                        }
                    }
                }
            }
        }
    }

void stk_track_efficiency(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev,
    const double simu_energy_gev,
    std::shared_ptr<histos> ps_histos,
    std::shared_ptr<config> cuts_config) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        
        best_track event_best_track;

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters);

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

                        if (trackselection_cut)
                            ps_histos->h_trackselection_lastcut_pass->Fill(simu_energy_gev);
                        ps_histos->h_trackselection_lastcut->Fill(simu_energy_gev);
                    }
                }
            }
        }
    }

void nbarlayer13_efficiency(
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
    const double simu_energy_gev, 
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

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;
                
            if (bgofiducial_cut) {
                
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

                            if (psdcharge_cut) {
                            
                                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);

                                if (nbarlayer13_cut)
                                    ps_histos->h_nbarlayer13_lastcut_pass->Fill(simu_energy_gev);
                                ps_histos->h_nbarlayer13_lastcut->Fill(simu_energy_gev);

                            }
                        }
                    }
                }
            }
        }
    }

void maxrms_efficiency(
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
    const double simu_energy_gev, 
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

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;
                
            if (bgofiducial_cut) {
                
                auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarNumber(13), evt_energy);

                if (nbarlayer13_cut) {

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

                            if (psdcharge_cut) {
                            
                                auto maxrms_cut = maxRms_cut(bgoVault->GetELayer(), bgoVault->GetRmsLayer(), evt_energy, (cuts_config->GetCutsConfig()).bgo_shower_width);

                                if (maxrms_cut)
                                    ps_histos->h_maxrms_lastcut_pass->Fill(simu_energy_gev);
                                ps_histos->h_maxrms_lastcut->Fill(simu_energy_gev);

                            }
                        }
                    }
                }
            }
        }
    }