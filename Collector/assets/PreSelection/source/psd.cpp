#include "psd.h"
#include "cuts.h"
#include "charge.h"
#include "trigger.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"

void psd_stk_match_distributions(
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
                            (cuts_config->GetCutsConfig()).track_Y_clusters);
                    
                        if (trackselection_cut) {
                            psd_fiducial_stk_match(
                                bgoVault->GetBGOslope(),
                                bgoVault->GetBGOintercept(),
                                psdVault->getPsdClusterIdxBegin(),
                                psdVault->getPsdClusterZ(),
                                psdVault->getPsdClusterMaxECoo(),
                                event_best_track,
                                clu_matching,
                                ps_histos,
                                evt_corr_energy_gev);

                            psd_stk_match(
                                bgoVault->GetBGOslope(),
                                bgoVault->GetBGOintercept(),
                                psdVault->getPsdClusterIdxBegin(),
                                psdVault->getPsdClusterZ(),
                                psdVault->getPsdClusterMaxECoo(),
                                event_best_track,
                                clu_matching,
                                ps_histos,
                                evt_corr_energy_gev);

                            psd_charge_explore(psdVault->getPsdClusterMaxE(), event_best_track, clu_matching, ps_histos);
                        }
                    }
                }
            }
        }
        
    }

void psd_stk_match(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::vector<std::vector<short>> psdCluster_idxBeg,
	const std::vector<std::vector<double>> psdCluster_Z,
	const std::vector<std::vector<double>> psdCluster_maxEcoordinate,
    const best_track &event_best_track,
    psd_cluster_match &clu_matching,
    std::shared_ptr<histos> ps_histos,
    const double evt_corr_energy_gev) {
        
        for (int nLayer = 0; nLayer < DAMPE_psd_nLayers; ++nLayer) {
            for (unsigned int iclu = 0; iclu < psdCluster_idxBeg[nLayer].size(); ++iclu) {
                bool IsMeasuringX = nLayer % 2;
                double hitZ = psdCluster_Z[nLayer][iclu];
                double thisCoord = psdCluster_maxEcoordinate[nLayer][iclu];

                // Get distance between actual coordinate and BGO rec coordinate
                double projCoord_bgoRec = IsMeasuringX ? bgoRec_slope[0] * hitZ + bgoRec_intercept[0] : bgoRec_slope[1] * hitZ + bgoRec_intercept[1];
                double dX_bgoRec = thisCoord - projCoord_bgoRec;
                if (fabs(dX_bgoRec) < fabs(clu_matching.dxCloPsdClu_bgoRec[nLayer])) {
                    clu_matching.dxCloPsdClu_bgoRec[nLayer] = dX_bgoRec;
                    clu_matching.icloPsdClu_bgoRec[nLayer] = iclu;
                }

                // Get distance between actual coordinate and best track coordinate
                double projCoord_track = IsMeasuringX ? event_best_track.track_slope[0] * hitZ + event_best_track.track_intercept[0] : event_best_track.track_slope[1] * hitZ + event_best_track.track_intercept[1];
                double dX_track = thisCoord - projCoord_track;

                if (fabs(dX_track) < fabs(clu_matching.dxCloPsdClu_track[nLayer])) {
                    clu_matching.dxCloPsdClu_track[nLayer] = dX_track;
                    clu_matching.icloPsdClu_track[nLayer] = iclu;
                }
            }
        }

        auto weight {ps_histos->GetWeight()};
        
        ps_histos->h_PSD_STK_X_match_energy_int->Fill(clu_matching.dxCloPsdClu_track[0], weight);
        ps_histos->h_PSD_STK_Y_match_energy_int->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        
        if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250) {
            ps_histos->h_PSD_STK_X_match_100_250->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_100_250->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500) {
            ps_histos->h_PSD_STK_X_match_250_500->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_250_500->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000) {
            ps_histos->h_PSD_STK_X_match_500_1000->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_500_1000->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<5000) {
            ps_histos->h_PSD_STK_X_match_1000_5000->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_1000_5000->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=5000) {
            ps_histos->h_PSD_STK_X_match_5000->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_5000->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        
        ps_histos->h_PSD_X_clusters->Fill(evt_corr_energy_gev, psdCluster_idxBeg[0].size(), weight);
        ps_histos->h_PSD_Y_clusters->Fill(evt_corr_energy_gev, psdCluster_idxBeg[1].size(), weight);
    }

void psd_fiducial_stk_match(
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept,
	const std::vector<std::vector<short>> psdCluster_idxBeg,
	const std::vector<std::vector<double>> psdCluster_Z,
	const std::vector<std::vector<double>> psdCluster_maxEcoordinate,
    const best_track &event_best_track,
    psd_cluster_match &clu_matching,
    std::shared_ptr<histos> ps_histos,
    const double evt_corr_energy_gev) {
        
        const double PSD_fiducial = 410;

        for (int nLayer = 0; nLayer < DAMPE_psd_nLayers; ++nLayer) {
            for (unsigned int iclu = 0; iclu < psdCluster_idxBeg[nLayer].size(); ++iclu) {
                bool IsMeasuringX = nLayer % 2;
                double hitZ = psdCluster_Z[nLayer][iclu];
                double thisCoord = psdCluster_maxEcoordinate[nLayer][iclu];

                // Get distance between actual coordinate and BGO rec coordinate
                double projCoord_bgoRec = IsMeasuringX ? bgoRec_slope[0] * hitZ + bgoRec_intercept[0] : bgoRec_slope[1] * hitZ + bgoRec_intercept[1];
                double dX_bgoRec = thisCoord - projCoord_bgoRec;
                if (fabs(dX_bgoRec) < fabs(clu_matching.dxCloPsdClu_bgoRec[nLayer])) {
                    clu_matching.dxCloPsdClu_bgoRec[nLayer] = dX_bgoRec;
                    clu_matching.icloPsdClu_bgoRec[nLayer] = iclu;
                }

                // Get distance between actual coordinate and best track coordinate
                double projCoord_track = IsMeasuringX ? event_best_track.track_slope[0] * hitZ + event_best_track.track_intercept[0] : event_best_track.track_slope[1] * hitZ + event_best_track.track_intercept[1];
                if (fabs(projCoord_track) > PSD_fiducial) continue;

                double dX_track = thisCoord - projCoord_track;

                if (fabs(dX_track) < fabs(clu_matching.dxCloPsdClu_track[nLayer])) {
                    clu_matching.dxCloPsdClu_track[nLayer] = dX_track;
                    clu_matching.icloPsdClu_track[nLayer] = iclu;
                }
            }
        }

        auto weight {ps_histos->GetWeight()};
        
        ps_histos->h_PSD_STK_X_match_energy_int_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[0], weight);
        ps_histos->h_PSD_STK_Y_match_energy_int_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        
        if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250) {
            ps_histos->h_PSD_STK_X_match_100_250_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_100_250_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500) {
            ps_histos->h_PSD_STK_X_match_250_500_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_250_500_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000) {
            ps_histos->h_PSD_STK_X_match_500_1000_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_500_1000_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<5000) {
            ps_histos->h_PSD_STK_X_match_1000_5000_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_1000_5000_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }
        else if (evt_corr_energy_gev>=5000) {
            ps_histos->h_PSD_STK_X_match_5000_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[0], weight);
            ps_histos->h_PSD_STK_Y_match_5000_psd_fiducial->Fill(clu_matching.dxCloPsdClu_track[1], weight);
        }   
    }