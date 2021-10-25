#include "psd.h"
#include "cuts.h"

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

inline void psd_stk_match(
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

        ps_histos->PSD_STK_X_match_energy_int->Fill(clu_matching.dxCloPsdClu_track[0]);
        ps_histos->PSD_STK_Y_match_energy_int->Fill(clu_matching.dxCloPsdClu_track[1]);
        
        if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<=250) {
            ps_histos->PSD_STK_X_match_100_250->Fill(clu_matching.dxCloPsdClu_track[0]);
            ps_histos->PSD_STK_Y_match_100_250->Fill(clu_matching.dxCloPsdClu_track[1]);
        }
        else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<=500) {
            ps_histos->PSD_STK_X_match_250_500->Fill(clu_matching.dxCloPsdClu_track[0]);
            ps_histos->PSD_STK_X_match_250_500->Fill(clu_matching.dxCloPsdClu_track[1]);
        }
        else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<=1000) {
            ps_histos->PSD_STK_X_match_500_1000->Fill(clu_matching.dxCloPsdClu_track[0]);
            ps_histos->PSD_STK_X_match_500_1000->Fill(clu_matching.dxCloPsdClu_track[1]);
        }
        else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<=5000) {
            ps_histos->PSD_STK_X_match_1000_5000->Fill(clu_matching.dxCloPsdClu_track[0]);
            ps_histos->PSD_STK_X_match_1000_5000->Fill(clu_matching.dxCloPsdClu_track[1]);
        }
        else {
            ps_histos->PSD_STK_X_match_5000->Fill(clu_matching.dxCloPsdClu_track[0]);
            ps_histos->PSD_STK_Y_match_5000->Fill(clu_matching.dxCloPsdClu_track[1]);
        }
        
    }

void psd_stk_distributions(
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
        double psd_min_energy           {0};    // Minimum energy per PSD bar

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

                        }
                    }
                }
            }
        }
        
    }

