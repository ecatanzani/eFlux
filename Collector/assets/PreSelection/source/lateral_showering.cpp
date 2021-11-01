#include "cuts.h"
#include "lateral_showering.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"
#include "Dmp/DmpStruct.h"

inline bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
    return trigger_general;
}

inline double get_mean_bar_energy(const std::vector<std::vector<double>> bar_energy) {
    double esum {0};
    double gev {0.001};
    double def_value = -999;
    for (auto&& layer_energy : bar_energy)
        for (auto && single_bar_energy : layer_energy)
            esum += single_bar_energy != def_value ? single_bar_energy : 0;
    return (esum/(bar_energy.size()*bar_energy[0].size()))*gev;
}

inline unsigned int count_bars_on_layer(const std::vector<double> layer_energy, const double energy_threshold) {
    unsigned int bars {0};
    for (auto&& energy : layer_energy)
        if (energy>=energy_threshold)
            bars += 1;
    return bars;
}

inline double get_max_rms(const std::vector<double> rms_layer, const std::vector<std::vector<short>> layerBarNumber, const double bgo_total_raw_energy) {
    auto max_rms = rms_layer[0];
    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx) {
        double layerTotEnergy = 0;
        std::accumulate(layerBarNumber[lIdx].begin(), layerBarNumber[lIdx].end(), layerTotEnergy);
        if (layerTotEnergy > bgo_total_raw_energy/100.)
            if (rms_layer[lIdx] > max_rms)
                max_rms = rms_layer[lIdx];
    }
    return max_rms;
}

void lateral_showering_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos) {

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();

        double bgo_layer_min_energy     {0};    // Minimum energy per BGO layer
        double bgo_max_energy_ratio     {0.35}; // Maximum energy ratio per layer
        double bgo_shower_axis_delta    {280};  // BGO maximum shower axis delta (mm)
        double bgo_shower_width         {100};  // BGO maximum shower width (mm)

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), bgo_layer_min_energy);

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

                        auto weight {ps_histos->GetWeight()};
                        
                        ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                        else if (evt_corr_energy_gev>=5000)
                            ps_histos->h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                    }
                }
            }
        }
    }

void lateral_showering_distributions_lastcut(
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

                        if (psdcharge_cut) {
                            auto stkcharge_cut = stk_charge_cut(stkclusters, event_best_track, STK_single_charge);

                            if (stkcharge_cut) {

                                auto weight {ps_histos->GetWeight()};
                                double gev {0.001};

                                ps_histos->h_bar_energy_lateral_showering_lastcut->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()), weight);

                                for (auto&& layer_energy : bgoVault->GetLayerBarEnergies())
                                    for (auto && single_bar_energy : layer_energy)
                                        ps_histos->h_bar_energy_2D_lateral_showering_lastcut->Fill(evt_corr_energy_gev, single_bar_energy*gev, weight);

                                ps_histos->h_bars_last_layer_10MeV_lateral_showering_lastcut->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);
                                ps_histos->h_bars_last_layer_10MeV_2D_lateral_showering_lastcut->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);

                                ps_histos->h_maxrms_lateral_showering_lastcut->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerBarNumber(), bgorec->GetTotalEnergy()), weight);
                                ps_histos->h_maxrms_2D_lateral_showering_lastcut->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerBarNumber(), bgorec->GetTotalEnergy()), weight);

                                ps_histos->h_BGOrec_sumRms_flast_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
                                    ps_histos->h_BGOrec_sumRms_flast_20_100_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
                                    ps_histos->h_BGOrec_sumRms_flast_100_250_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
                                    ps_histos->h_BGOrec_sumRms_flast_250_500_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
                                    ps_histos->h_BGOrec_sumRms_flast_500_1000_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
                                    ps_histos->h_BGOrec_sumRms_flast_1000_3000_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
                                    ps_histos->h_BGOrec_sumRms_flast_3000_5000_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                else if (evt_corr_energy_gev>=5000)
                                    ps_histos->h_BGOrec_sumRms_flast_5000_lateral_showering_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);

                            }
                        }
                    }
                }
            }
        }

    }