#include "cuts.h"
#include "lateral_showering.h"

#include "Dmp/DmpBgoContainer.h"

inline bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
    return trigger_general;
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