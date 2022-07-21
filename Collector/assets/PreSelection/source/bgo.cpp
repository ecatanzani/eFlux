#include "bgo.h"
#include "cuts.h"
#include "trigger.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStkContainer.h"
#include "Dmp/DmpPsdContainer.h"
#include "Dmp/DmpGeoStruct.h"

#include "TMath.h"
#include "TVector3.h"

void bgo_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    std::shared_ptr<DmpEvtHeader> evt_header,
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos,
    std::shared_ptr<config> cuts_config) {

    double gev {0.001};
    auto weight {ps_histos->GetWeight()}; // Weight for histos

    std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
    bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);

    if (check_trigger(evt_header)) {
        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer()) {
            ps_histos->h_energy_fraction->Fill(elm_energy_fraction, weight);
            ps_histos->h_energy_fraction_2D->Fill(evt_corr_energy_gev, elm_energy_fraction, weight);
        }

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(56*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_56->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(60*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_60->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(61*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_61->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(62*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_62->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(63*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_63->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(64*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_64->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(65*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_65->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(70*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_70->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(80*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_80->Fill(elm_energy_fraction, weight);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(85*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_85->Fill(elm_energy_fraction, weight);

        ps_histos->h_bar_energy->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()), weight);

        for (auto&& layer_energy : bgoVault->GetLayerBarEnergies())
            for (auto && single_bar_energy : layer_energy)
                ps_histos->h_bar_energy_2D->Fill(evt_corr_energy_gev, single_bar_energy*gev, weight);

        ps_histos->h_bars_last_layer_10MeV->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);
        ps_histos->h_bars_last_layer_10MeV_2D->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);

        ps_histos->h_maxrms->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
        ps_histos->h_maxrms_2D->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);

        std::vector<double> bgorec_proj_X, bgorec_proj_Y;
        std::tie(bgorec_proj_X, bgorec_proj_Y) = get_shorew_axis_reco_projections(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
        ps_histos->h_bgoshower_top_X->Fill(bgorec_proj_X[0], weight);
        ps_histos->h_bgoshower_bottom_X->Fill(bgorec_proj_X[1], weight);
        ps_histos->h_bgoshower_top_Y->Fill(bgorec_proj_Y[0], weight);
        ps_histos->h_bgoshower_bottom_Y->Fill(bgorec_proj_Y[1], weight);

        if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY)) {
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_sh_axis_contained->Fill(elm_energy_fraction, weight);
        }
        else {
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_sh_axis_not_contained->Fill(elm_energy_fraction, weight);
        }       
        
        ps_histos->h_BGOrec_sumRms_flast->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
            ps_histos->h_BGOrec_sumRms_flast_20_100->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
            ps_histos->h_BGOrec_sumRms_flast_100_250->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
            ps_histos->h_BGOrec_sumRms_flast_250_500->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
            ps_histos->h_BGOrec_sumRms_flast_500_1000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
            ps_histos->h_BGOrec_sumRms_flast_1000_3000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
            ps_histos->h_BGOrec_sumRms_flast_3000_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
        else if (evt_corr_energy_gev>=5000)
            ps_histos->h_BGOrec_sumRms_flast_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);

        if (simu_primaries!=nullptr) {

            auto energy_diff = (simu_primaries->pvpart_ekin - evt_corr_energy)/simu_primaries->pvpart_ekin;
            auto idxBarMaxLayer = bgoVault->GetIdxBarMaxLayer();

            ps_histos->h_bgoshower_top_X_simu_reco_energy_diff->Fill(bgorec_proj_X[0], energy_diff, weight);
            ps_histos->h_bgoshower_bottom_X_simu_reco_energy_diff->Fill(bgorec_proj_X[1], energy_diff, weight);
            ps_histos->h_bgoshower_top_Y_simu_reco_energy_diff->Fill(bgorec_proj_Y[0], energy_diff, weight);
            ps_histos->h_bgoshower_bottom_Y_simu_reco_energy_diff->Fill(bgorec_proj_Y[1], energy_diff, weight);

            for (int lIdx=0; lIdx<14; ++lIdx) {
                if (bgoVault->GetLayerBarNumber()[lIdx].size()) {
                    if (bgoVault->GetiMaxLayer()[lIdx] > -1) {
                        ps_histos->h_max_bar_position_simu_reco_energy_diff->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                        switch (lIdx) {
                            case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                            case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    break;
                        }
                    
                        if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY)) {
                    
                            ps_histos->h_max_bar_position_simu_reco_energy_diff_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                            switch (lIdx) {
                                case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                            }

                        }
                        else {

                            ps_histos->h_max_bar_position_simu_reco_energy_diff_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                            switch (lIdx) {
                                case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                        break;
                                        
                            }
                        }
                    }
                }
            }

            // Get info on BGO angular and spacial resolution (on the extrapolation)
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

            ps_histos->h_diff_bgo_simu_particle_direction_X->Fill(slope[0]-(bgoVault->GetBGOslope())[0], weight);
            ps_histos->h_diff_bgo_simu_particle_direction_Y->Fill(slope[1]-(bgoVault->GetBGOslope())[1], weight);
            ps_histos->h_diff_bgo_simu_extr_top_position_X->Fill(intercept[0]-(bgoVault->GetBGOintercept())[0], weight);
            ps_histos->h_diff_bgo_simu_extr_top_position_Y->Fill(intercept[1]-(bgoVault->GetBGOintercept())[1], weight);
        }
    }
    else {
        if (evt_corr_energy_gev) {
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_no_trigger->Fill(elm_energy_fraction, weight);

            ps_histos->h_maxrms_no_trigger->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
            ps_histos->h_maxrms_2D_no_trigger->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
        }
    }
}

void bgofiducial_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    std::shared_ptr<DmpEvtHeader> evt_header,
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos,
    std::shared_ptr<config> cuts_config) {

        double gev {0.001};
        auto weight {ps_histos->GetWeight()}; // Weight for histos

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        
        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);

        if (check_trigger(evt_header)) {

            auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
            auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
            auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);

            auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

            if (bgofiducial_cut) {

                for(auto&& elm_energy_fraction : bgoVault->GetFracLayer()) {
                    ps_histos->h_energy_fraction_after_bgofiducial->Fill(elm_energy_fraction, weight);
                    ps_histos->h_energy_fraction_2D_after_bgofiducial->Fill(evt_corr_energy_gev, elm_energy_fraction, weight);
                }
                
                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(56*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_56_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(60*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_60_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(61*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_61_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(62*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_62_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(63*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_63_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(64*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_64_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(65*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_65_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(70*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_70_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(80*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_80_after_bgofiducial->Fill(elm_energy_fraction, weight);

                if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(85*TMath::DegToRad()))
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_large_angles_85_after_bgofiducial->Fill(elm_energy_fraction, weight);

                ps_histos->h_bar_energy_after_bgofiducial->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()), weight);

                for (auto&& layer_energy : bgoVault->GetLayerBarEnergies())
                    for (auto && single_bar_energy : layer_energy)
                        ps_histos->h_bar_energy_2D_after_bgofiducial->Fill(evt_corr_energy_gev, single_bar_energy*gev, weight);

                    
                ps_histos->h_bars_last_layer_10MeV_after_bgofiducial->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);
                ps_histos->h_bars_last_layer_10MeV_2D_after_bgofiducial->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);
                ps_histos->h_maxrms_after_bgofiducial->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
                ps_histos->h_maxrms_2D_after_bgofiducial->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);

                std::vector<double> bgorec_proj_X, bgorec_proj_Y;
                std::tie(bgorec_proj_X, bgorec_proj_Y) = get_shorew_axis_reco_projections(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
                ps_histos->h_bgoshower_top_X_after_bgofiducial->Fill(bgorec_proj_X[0], weight);
                ps_histos->h_bgoshower_bottom_X_after_bgofiducial->Fill(bgorec_proj_X[1], weight);
                ps_histos->h_bgoshower_top_Y_after_bgofiducial->Fill(bgorec_proj_Y[0]);
                ps_histos->h_bgoshower_bottom_Y_after_bgofiducial->Fill(bgorec_proj_Y[1], weight);

                if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY)) {
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_sh_axis_contained_after_bgofiducial->Fill(elm_energy_fraction, weight);
                }
                else {
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        ps_histos->h_energy_fraction_sh_axis_not_contained_after_bgofiducial->Fill(elm_energy_fraction, weight);
                }
                
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_20_100->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_100_250->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_250_500->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_500_1000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                else if (evt_corr_energy_gev>=5000)
                    ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);

                if (simu_primaries!=nullptr) {

                    auto energy_diff = (simu_primaries->pvpart_ekin - evt_corr_energy)/simu_primaries->pvpart_ekin;
                    auto idxBarMaxLayer = bgoVault->GetIdxBarMaxLayer();

                    ps_histos->h_bgoshower_top_X_simu_reco_energy_diff_after_bgofiducial->Fill(bgorec_proj_X[0], energy_diff, weight);
                    ps_histos->h_bgoshower_bottom_X_simu_reco_energy_diff_after_bgofiducial->Fill(bgorec_proj_X[1], energy_diff, weight);
                    ps_histos->h_bgoshower_top_Y_simu_reco_energy_diff_after_bgofiducial->Fill(bgorec_proj_Y[0], energy_diff, weight);
                    ps_histos->h_bgoshower_bottom_Y_simu_reco_energy_diff_after_bgofiducial->Fill(bgorec_proj_Y[1], energy_diff, weight);

                    for (int lIdx=0; lIdx<14; ++lIdx) {
                        if (bgoVault->GetLayerBarNumber()[lIdx].size()) {
                            if (bgoVault->GetiMaxLayer()[lIdx] > -1) {
                                ps_histos->h_max_bar_position_simu_reco_energy_diff_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                switch (lIdx) {
                                    case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                    case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                            break;
                                }

                                if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY)) {

                                    ps_histos->h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    switch (lIdx) {
                                        case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                    }

                                }
                                else {

                                    ps_histos->h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                    switch (lIdx) {
                                        case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                        case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                break;
                                    }

                                }
                            }
                        }
                    }

                    // Get info on BGO angular and spacial resolution (on the extrapolation)
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

                    ps_histos->h_diff_bgo_simu_particle_direction_X_after_bgofiducial->Fill(slope[0]-(bgoVault->GetBGOslope())[0], weight);
                    ps_histos->h_diff_bgo_simu_particle_direction_Y_after_bgofiducial->Fill(slope[1]-(bgoVault->GetBGOslope())[1], weight);
                    ps_histos->h_diff_bgo_simu_extr_top_position_X_after_bgofiducial->Fill(intercept[0]-(bgoVault->GetBGOintercept())[0], weight);
                    ps_histos->h_diff_bgo_simu_extr_top_position_Y_after_bgofiducial->Fill(intercept[1]-(bgoVault->GetBGOintercept())[1], weight);
                }
            }
        }
        else {
            if (evt_corr_energy_gev) {
                for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                    ps_histos->h_energy_fraction_no_trigger->Fill(elm_energy_fraction, weight);

                ps_histos->h_maxrms_no_trigger->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
                ps_histos->h_maxrms_2D_no_trigger->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
            }
        }
    }

void bgofiducial_distributions_lastcut(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
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
        stkVault->scanSTKHits(stkclusters, stktracks, bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
        psdVault->scanPSDHits(psdhits, (cuts_config->GetCutsConfig()).psd_min_energy);

        if (check_trigger(evt_header)) {
            
            auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarIndex(13), evt_energy);
            
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

                            if (psdcharge_cut) {
                                auto stkcharge_cut = stk_charge_cut(stkclusters, event_best_track, (cuts_config->GetCutsConfig()).STK_single_charge);

                                if (stkcharge_cut) {

                                    double gev {0.001};
                                    auto weight {ps_histos->GetWeight()};
                                    
                                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer()) {
                                        ps_histos->h_energy_fraction_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);
                                        ps_histos->h_energy_fraction_2D_bgofiducial_lastcut->Fill(evt_corr_energy_gev, elm_energy_fraction, weight);
                                    }

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(56*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_56_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(60*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_60_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(61*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_61_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(62*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_62_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(63*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_63_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(64*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_64_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(65*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_65_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(70*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_70_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(80*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_80_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(85*TMath::DegToRad()))
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_large_angles_85_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);

                                    ps_histos->h_bar_energy_bgofiducial_lastcut->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()), weight);

                                    for (auto&& layer_energy : bgoVault->GetLayerBarEnergies())
                                        for (auto && single_bar_energy : layer_energy)
                                            ps_histos->h_bar_energy_2D_bgofiducial_lastcut->Fill(evt_corr_energy_gev, single_bar_energy*gev, weight);

                                    ps_histos->h_bars_last_layer_10MeV_bgofiducial_lastcut->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);
                                    ps_histos->h_bars_last_layer_10MeV_2D_bgofiducial_lastcut->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10), weight);

                                    ps_histos->h_maxrms_bgofiducial_lastcut->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);
                                    ps_histos->h_maxrms_2D_bgofiducial_lastcut->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetELayer(), bgorec->GetTotalEnergy()), weight);

                                    std::vector<double> bgorec_proj_X, bgorec_proj_Y;
                                    std::tie(bgorec_proj_X, bgorec_proj_Y) = get_shorew_axis_reco_projections(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
                                    ps_histos->h_bgoshower_top_X_bgofiducial_lastcut->Fill(bgorec_proj_X[0], weight);
                                    ps_histos->h_bgoshower_bottom_X_bgofiducial_lastcut->Fill(bgorec_proj_X[1], weight);
                                    ps_histos->h_bgoshower_top_Y_bgofiducial_lastcut->Fill(bgorec_proj_Y[0], weight);
                                    ps_histos->h_bgoshower_bottom_Y_bgofiducial_lastcut->Fill(bgorec_proj_Y[1], weight);

                                    if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY)) {
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_sh_axis_contained_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);
                                    }
                                    else {
                                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                            ps_histos->h_energy_fraction_sh_axis_not_contained_bgofiducial_lastcut->Fill(elm_energy_fraction, weight);
                                    }

                                    ps_histos->h_BGOrec_sumRms_flast_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
                                        ps_histos->h_BGOrec_sumRms_flast_20_100_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
                                        ps_histos->h_BGOrec_sumRms_flast_100_250_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
                                        ps_histos->h_BGOrec_sumRms_flast_250_500_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
                                        ps_histos->h_BGOrec_sumRms_flast_500_1000_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
                                        ps_histos->h_BGOrec_sumRms_flast_1000_3000_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
                                        ps_histos->h_BGOrec_sumRms_flast_3000_5000_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);
                                    else if (evt_corr_energy_gev>=5000)
                                        ps_histos->h_BGOrec_sumRms_flast_5000_bgofiducial_lastcut->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13), weight);

                                    if (simu_primaries!=nullptr) {

                                        auto energy_diff = (simu_primaries->pvpart_ekin - evt_corr_energy)/simu_primaries->pvpart_ekin;
                                        auto idxBarMaxLayer = bgoVault->GetIdxBarMaxLayer();

                                        ps_histos->h_bgoshower_top_X_simu_reco_energy_diff_bgofiducial_lastcut->Fill(bgorec_proj_X[0], energy_diff, weight);
                                        ps_histos->h_bgoshower_bottom_X_simu_reco_energy_diff_bgofiducial_lastcut->Fill(bgorec_proj_X[1], energy_diff, weight);
                                        ps_histos->h_bgoshower_top_Y_simu_reco_energy_diff_bgofiducial_lastcut->Fill(bgorec_proj_Y[0], energy_diff, weight);
                                        ps_histos->h_bgoshower_bottom_Y_simu_reco_energy_diff_bgofiducial_lastcut->Fill(bgorec_proj_Y[1], energy_diff, weight);

                                        for (int lIdx=0; lIdx<14; ++lIdx) {
                                            if (bgoVault->GetLayerBarNumber()[lIdx].size()) {
                                                if (bgoVault->GetiMaxLayer()[lIdx] > -1) {
                                                    ps_histos->h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                    switch (lIdx) {
                                                        case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                        case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                break;
                                                    }

                                                    if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY)) {
                                                        
                                                        ps_histos->h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                        switch (lIdx) {
                                                            case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_in->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                        }

                                                    }
                                                    else {
                                                        
                                                        ps_histos->h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                        switch (lIdx) {
                                                            case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                            case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_out->Fill(idxBarMaxLayer[lIdx], energy_diff, weight);
                                                                    break;
                                                        }

                                                    }
                                                }
                                            }
                                        }

                                        // Get info on BGO angular and spacial resolution (on the extrapolation)
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

                                        ps_histos->h_diff_bgo_simu_particle_direction_X_bgofiducial_lastcut->Fill(slope[0]-(bgoVault->GetBGOslope())[0], weight);
                                        ps_histos->h_diff_bgo_simu_particle_direction_Y_bgofiducial_lastcut->Fill(slope[1]-(bgoVault->GetBGOslope())[1], weight);
                                        ps_histos->h_diff_bgo_simu_extr_top_position_X_bgofiducial_lastcut->Fill(intercept[0]-(bgoVault->GetBGOintercept())[0], weight);
                                        ps_histos->h_diff_bgo_simu_extr_top_position_Y_bgofiducial_lastcut->Fill(intercept[1]-(bgoVault->GetBGOintercept())[1], weight);
                                    }

                                    // Fill efficiency histos
                                    auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);
                                    auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
                                    auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

double get_mean_bar_energy(const std::vector<std::vector<double>> bar_energy) {
    double esum {0};
    double gev {0.001};
    double def_value = -999;
    for (auto&& layer_energy : bar_energy)
        for (auto && single_bar_energy : layer_energy)
            esum += single_bar_energy != def_value ? single_bar_energy : 0;
    return (esum/(bar_energy.size()*bar_energy[0].size()))*gev;
}

unsigned int count_bars_on_layer(const std::vector<double> layer_energy, const double energy_threshold) {
        unsigned int bars {0};
        for (auto&& energy : layer_energy)
        if (energy>=energy_threshold)
            bars += 1;
        return bars;
}

double get_max_rms(const std::vector<double> rms_layer, const std::vector<double> layer_energy, const double bgo_total_raw_energy) {
    double max_rms = rms_layer[0];
    double e_threshold = bgo_total_raw_energy/100.;
    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx) {
        if (layer_energy[lIdx] > e_threshold)
            if (rms_layer[lIdx] > max_rms)
                max_rms = rms_layer[lIdx];
    }
    return max_rms;
}

std::tuple<std::vector<double>, std::vector<double>> get_shorew_axis_reco_projections(const std::vector<double> brorec_slope, const std::vector<double> bgorec_intercept) {
        double topX = brorec_slope[0] * BGO_TopZ + bgorec_intercept[0];
        double topY = brorec_slope[1] * BGO_TopZ + bgorec_intercept[1];
        double bottomX = brorec_slope[0] * BGO_BottomZ + bgorec_intercept[0];
        double bottomY = brorec_slope[1] * BGO_BottomZ + bgorec_intercept[1];
        return std::tuple<std::vector<double>, std::vector<double>>(std::vector<double>{topX, bottomX}, std::vector<double>{topY, bottomY});
}

const bool isEventHittingEdgeBars(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header, 
    std::shared_ptr<TClonesArray> stkclusters, 
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    const double evt_energy,
    std::shared_ptr<config> cuts_config) {
        
        bool is_hitting_edge_bars {false};

        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        std::unique_ptr<DmpStkContainer> stkVault = std::make_unique<DmpStkContainer>();
        std::unique_ptr<DmpPsdContainer> psdVault = std::make_unique<DmpPsdContainer>();

        best_track event_best_track;
        psd_cluster_match clu_matching;

        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);
        stkVault->scanSTKHits(stkclusters, stktracks, bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
        psdVault->scanPSDHits(psdhits, (cuts_config->GetCutsConfig()).psd_min_energy);

        if (check_trigger(evt_header)) {
            
            auto nbarlayer13_cut = nBarLayer13_cut(bgohits, bgoVault->GetSingleLayerBarIndex(13), evt_energy);
            
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
                                is_hitting_edge_bars = ! maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
                        }
                    }
                }
            }
        }

        return is_hitting_edge_bars;
    }

const bool isEventOutsideBGOFiducial(
    std::shared_ptr<DmpEvtBgoHits> bgohits, 
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header,
    const double evt_energy,
    std::shared_ptr<config> cuts_config) {

        bool is_outside_bgo_fiducial {false};
        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);

        if (check_trigger(evt_header))
            is_outside_bgo_fiducial = ! BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), (cuts_config->GetCutsConfig()).bgo_shower_axis_delta);

        return is_outside_bgo_fiducial;
    
    }

const bool isEnergyPercentageGT35(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    std::shared_ptr<DmpEvtBgoRec> bgorec, 
    std::shared_ptr<DmpEvtHeader> evt_header,
    std::shared_ptr<TClonesArray> stkclusters,
    std::shared_ptr<TClonesArray> stktracks,
    std::shared_ptr<DmpEvtPsdHits> psdhits,
    const double evt_energy,
    std::shared_ptr<config> cuts_config) {

        bool is_energy_percentage_gt_35 {false};
        std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
        bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), (cuts_config->GetCutsConfig()).bgo_layer_min_energy);
        
        if (check_trigger(evt_header))
            is_energy_percentage_gt_35 = ! maxElayer_cut(bgoVault->GetLayerEnergies(), (cuts_config->GetCutsConfig()).bgo_max_energy_ratio, evt_energy);

        return is_energy_percentage_gt_35;
    }