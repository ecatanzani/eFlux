#include "bgo.h"
#include "cuts.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStruct.h"

#include "TMath.h"
#include "TVector3.h"

#include <vector>
#include <tuple>

inline bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
    return trigger_general;
}

inline bool check_HE_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	return evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
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

inline double get_max_rms(const std::vector<double> rms_layer, const std::vector<double> energy_layer, const double bgo_total_raw_energy) {
    double maxrms = -999;
    for (auto it=std::begin(energy_layer); it!=std::end(energy_layer); ++it) {
        if (*it>bgo_total_raw_energy/100) {
            if (maxrms == -999) maxrms = rms_layer[std::distance(std::begin(energy_layer), it)];
            else {
                if (rms_layer[std::distance(std::begin(energy_layer), it)] > maxrms)
                    maxrms = rms_layer[std::distance(std::begin(energy_layer), it)];
            }
        }
    }
    return maxrms;
}

inline std::tuple<std::vector<double>, std::vector<double>> get_shorew_axis_reco_projections(const std::vector<double> brorec_slope, const std::vector<double> bgorec_intercept) {
    double topX = brorec_slope[0] * BGO_TopZ + bgorec_intercept[0];
	double topY = brorec_slope[1] * BGO_TopZ + bgorec_intercept[1];
	double bottomX = brorec_slope[0] * BGO_BottomZ + bgorec_intercept[0];
	double bottomY = brorec_slope[1] * BGO_BottomZ + bgorec_intercept[1];
    return std::tuple<std::vector<double>, std::vector<double>>(std::vector<double>{topX, bottomX}, std::vector<double>{topY, bottomY});
}

void bgo_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    std::shared_ptr<DmpEvtHeader> evt_header,
    std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const double evt_energy, 
    const double evt_corr_energy,
    const double evt_energy_gev, 
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos) {

    double gev {0.001};
    double layer_min_energy {0}; //Minimum energy per layer
    auto weight {ps_histos->GetWeight()}; // Weight for histos

    std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
    bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), layer_min_energy);

    if (check_trigger(evt_header)) {
        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
            ps_histos->h_energy_fraction->Fill(elm_energy_fraction, weight);
        
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

        ps_histos->h_bar_energy->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()));

        for (auto&& layer_energy : bgoVault->GetLayerBarEnergies())
            for (auto && single_bar_energy : layer_energy)
                ps_histos->h_bar_energy_2D->Fill(evt_corr_energy_gev, single_bar_energy*gev);

        ps_histos->h_bars_last_layer_10MeV->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10));
        ps_histos->h_bars_last_layer_10MeV_2D->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10));
        ps_histos->h_maxrms->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()), weight);
        ps_histos->h_maxrms_2D->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()), weight);

        std::vector<double> bgorec_proj_X, bgorec_proj_Y;
        std::tie(bgorec_proj_X, bgorec_proj_Y) = get_shorew_axis_reco_projections(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
        ps_histos->h_bgoshower_top_X->Fill(bgorec_proj_X[0]);
        ps_histos->h_bgoshower_bottom_X->Fill(bgorec_proj_X[1]);
        ps_histos->h_bgoshower_top_Y->Fill(bgorec_proj_Y[0]);
        ps_histos->h_bgoshower_bottom_Y->Fill(bgorec_proj_Y[1]);

        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
            if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY))
                ps_histos->h_energy_fraction_sh_axis_contained->Fill(elm_energy_fraction, weight);
            else
                ps_histos->h_energy_fraction_sh_axis_not_contained->Fill(elm_energy_fraction, weight);
        
        ps_histos->h_BGOrec_sumRms_flast->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
            ps_histos->h_BGOrec_sumRms_flast_20_100->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
            ps_histos->h_BGOrec_sumRms_flast_100_250->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
            ps_histos->h_BGOrec_sumRms_flast_250_500->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
            ps_histos->h_BGOrec_sumRms_flast_500_1000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
            ps_histos->h_BGOrec_sumRms_flast_1000_3000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
            ps_histos->h_BGOrec_sumRms_flast_3000_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        else if (evt_corr_energy_gev>=5000)
            ps_histos->h_BGOrec_sumRms_flast_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));

        double bgo_layer_min_energy     {0};    // Minimum energy per BGO layer
        double bgo_max_energy_ratio     {0.35}; // Maximum energy ratio per layer
        double bgo_shower_axis_delta    {280};  // BGO maximum shower axis delta (mm)
        double bgo_shower_width         {100};  // BGO maximum shower width (mm)

        auto maxelayer_cut = maxElayer_cut(bgoVault->GetLayerEnergies(), bgo_max_energy_ratio, evt_energy);
        auto maxbarlayer_cut = maxBarLayer_cut(bgoVault->GetLayerBarNumber(), bgoVault->GetiMaxLayer(), bgoVault->GetIdxBarMaxLayer());
        auto bgotrack_cut = BGOTrackContainment_cut(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept(), bgo_shower_axis_delta);

        auto bgofiducial_cut = maxelayer_cut && maxbarlayer_cut && bgotrack_cut;

        if (bgofiducial_cut) {
            ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            if (evt_corr_energy_gev>=20 && evt_corr_energy_gev<100)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_20_100->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            else if (evt_corr_energy_gev>=100 && evt_corr_energy_gev<250)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_100_250->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            else if (evt_corr_energy_gev>=250 && evt_corr_energy_gev<500)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_250_500->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            else if (evt_corr_energy_gev>=500 && evt_corr_energy_gev<1000)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_500_1000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            else if (evt_corr_energy_gev>=1000 && evt_corr_energy_gev<3000)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            else if (evt_corr_energy_gev>=3000 && evt_corr_energy_gev<5000)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
            else if (evt_corr_energy_gev>=5000)
                ps_histos->h_BGOrec_sumRms_flast_after_bgofiducial_5000->Fill(bgoVault->GetSumRMS(), bgoVault->GetSingleFracLayer(13));
        }

        if (simu_primaries!=nullptr) {

            auto energy_diff = (simu_primaries->pvpart_ekin - evt_corr_energy)/simu_primaries->pvpart_ekin;
            auto idxBarMaxLayer = bgoVault->GetIdxBarMaxLayer();

            ps_histos->h_bgoshower_top_X_simu_reco_energy_diff->Fill(bgorec_proj_X[0], energy_diff);
            ps_histos->h_bgoshower_bottom_X_simu_reco_energy_diff->Fill(bgorec_proj_X[1], energy_diff);
            ps_histos->h_bgoshower_top_Y_simu_reco_energy_diff->Fill(bgorec_proj_Y[0], energy_diff);
            ps_histos->h_bgoshower_bottom_Y_simu_reco_energy_diff->Fill(bgorec_proj_Y[1], energy_diff);

            for (int lIdx=0; lIdx<14; ++lIdx) {
                if (bgoVault->GetLayerBarNumber()[lIdx].size()) {
                    if (bgoVault->GetiMaxLayer()[lIdx] > -1) {
                        ps_histos->h_max_bar_position_simu_reco_energy_diff->Fill(idxBarMaxLayer[lIdx], energy_diff);
                        switch (lIdx) {
                            case 0: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_0->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 1: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_1->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 2: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_2->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 3: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_3->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 4: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_4->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 5: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_5->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 6: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_6->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 7: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_7->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 8: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_8->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 9: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_9->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 10: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_10->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 11: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_11->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 12: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_12->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
                            case 13: ps_histos->h_max_bar_position_simu_reco_energy_diff_ly_13->Fill(idxBarMaxLayer[lIdx], energy_diff);
                                    break;
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

            ps_histos->h_diff_bgo_simu_particle_direction_X->Fill(slope[0]-(bgoVault->GetBGOslope())[0]);
            ps_histos->h_diff_bgo_simu_particle_direction_Y->Fill(slope[1]-(bgoVault->GetBGOslope())[1]);
            ps_histos->h_diff_bgo_simu_extr_top_position_X->Fill(intercept[0]-(bgoVault->GetBGOintercept())[0]);
            ps_histos->h_diff_bgo_simu_extr_top_position_Y->Fill(intercept[1]-(bgoVault->GetBGOintercept())[1]);
        }
    }
    else {
        if (evt_corr_energy_gev) {
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_no_trigger->Fill(elm_energy_fraction, weight);

            ps_histos->h_maxrms_no_trigger->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
            ps_histos->h_maxrms_2D_no_trigger->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
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
    std::shared_ptr<histos> ps_histos) {



    }