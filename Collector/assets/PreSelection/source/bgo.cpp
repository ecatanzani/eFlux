#include "bgo.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStruct.h"

#include "TMath.h"

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

void bgofiducial_distributions(
    std::shared_ptr<DmpEvtBgoHits> bgohits,
    std::shared_ptr<DmpEvtBgoRec> bgorec,
    std::shared_ptr<DmpEvtHeader> evt_header,
    const double evt_corr_energy_gev, 
    std::shared_ptr<histos> ps_histos) {

    double layer_min_energy {0}; //Minimum energy per layer
    std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
    bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), layer_min_energy);

    if (check_trigger(evt_header)) {
        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
            ps_histos->h_energy_fraction->Fill(elm_energy_fraction);
        
        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(56*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_56->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(60*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_60->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(61*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_61->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(62*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_62->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(63*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_63->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(64*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_64->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(65*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_65->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(70*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_70->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(80*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_80->Fill(elm_energy_fraction);

        if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(85*TMath::DegToRad()))
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_large_angles_85->Fill(elm_energy_fraction);

        ps_histos->h_bar_energy->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()));
        ps_histos->h_bars_last_layer_10MeV->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10));
        ps_histos->h_bars_last_layer_10MeV_2D->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10));
        ps_histos->h_maxrms->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
        ps_histos->h_maxrms_2D->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));

        std::vector<double> bgorec_proj_X, bgorec_proj_Y;
        std::tie(bgorec_proj_X, bgorec_proj_Y) = get_shorew_axis_reco_projections(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
        ps_histos->h_bgoshower_top_X->Fill(bgorec_proj_X[0]);
        ps_histos->h_bgoshower_bottom_X->Fill(bgorec_proj_X[1]);
        ps_histos->h_bgoshower_top_Y->Fill(bgorec_proj_Y[0]);
        ps_histos->h_bgoshower_bottom_Y->Fill(bgorec_proj_Y[1]);


        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
            if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY))
                ps_histos->h_energy_fraction_sh_axis_contained->Fill(elm_energy_fraction);
            else
                ps_histos->h_energy_fraction_sh_axis_not_contained->Fill(elm_energy_fraction);
    }
    else {
        if (evt_corr_energy_gev) {
            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                ps_histos->h_energy_fraction_no_trigger->Fill(elm_energy_fraction);

            ps_histos->h_maxrms_no_trigger->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
            ps_histos->h_maxrms_2D_no_trigger->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
        }
    }
}