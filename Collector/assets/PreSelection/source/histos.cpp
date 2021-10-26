#include "histos.h"
#include "binning.h"

#include <iostream>

#include "Dmp/DmpStruct.h"

#include "TFile.h"

histos::histos(std::shared_ptr<energy_config> econfig, const bool mc) {
    
    h_simu = mc;
    energy_binning = econfig->GetEnergyBinning();
    energy_nbins = (int)energy_binning.size() - 1;

    bars_energy_bins = createLogBinning(1e-2, 1e+5, 1000);
    bars_2d_energy_bins = createLogBinning(1e-5, 1e+5, 2000);
    number_of_bars_last_layer = createLinearBinning(0, 22, 23);
    max_rms_bins = createLinearBinning(0, 3000, 1000);
    sumRms_binning = createLogBinning(10, 2e+3, 1e+2);
    flast_binning = createLogBinning(1e-5, 2e-1, 1e+3);

    h_energy_fraction = std::make_shared<TH1D>("h_energy_fraction", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_no_trigger = std::make_shared<TH1D>("h_energy_fraction_no_trigger", "Energy fraction - No triggered events; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_56 = std::make_shared<TH1D>("h_energy_fraction_large_angles_56", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_60 = std::make_shared<TH1D>("h_energy_fraction_large_angles_60", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_61 = std::make_shared<TH1D>("h_energy_fraction_large_angles_61", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_62 = std::make_shared<TH1D>("h_energy_fraction_large_angles_62", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_63 = std::make_shared<TH1D>("h_energy_fraction_large_angles_63", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_64 = std::make_shared<TH1D>("h_energy_fraction_large_angles_64", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_65 = std::make_shared<TH1D>("h_energy_fraction_large_angles_65", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_70 = std::make_shared<TH1D>("h_energy_fraction_large_angles_70", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_80 = std::make_shared<TH1D>("h_energy_fraction_large_angles_80", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_85 = std::make_shared<TH1D>("h_energy_fraction_large_angles_85", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);

    h_bar_energy = std::make_shared<TH2D>("h_bar_energy", "Mean energy bar; Reconstructed Energy [GeV]; Mean Bar Energy [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_energy_bins.size() -1, &bars_energy_bins[0]);
    h_bar_energy_2D = std::make_shared<TH2D>("h_bar_energy_2D", "Energy bar; Reconstructed Energy [GeV]; Energy Bar [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_2d_energy_bins.size() -1, &bars_2d_energy_bins[0]);
    h_bars_last_layer_10MeV = std::make_shared<TH1D>("h_bars_last_layer_10MeV", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 22, 0, 22);
    h_bars_last_layer_10MeV_2D = std::make_shared<TH2D>("h_bars_last_layer_10MeV_2D", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);
    h_maxrms = std::make_shared<TH1D>("h_maxrms", "Max RMS Layers", 1000, 0, 3000);
    h_maxrms_2D = std::make_shared<TH2D>("h_maxrms_2D", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);
    h_maxrms_no_trigger = std::make_shared<TH1D>("h_maxrms_no_trigger", "Max RMS Layers - No Trigger", 1000, 0, 3000);
    h_maxrms_2D_no_trigger = std::make_shared<TH2D>("h_maxrms_2D_no_trigger", "Max RMS - [> 1% energy content] - No Trigger; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

    h_bgoshower_top_X = std::make_shared<TH1D>("h_bgoshower_top_X", "BGO Shower X TOP Projection; Top X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_X = std::make_shared<TH1D>("h_bgoshower_bottom_X", "BGO Shower X BOTTOM Projection; BOTTOM X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_top_Y = std::make_shared<TH1D>("h_bgoshower_top_Y", "BGO Shower X TOP Projection; Top Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_Y = std::make_shared<TH1D>("h_bgoshower_bottom_Y", "BGO Shower Y BOTTOM Projection; BOTTOM Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    
    h_energy_fraction_sh_axis_contained = std::make_shared<TH1D>("h_energy_fraction_sh_axis_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_sh_axis_not_contained = std::make_shared<TH1D>("h_energy_fraction_sh_axis_not_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);

    h_PSD_STK_X_match_energy_int = std::make_shared<TH1D>("h_PSD_STK_X_match_energy_int", "PSD - STK match X view- All Energies; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_energy_int = std::make_shared<TH1D>("h_PSD_STK_Y_match_energy_int", "PSD - STK match Y view - All Energies; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_100_250 = std::make_shared<TH1D>("h_PSD_STK_X_match_100_250", "PSD - STK match X view - 100-250 GeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_100_250 = std::make_shared<TH1D>("h_PSD_STK_Y_match_100_250", "PSD - STK match Y view - 100-250 GeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_250_500 = std::make_shared<TH1D>("h_PSD_STK_X_match_250_500", "PSD - STK match X view - 250-500 GeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_250_500 = std::make_shared<TH1D>("h_PSD_STK_Y_match_250_500", "PSD - STK match Y view - 250-500 GeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_500_1000 = std::make_shared<TH1D>("h_PSD_STK_X_match_500_1000", "PSD - STK match X view - 500GeV-1TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_500_1000 = std::make_shared<TH1D>("h_PSD_STK_Y_match_500_1000", "PSD - STK match Y view - 500GeV-1TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_1000_5000 = std::make_shared<TH1D>("h_PSD_STK_X_match_1000_5000", "PSD - STK match X view - 1TeV-5TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_1000_5000 = std::make_shared<TH1D>("h_PSD_STK_Y_match_1000_5000", "PSD - STK match Y view - 1TeV-5TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_5000 = std::make_shared<TH1D>("h_PSD_STK_X_match_5000", "PSD - STK match X view - >5TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_5000 = std::make_shared<TH1D>("h_PSD_STK_Y_match_5000", "PSD - STK match Y view - >5TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);

    h_BGOrec_sumRms_flast = std::make_shared<TH2D>("h_BGOrec_sumRms_flast", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_20_100 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_20_100", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_100_250 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_100_250", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_250_500 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_250_500", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_500_1000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_500_1000", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_1000_3000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_1000_3000", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_3000_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_3000_5000", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_5000", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);

    h_BGOrec_sumRms_flast_after_bgofiducial = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_20_100 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_20_100", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_100_250 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_100_250", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_250_500 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_250_500", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_500_1000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_500_1000", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_bgofiducial_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_bgofiducial_5000", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);

    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);

    if (h_simu) {
        h_max_bar_position_simu_reco_energy_diff = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff", "Energy diff vs Macx Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13 = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_bgoshower_top_X_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_top_X_simu_reco_energy_diff", "Energy diff vs TOP X BGO position; TOP X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_X_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_bottom_X_simu_reco_energy_diff", "Energy diff vs BOTTOM X BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_top_Y_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_top_Y_simu_reco_energy_diff", "Energy diff vs TOP Y BGO position; TOP Y BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_Y_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_bottom_Y_simu_reco_energy_diff", "Energy diff vs BOTTOM Y BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);

        h_diff_bgo_simu_particle_direction_X = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_X", "X view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_particle_direction_Y = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_Y", "Y view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_extr_top_position_X = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_X", "TOP X ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
        h_diff_bgo_simu_extr_top_position_Y = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_Y", "TOP Y ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
    }
}

void histos::Write(const std::string output_wd, const bool verbose) {

    auto final_output_path = output_wd+std::string("/preselection.root");
    TFile* outfile = TFile::Open(final_output_path.c_str(), "RECREATE");
    if (!outfile->IsOpen()) {
        std::cerr << "\n\nError writing output file [" << final_output_path << "]\n\n";
        exit(100);
    }

    outfile->mkdir("BGO");
    outfile->cd("BGO");

    h_energy_fraction->Write();
    h_energy_fraction_no_trigger->Write();
    h_energy_fraction_large_angles_56->Write();
    h_energy_fraction_large_angles_60->Write();
    h_energy_fraction_large_angles_61->Write();
    h_energy_fraction_large_angles_62->Write();
    h_energy_fraction_large_angles_63->Write();
    h_energy_fraction_large_angles_64->Write();
    h_energy_fraction_large_angles_65->Write();
    h_energy_fraction_large_angles_70->Write();
    h_energy_fraction_large_angles_80->Write();
    h_energy_fraction_large_angles_85->Write();

    h_bar_energy->Write();
    h_bar_energy_2D->Write();
    h_bars_last_layer_10MeV->Write();
    h_bars_last_layer_10MeV_2D->Write();

    h_maxrms->Write();
    h_maxrms_2D->Write();
    h_maxrms_no_trigger->Write();
    h_maxrms_2D_no_trigger->Write();

    h_bgoshower_top_X->Write();
    h_bgoshower_bottom_X->Write();
    h_bgoshower_top_Y->Write();
    h_bgoshower_bottom_Y->Write();

    h_energy_fraction_sh_axis_contained->Write();
    h_energy_fraction_sh_axis_not_contained->Write();

    h_BGOrec_sumRms_flast->Write();
    h_BGOrec_sumRms_flast_20_100->Write();
    h_BGOrec_sumRms_flast_100_250->Write();
    h_BGOrec_sumRms_flast_250_500->Write();
    h_BGOrec_sumRms_flast_500_1000->Write();
    h_BGOrec_sumRms_flast_1000_3000->Write();
    h_BGOrec_sumRms_flast_3000_5000->Write();
    h_BGOrec_sumRms_flast_5000->Write();

    if (h_simu) {
        h_max_bar_position_simu_reco_energy_diff->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13->Write();

        h_bgoshower_top_X_simu_reco_energy_diff->Write();
        h_bgoshower_bottom_X_simu_reco_energy_diff->Write();
        h_bgoshower_top_Y_simu_reco_energy_diff->Write();
        h_bgoshower_bottom_Y_simu_reco_energy_diff->Write();

        h_diff_bgo_simu_particle_direction_X->Write();
        h_diff_bgo_simu_particle_direction_Y->Write();
        h_diff_bgo_simu_extr_top_position_X->Write();
        h_diff_bgo_simu_extr_top_position_Y->Write();
    }

    outfile->mkdir("BGO_fiducial");
    outfile->cd("BGO_fiducial");

    h_BGOrec_sumRms_flast_after_bgofiducial->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_20_100->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_100_250->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_250_500->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_500_1000->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_5000->Write();

    outfile->mkdir("lateral_showering");
    outfile->cd("lateral_showering");
    
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000->Write();
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000->Write();

    outfile->mkdir("PSD_STK");
    outfile->cd("PSD_STK");

    h_PSD_STK_X_match_energy_int->Write();
    h_PSD_STK_Y_match_energy_int->Write();
    h_PSD_STK_X_match_100_250->Write();
    h_PSD_STK_Y_match_100_250->Write();
    h_PSD_STK_X_match_250_500->Write();
    h_PSD_STK_Y_match_250_500->Write();
    h_PSD_STK_X_match_500_1000->Write();
    h_PSD_STK_Y_match_500_1000->Write();
    h_PSD_STK_X_match_1000_5000->Write();
    h_PSD_STK_Y_match_1000_5000->Write();
    h_PSD_STK_X_match_5000->Write();
    h_PSD_STK_Y_match_5000->Write();

    outfile->Close();

    if (verbose) std::cout << "\n\nOutput file has been written [" << final_output_path << "]\n\n";  
}