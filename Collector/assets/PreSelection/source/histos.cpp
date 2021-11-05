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
    energy_fraction_bins = createLinearBinning(0, 1, 100);
    bars_2d_energy_bins = createLogBinning(1e-5, 1e+5, 2000);
    number_of_bars_last_layer = createLinearBinning(0, 22, 22);
    max_rms_bins = createLinearBinning(0, 3000, 1000);
    sumRms_binning = createLogBinning(10, 2e+3, 1e+2);
    flast_binning = createLogBinning(1e-5, 2e-1, 1e+3);
    stk_track_binning = createLinearBinning(0, 10, 10);
    stk_ntracks_binning = createLinearBinning(0, 1000, 1000);
    psd_clusters_binning = createLinearBinning(0, 50, 50);

    h_energy_fraction = std::make_shared<TH1D>("h_energy_fraction", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_2D = std::make_shared<TH2D>("h_energy_fraction_2D", "Energy Fraction; Reconstructed Energy [GeV]; Energy fraction", energy_nbins -1, &energy_binning[0], (int)energy_fraction_bins.size() -1, &energy_fraction_bins[0]); 
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

    h_energy_fraction_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_after_bgofiducial", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_2D_after_bgofiducial = std::make_shared<TH2D>("h_energy_fraction_2D_after_bgofiducial", "Energy Fraction; Reconstructed Energy [GeV]; Energy fraction", energy_nbins -1, &energy_binning[0], (int)energy_fraction_bins.size() -1, &energy_fraction_bins[0]);
    h_energy_fraction_no_trigger_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_no_trigger_after_bgofiducial", "Energy fraction - No triggered events; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_56_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_56_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_60_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_60_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_61_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_61_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_62_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_62_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_63_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_63_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_64_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_64_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_65_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_65_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_70_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_70_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_80_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_80_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_85_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_large_angles_85_after_bgofiducial", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);

    h_energy_fraction_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_bgofiducial_lastcut", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_2D_bgofiducial_lastcut = std::make_shared<TH2D>("h_energy_fraction_2D_bgofiducial_lastcut", "Energy Fraction; Reconstructed Energy [GeV]; Energy fraction", energy_nbins -1, &energy_binning[0], (int)energy_fraction_bins.size() -1, &energy_fraction_bins[0]);
    h_energy_fraction_no_trigger_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_no_trigger_bgofiducial_lastcut", "Energy fraction - No triggered events; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_56_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_56_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_60_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_60_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_61_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_61_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_62_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_62_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_63_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_63_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_64_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_64_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_65_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_65_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_70_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_70_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_80_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_80_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_large_angles_85_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_large_angles_85_bgofiducial_lastcut", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);

    h_bar_energy = std::make_shared<TH2D>("h_bar_energy", "Mean energy bar; Reconstructed Energy [GeV]; Mean Bar Energy [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_energy_bins.size() -1, &bars_energy_bins[0]);
    h_bar_energy_2D = std::make_shared<TH2D>("h_bar_energy_2D", "Energy bar; Reconstructed Energy [GeV]; Energy Bar [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_2d_energy_bins.size() -1, &bars_2d_energy_bins[0]);
    h_bars_last_layer_10MeV = std::make_shared<TH1D>("h_bars_last_layer_10MeV", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 22, 0, 22);
    h_bars_last_layer_10MeV_2D = std::make_shared<TH2D>("h_bars_last_layer_10MeV_2D", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);

    h_bar_energy_after_bgofiducial = std::make_shared<TH2D>("h_bar_energy_after_bgofiducial", "Mean energy bar; Reconstructed Energy [GeV]; Mean Bar Energy [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_energy_bins.size() -1, &bars_energy_bins[0]);
    h_bar_energy_2D_after_bgofiducial = std::make_shared<TH2D>("h_bar_energy_2D_after_bgofiducial", "Energy bar; Reconstructed Energy [GeV]; Energy Bar [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_2d_energy_bins.size() -1, &bars_2d_energy_bins[0]);
    h_bars_last_layer_10MeV_after_bgofiducial = std::make_shared<TH1D>("h_bars_last_layer_10MeV_after_bgofiducial", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 22, 0, 22);
    h_bars_last_layer_10MeV_2D_after_bgofiducial = std::make_shared<TH2D>("h_bars_last_layer_10MeV_2D_after_bgofiducial", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);

    h_bar_energy_bgofiducial_lastcut = std::make_shared<TH2D>("h_bar_energy_bgofiducial_lastcut", "Mean energy bar; Reconstructed Energy [GeV]; Mean Bar Energy [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_energy_bins.size() -1, &bars_energy_bins[0]);
    h_bar_energy_2D_bgofiducial_lastcut = std::make_shared<TH2D>("h_bar_energy_2D_bgofiducial_lastcut", "Energy bar; Reconstructed Energy [GeV]; Energy Bar [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_2d_energy_bins.size() -1, &bars_2d_energy_bins[0]);
    h_bars_last_layer_10MeV_bgofiducial_lastcut = std::make_shared<TH1D>("h_bars_last_layer_10MeV_bgofiducial_lastcut", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 22, 0, 22);
    h_bars_last_layer_10MeV_2D_bgofiducial_lastcut = std::make_shared<TH2D>("h_bars_last_layer_10MeV_2D_bgofiducial_lastcut", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);

    h_bar_energy_lateral_showering_lastcut = std::make_shared<TH2D>("h_bar_energy_lateral_showering_lastcut", "Mean energy bar; Reconstructed Energy [GeV]; Mean Bar Energy [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_energy_bins.size() -1, &bars_energy_bins[0]);
    h_bar_energy_2D_lateral_showering_lastcut = std::make_shared<TH2D>("h_bar_energy_2D_lateral_showering_lastcut", "Energy bar; Reconstructed Energy [GeV]; Energy Bar [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_2d_energy_bins.size() -1, &bars_2d_energy_bins[0]);
    h_bars_last_layer_10MeV_lateral_showering_lastcut = std::make_shared<TH1D>("h_bars_last_layer_10MeV_lateral_showering_lastcut", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 22, 0, 22);
    h_bars_last_layer_10MeV_2D_lateral_showering_lastcut = std::make_shared<TH2D>("h_bars_last_layer_10MeV_2D_lateral_showering_lastcut", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);

    h_maxrms = std::make_shared<TH1D>("h_maxrms", "Max RMS Layers", 1000, 0, 3000);
    h_maxrms_2D = std::make_shared<TH2D>("h_maxrms_2D", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);
    h_maxrms_no_trigger = std::make_shared<TH1D>("h_maxrms_no_trigger", "Max RMS Layers - No Trigger", 1000, 0, 3000);
    h_maxrms_2D_no_trigger = std::make_shared<TH2D>("h_maxrms_2D_no_trigger", "Max RMS - [> 1% energy content] - No Trigger; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

    h_maxrms_after_bgofiducial = std::make_shared<TH1D>("h_maxrms_after_bgofiducial", "Max RMS Layers", 1000, 0, 3000);
    h_maxrms_2D_after_bgofiducial = std::make_shared<TH2D>("h_maxrms_2D_after_bgofiducial", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);
    h_maxrms_no_trigger_after_bgofiducial = std::make_shared<TH1D>("h_maxrms_no_trigger_after_bgofiducial", "Max RMS Layers - No Trigger", 1000, 0, 3000);
    h_maxrms_2D_no_trigger_after_bgofiducial = std::make_shared<TH2D>("h_maxrms_2D_no_trigger_after_bgofiducial", "Max RMS - [> 1% energy content] - No Trigger; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

    h_maxrms_bgofiducial_lastcut = std::make_shared<TH1D>("h_maxrms_bgofiducial_lastcut", "Max RMS Layers", 1000, 0, 3000);
    h_maxrms_2D_bgofiducial_lastcut = std::make_shared<TH2D>("h_maxrms_2D_bgofiducial_lastcut", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);
    h_maxrms_no_trigger_bgofiducial_lastcut = std::make_shared<TH1D>("h_maxrms_no_trigger_bgofiducial_lastcut", "Max RMS Layers - No Trigger", 1000, 0, 3000);
    h_maxrms_2D_no_trigger_bgofiducial_lastcut = std::make_shared<TH2D>("h_maxrms_2D_no_trigger_bgofiducial_lastcut", "Max RMS - [> 1% energy content] - No Trigger; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

    h_maxrms_lateral_showering_lastcut = std::make_shared<TH1D>("h_maxrms_lateral_showering_lastcut", "Max RMS Layers", 1000, 0, 3000);
    h_maxrms_2D_lateral_showering_lastcut = std::make_shared<TH2D>("h_maxrms_2D_lateral_showering_lastcut", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

    h_bgoshower_top_X = std::make_shared<TH1D>("h_bgoshower_top_X", "BGO Shower X TOP Projection; Top X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_X = std::make_shared<TH1D>("h_bgoshower_bottom_X", "BGO Shower X BOTTOM Projection; BOTTOM X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_top_Y = std::make_shared<TH1D>("h_bgoshower_top_Y", "BGO Shower X TOP Projection; Top Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_Y = std::make_shared<TH1D>("h_bgoshower_bottom_Y", "BGO Shower Y BOTTOM Projection; BOTTOM Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);

    h_bgoshower_top_X_after_bgofiducial = std::make_shared<TH1D>("h_bgoshower_top_X_after_bgofiducial", "BGO Shower X TOP Projection; Top X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_X_after_bgofiducial = std::make_shared<TH1D>("h_bgoshower_bottom_X_after_bgofiducial", "BGO Shower X BOTTOM Projection; BOTTOM X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_top_Y_after_bgofiducial = std::make_shared<TH1D>("h_bgoshower_top_Y_after_bgofiducial", "BGO Shower X TOP Projection; Top Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_Y_after_bgofiducial = std::make_shared<TH1D>("h_bgoshower_bottom_Y_after_bgofiducial", "BGO Shower Y BOTTOM Projection; BOTTOM Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    
    h_bgoshower_top_X_bgofiducial_lastcut = std::make_shared<TH1D>("h_bgoshower_top_X_bgofiducial_lastcut", "BGO Shower X TOP Projection; Top X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_X_bgofiducial_lastcut = std::make_shared<TH1D>("h_bgoshower_bottom_X_bgofiducial_lastcut", "BGO Shower X BOTTOM Projection; BOTTOM X Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_top_Y_bgofiducial_lastcut = std::make_shared<TH1D>("h_bgoshower_top_Y_bgofiducial_lastcut", "BGO Shower X TOP Projection; Top Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_Y_bgofiducial_lastcut = std::make_shared<TH1D>("h_bgoshower_bottom_Y_bgofiducial_lastcut", "BGO Shower Y BOTTOM Projection; BOTTOM Y Projection [mm]", 150, -BGO_SideXY, BGO_SideXY);

    h_energy_fraction_sh_axis_contained = std::make_shared<TH1D>("h_energy_fraction_sh_axis_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_sh_axis_not_contained = std::make_shared<TH1D>("h_energy_fraction_sh_axis_not_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);

    h_energy_fraction_sh_axis_contained_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_sh_axis_contained_after_bgofiducial", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_sh_axis_not_contained_after_bgofiducial = std::make_shared<TH1D>("h_energy_fraction_sh_axis_not_contained_after_bgofiducial", "Energy fraction; Energy fraction; counts", 100, 0, 1);

    h_energy_fraction_sh_axis_contained_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_sh_axis_contained_bgofiducial_lastcut", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_sh_axis_not_contained_bgofiducial_lastcut = std::make_shared<TH1D>("h_energy_fraction_sh_axis_not_contained_bgofiducial_lastcut", "Energy fraction; Energy fraction; counts", 100, 0, 1);

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

    h_PSD_STK_X_match_energy_int_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_X_match_energy_int_psd_fiducial", "PSD - STK match X view + PSD fiducial - All Energies; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_energy_int_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_Y_match_energy_int_psd_fiducial", "PSD - STK match Y view + PSD fiducial - All Energies; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_100_250_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_X_match_100_250_psd_fiducial", "PSD - STK match X view + PSD fiducial - 100-250 GeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_100_250_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_Y_match_100_250_psd_fiducial", "PSD - STK match Y view + PSD fiducial - 100-250 GeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_250_500_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_X_match_250_500_psd_fiducial", "PSD - STK match X view + PSD fiducial - 250-500 GeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_250_500_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_Y_match_250_500_psd_fiducial", "PSD - STK match Y view + PSD fiducial - 250-500 GeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_500_1000_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_X_match_500_1000_psd_fiducial", "PSD - STK match X view + PSD fiducial - 500GeV-1TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_500_1000_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_Y_match_500_1000_psd_fiducial", "PSD - STK match Y view + PSD fiducial - 500GeV-1TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_1000_5000_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_X_match_1000_5000_psd_fiducial", "PSD - STK match X view + PSD fiducial - 1TeV-5TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_1000_5000_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_Y_match_1000_5000_psd_fiducial", "PSD - STK match Y view + PSD fiducial - 1TeV-5TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_X_match_5000_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_X_match_5000_psd_fiducial", "PSD - STK match X view + PSD fiducial - >5TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    h_PSD_STK_Y_match_5000_psd_fiducial = std::make_shared<TH1D>("h_PSD_STK_Y_match_5000_psd_fiducial", "PSD - STK match Y view + PSD fiducial - >5TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    
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

    h_BGOrec_sumRms_flast_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_bgofiducial_lastcut", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_20_100_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_20_100_bgofiducial_lastcut", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_100_250_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_100_250_bgofiducial_lastcut", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_250_500_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_250_500_bgofiducial_lastcut", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_500_1000_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_500_1000_bgofiducial_lastcut", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_1000_3000_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_1000_3000_bgofiducial_lastcut", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_3000_5000_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_3000_5000_bgofiducial_lastcut", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_5000_bgofiducial_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_5000_bgofiducial_lastcut", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);

    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_20_100", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_100_250", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_250_500", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_500_1000", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_1000_3000", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_3000_5000", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_remove_lateral_and_showering_5000", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);

    h_BGOrec_sumRms_flast_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_lateral_showering_lastcut", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_20_100_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_20_100_lateral_showering_lastcut", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_100_250_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_100_250_lateral_showering_lastcut", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_250_500_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_250_500_lateral_showering_lastcut", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_500_1000_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_500_1000_lateral_showering_lastcut", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_1000_3000_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_1000_3000_lateral_showering_lastcut", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_3000_5000_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_3000_5000_lateral_showering_lastcut", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_5000_lateral_showering_lastcut = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_5000_lateral_showering_lastcut", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    
    h_STK_tracks = std::make_shared<TH1D>("h_STK_tracks", "STK tracks", 1000, 0, 1000);
    h_STK_tracks_vs_energy = std::make_shared<TH2D>("h_STK_tracks_vs_energy", "STK tracks vs Energy; Reconstructed Energy [GeV]; STK Tracks", energy_nbins -1, &energy_binning[0], (int)stk_ntracks_binning.size() -1, &stk_ntracks_binning[0]);
    
    h_STK_X_clusters = std::make_shared<TH1D>("h_STK_X_clusters", "STK X clusters", 10, 0, 10);
    h_STK_Y_clusters = std::make_shared<TH1D>("h_STK_Y_clusters", "STK Y clusters", 10, 0, 10);
    h_STK_X_holes = std::make_shared<TH1D>("h_STK_X_holes", "STK X holes", 10, 0, 10);
    h_STK_Y_holes = std::make_shared<TH1D>("h_STK_Y_holes", "STK Y holes", 10, 0, 10);

    h_STK_X_clusters_best_track = std::make_shared<TH1D>("h_STK_X_clusters_best_track", "STK X clusters", 10, 0, 10);
    h_STK_Y_clusters_best_track = std::make_shared<TH1D>("h_STK_Y_clusters_best_track", "STK Y clusters", 10, 0, 10);
    h_STK_X_holes_best_track = std::make_shared<TH1D>("h_STK_X_holes_best_track", "STK X holes", 10, 0, 10);
    h_STK_Y_holes_best_track = std::make_shared<TH1D>("h_STK_Y_holes_best_track", "STK Y holes", 10, 0, 10);

    h_STK_best_track_clusters = std::make_shared<TH1D>("h_STK_best_track_clusters", "STK - best track XY clusters", 10, 0, 10);

    h_STK_X_clusters_vs_energy = std::make_shared<TH2D>("h_STK_X_clusters_vs_energy", "STK X clusters vs Particle Energy; Reconstructed Energy [GeV]; STK X clusters", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);
    h_STK_Y_clusters_vs_energy = std::make_shared<TH2D>("h_STK_Y_clusters_vs_energy", "STK Y clusters vs Particle Energy; Reconstructed Energy [GeV]; STK Y clusters", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);
    h_STK_X_holes_vs_energy = std::make_shared<TH2D>("h_STK_X_holes_vs_energy", "STK X holes vs Particle Energy; Reconstructed Energy [GeV]; STK X holes", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);
    h_STK_Y_holes_vs_energy = std::make_shared<TH2D>("h_STK_Y_holes_vs_energy", "STK Y holes vs Particle Energy; Reconstructed Energy [GeV]; STK Y holes", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);

    h_STK_X_clusters_vs_energy_best_track = std::make_shared<TH2D>("h_STK_X_clusters_vs_energy_best_track", "STK X clusters vs Particle Energy; Reconstructed Energy [GeV]; STK X clusters", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);
    h_STK_Y_clusters_vs_energy_best_track = std::make_shared<TH2D>("h_STK_Y_clusters_vs_energy_best_track", "STK Y clusters vs Particle Energy; Reconstructed Energy [GeV]; STK Y clusters", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);
    h_STK_X_holes_vs_energy_best_track = std::make_shared<TH2D>("h_STK_X_holes_vs_energy_best_track", "STK X holes vs Particle Energy; Reconstructed Energy [GeV]; STK X holes", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);
    h_STK_Y_holes_vs_energy_best_track = std::make_shared<TH2D>("h_STK_Y_holes_vs_energy_best_track", "STK Y holes vs Particle Energy; Reconstructed Energy [GeV]; STK Y holes", energy_nbins -1, &energy_binning[0], (int)stk_track_binning.size() -1, &stk_track_binning[0]);

    h_STK_BGO_TOP_spatial_difference = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_3_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_3_clusters", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_3_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_3_clusters", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_3_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_3_clusters", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_3_clusters = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_3_clusters", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_4_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_4_clusters", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_4_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_4_clusters", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_4_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_4_clusters", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_4_clusters = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_4_clusters", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_5_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_5_clusters", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_5_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_5_clusters", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_5_clusters = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_5_clusters", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_5_clusters = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_5_clusters", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_best_track", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_best_track", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_best_track", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_best_track = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_best_track", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_3_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_3_clusters_best_track", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_3_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_3_clusters_best_track", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_3_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_3_clusters_best_track", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_3_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_3_clusters_best_track", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_4_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_4_clusters_best_track", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_4_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_4_clusters_best_track", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_4_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_4_clusters_best_track", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_4_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_4_clusters_best_track", "STK - BGO track angular difference", 100, 0, 40);

    h_STK_BGO_TOP_spatial_difference_5_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_difference_5_clusters_best_track", "STK - BGO TOP spatial difference", 80, 0, BGO_SideXY);
    h_STK_BGO_TOP_spatial_X_difference_5_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_X_difference_5_clusters_best_track", "STK - BGO TOP spatial difference - X view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_TOP_spatial_Y_difference_5_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_TOP_spatial_Y_difference_5_clusters_best_track", "STK - BGO TOP spatial difference - Y view", 100, -BGO_SideXY, BGO_SideXY);
    h_STK_BGO_track_angular_difference_5_clusters_best_track = std::make_shared<TH1D>("h_STK_BGO_track_angular_difference_5_clusters_best_track", "STK - BGO track angular difference", 100, 0, 40);

    h_BGOrec_sumRms_flast_after_track_selection = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection", "F_{last} vs sumRms correlation; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_20_100 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_20_100", "F_{last} vs sumRms correlation - 20 GeV - 100 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_100_250 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_100_250", "F_{last} vs sumRms correlation - 100 GeV - 250 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_250_500 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_250_500", "F_{last} vs sumRms correlation - 250 GeV - 500 GeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_500_1000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_500_1000", "F_{last} vs sumRms correlation - 500 GeV - 1 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_1000_3000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_1000_3000", "F_{last} vs sumRms correlation - 1 TeV - 3 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_3000_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_3000_5000", "F_{last} vs sumRms correlation - 3 TeV - 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);
    h_BGOrec_sumRms_flast_after_track_selection_5000 = std::make_shared<TH2D>("h_BGOrec_sumRms_flast_after_track_selection_5000", "F_{last} vs sumRms correlation - > 5 TeV; sumRMS [mm]; F_{last}", (int)sumRms_binning.size() - 1, &sumRms_binning[0], (int)flast_binning.size() - 1, &flast_binning[0]);

    h_STK_charge_X = std::make_shared<TH1D>("h_STK_charge_X", "STK Charge - X view; STK Charge X; entries", 400, 0, 100);
    h_STK_charge_Y = std::make_shared<TH1D>("h_STK_charge_Y", "STK Charge - Y view; STK Charge Y; entries", 400, 0, 100);
    h_STK_charge = std::make_shared<TH1D>("h_STK_charge", "STK Charge; STK Charge; entries", 400, 0, 100);
    h_STK_charge_2D = std::make_shared<TH2D>("h_STK_charge_2D", "STK Charge; STK Charge X; STK Charge Y", 400, 0, 100, 400, 0, 100);

    h_STK_charge_X_3_clusters = std::make_shared<TH1D>("h_STK_charge_X_3_clusters", "STK Charge - X view; STK Charge X; entries", 400, 0, 100);
    h_STK_charge_Y_3_clusters = std::make_shared<TH1D>("h_STK_charge_Y_3_clusters", "STK Charge - Y view; STK Charge Y; entries", 400, 0, 100);
    h_STK_charge_3_clusters = std::make_shared<TH1D>("h_STK_charge_3_clusters", "STK Charge; STK Charge; entries", 400, 0, 100);
    h_STK_charge_2D_3_clusters = std::make_shared<TH2D>("h_STK_charge_2D_3_clusters", "STK Charge; STK Charge X; STK Charge Y", 400, 0, 100, 400, 0, 100);

    h_STK_charge_X_4_clusters = std::make_shared<TH1D>("h_STK_charge_X_4_clusters", "STK Charge - X view; STK Charge X; entries", 400, 0, 100);
    h_STK_charge_Y_4_clusters = std::make_shared<TH1D>("h_STK_charge_Y_4_clusters", "STK Charge - Y view; STK Charge Y; entries", 400, 0, 100);
    h_STK_charge_4_clusters = std::make_shared<TH1D>("h_STK_charge_4_clusters", "STK Charge; STK Charge; entries", 400, 0, 100);
    h_STK_charge_2D_4_clusters = std::make_shared<TH2D>("h_STK_charge_2D_4_clusters", "STK Charge; STK Charge X; STK Charge Y", 400, 0, 100, 400, 0, 100);

    h_STK_charge_X_5_clusters = std::make_shared<TH1D>("h_STK_charge_X_5_clusters", "STK Charge - X view; STK Charge X; entries", 400, 0, 100);
    h_STK_charge_Y_5_clusters = std::make_shared<TH1D>("h_STK_charge_Y_5_clusters", "STK Charge - Y view; STK Charge Y; entries", 400, 0, 100);
    h_STK_charge_5_clusters = std::make_shared<TH1D>("h_STK_charge_5_clusters", "STK Charge; STK Charge; entries", 400, 0, 100);
    h_STK_charge_2D_5_clusters = std::make_shared<TH2D>("h_STK_charge_2D_5_clusters", "STK Charge; STK Charge X; STK Charge Y", 400, 0, 100, 400, 0, 100);

    h_PSD_charge_X = std::make_shared<TH1D>("h_PSD_charge_X", "PSD Charge - X view; PSD Charge X; entries", 300, 0, 40);
    h_PSD_charge_Y = std::make_shared<TH1D>("h_PSD_charge_Y", "PSD Charge - Y view; PSD Charge Y; entries", 300, 0, 40);
    h_PSD_charge = std::make_shared<TH1D>("h_PSD_charge", "PSD Charge; PSD Charge; entries", 300, 0, 40);
    h_PSD_charge_2D = std::make_shared<TH2D>("h_PSD_charge_2D", "PSD Charge; PSD Charge X; PSD Charge Y", 300, 0, 40, 300, 0, 40);
    h_PSD_sum_of_XY_charges = std::make_shared<TH1D>("h_PSD_sum_of_XY_charges", "PSD Charge (X+Y); PSD Charge (X+Y); entries", 800, 0, 200);

    h_PSD_X_clusters = std::make_shared<TH2D>("h_PSD_X_clusters", "PSD X Clusters vs Reconstructed Energy; Reconstructed Energy [GeV]; PSD X Clusters;", energy_nbins -1, &energy_binning[0], (int)psd_clusters_binning.size() -1, &psd_clusters_binning[0]);
    h_PSD_Y_clusters = std::make_shared<TH2D>("h_PSD_Y_clusters", "PSD Y Clusters vs Reconstructed Energy; Reconstructed Energy [GeV]; PSD Y Clusters;", energy_nbins -1, &energy_binning[0], (int)psd_clusters_binning.size() -1, &psd_clusters_binning[0]);

    h_STK_charge_X_nocut = std::make_shared<TH1D>("h_STK_charge_X_nocut", "STK Charge - X view; STK Charge X; entries", 400, 0, 100);
    h_STK_charge_Y_nocut = std::make_shared<TH1D>("h_STK_charge_Y_nocut", "STK Charge - Y view; STK Charge Y; entries", 400, 0, 100);
    h_STK_charge_nocut = std::make_shared<TH1D>("h_STK_charge_nocut", "STK Charge; STK Charge; entries", 400, 0, 100);
    h_STK_charge_2D_nocut = std::make_shared<TH2D>("h_STK_charge_2D_nocut", "STK Charge; STK Charge X; STK Charge Y", 400, 0, 100, 400, 0, 100);

    h_STK_charge_X_PSD_charge_cut = std::make_shared<TH1D>("h_STK_charge_X_PSD_charge_cut", "STK Charge - X view; STK Charge X; entries", 400, 0, 100);
    h_STK_charge_Y_PSD_charge_cut = std::make_shared<TH1D>("h_STK_charge_Y_PSD_charge_cut", "STK Charge - Y view; STK Charge Y; entries", 400, 0, 100);
    h_STK_charge_PSD_charge_cut = std::make_shared<TH1D>("h_STK_charge_PSD_charge_cut", "STK Charge; STK Charge; entries", 400, 0, 100);
    h_STK_charge_2D_PSD_charge_cut = std::make_shared<TH2D>("h_STK_charge_2D_PSD_charge_cut", "STK Charge; STK Charge X; STK Charge Y", 400, 0, 100, 400, 0, 100);

    h_PSD_charge_X_nocut = std::make_shared<TH1D>("h_PSD_charge_X_nocut", "PSD Charge - X view; PSD Charge X; entries", 300, 0, 40);
    h_PSD_charge_Y_nocut = std::make_shared<TH1D>("h_PSD_charge_Y_nocut", "PSD Charge - Y view; PSD Charge Y; entries", 300, 0, 40);
    h_PSD_charge_nocut = std::make_shared<TH1D>("h_PSD_charge_nocut", "PSD Charge; PSD Charge; entries", 300, 0, 40);
    h_PSD_charge_2D_nocut = std::make_shared<TH2D>("h_PSD_charge_2D_nocut", "PSD Charge; PSD Charge X; PSD Charge Y", 300, 0, 40, 300, 0, 40);
    h_PSD_sum_of_XY_charges_nocut = std::make_shared<TH1D>("h_PSD_sum_of_XY_charges_nocut", "PSD Charge (X+Y); PSD Charge (X+Y); entries", 1000, 0, 200);

    h_PSD_charge_X_STK_charge_cut = std::make_shared<TH1D>("h_PSD_charge_X_STK_charge_cut", "PSD Charge - X view; PSD Charge X; entries", 300, 0, 40);
    h_PSD_charge_Y_STK_charge_cut = std::make_shared<TH1D>("h_PSD_charge_Y_STK_charge_cut", "PSD Charge - Y view; PSD Charge Y; entries", 300, 0, 40);
    h_PSD_charge_STK_charge_cut = std::make_shared<TH1D>("h_PSD_charge_STK_charge_cut", "PSD Charge; PSD Charge; entries", 300, 0, 40);
    h_PSD_charge_2D_STK_charge_cut = std::make_shared<TH2D>("h_PSD_charge_2D_STK_charge_cut", "PSD Charge; PSD Charge X; PSD Charge Y", 300, 0, 40, 300, 0, 40);
    h_PSD_sum_of_XY_charges_STK_charge_cut = std::make_shared<TH1D>("h_PSD_sum_of_XY_charges_STK_charge_cut", "PSD Charge (X+Y); PSD Charge (X+Y); entries", 800, 0, 200);

    if (h_simu) {
        h_max_bar_position_simu_reco_energy_diff = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
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

        h_max_bar_position_simu_reco_energy_diff_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_bgoshower_in", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_in", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_in", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_in", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_in", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_in", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_in", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_in", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_in", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_in", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_in", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_in", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_in", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_in", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_in", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_bgoshower_out", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_out", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_out", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_out", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_out", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_out", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_out", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_out", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_out", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_out", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_out", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_out", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_out", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_out", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_out", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_after_bgofiducial", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_in", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_out", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_in = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_in", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 0; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 1; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 2; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 3; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 4; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 5; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 6; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 7; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 8; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 9; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 10; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 11; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 12; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_out = std::make_shared<TH2D>("h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_out", "Energy diff vs Max Bar Position - layer 13; Max Bar; Energy_{simu} - Energy_{reco} / Energy_{simu}", 22, 0, 22, 100, -0.2, 1);

        h_bgoshower_top_X_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_top_X_simu_reco_energy_diff", "Energy diff vs TOP X BGO position; TOP X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_X_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_bottom_X_simu_reco_energy_diff", "Energy diff vs BOTTOM X BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_top_Y_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_top_Y_simu_reco_energy_diff", "Energy diff vs TOP Y BGO position; TOP Y BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_Y_simu_reco_energy_diff = std::make_shared<TH2D>("h_bgoshower_bottom_Y_simu_reco_energy_diff", "Energy diff vs BOTTOM Y BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);

        h_bgoshower_top_X_simu_reco_energy_diff_after_bgofiducial = std::make_shared<TH2D>("h_bgoshower_top_X_simu_reco_energy_diff_after_bgofiducial", "Energy diff vs TOP X BGO position; TOP X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_X_simu_reco_energy_diff_after_bgofiducial = std::make_shared<TH2D>("h_bgoshower_bottom_X_simu_reco_energy_diff_after_bgofiducial", "Energy diff vs BOTTOM X BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_top_Y_simu_reco_energy_diff_after_bgofiducial = std::make_shared<TH2D>("h_bgoshower_top_Y_simu_reco_energy_diff_after_bgofiducial", "Energy diff vs TOP Y BGO position; TOP Y BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_Y_simu_reco_energy_diff_after_bgofiducial = std::make_shared<TH2D>("h_bgoshower_bottom_Y_simu_reco_energy_diff_after_bgofiducial", "Energy diff vs BOTTOM Y BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);

        h_bgoshower_top_X_simu_reco_energy_diff_bgofiducial_lastcut = std::make_shared<TH2D>("h_bgoshower_top_X_simu_reco_energy_diff_bgofiducial_lastcut", "Energy diff vs TOP X BGO position; TOP X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_X_simu_reco_energy_diff_bgofiducial_lastcut = std::make_shared<TH2D>("h_bgoshower_bottom_X_simu_reco_energy_diff_bgofiducial_lastcut", "Energy diff vs BOTTOM X BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_top_Y_simu_reco_energy_diff_bgofiducial_lastcut = std::make_shared<TH2D>("h_bgoshower_top_Y_simu_reco_energy_diff_bgofiducial_lastcut", "Energy diff vs TOP Y BGO position; TOP Y BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);
        h_bgoshower_bottom_Y_simu_reco_energy_diff_bgofiducial_lastcut = std::make_shared<TH2D>("h_bgoshower_bottom_Y_simu_reco_energy_diff_bgofiducial_lastcut", "Energy diff vs BOTTOM Y BGO position; BOTTOM X BGO position [mm]; Energy_{simu} - Energy_{reco} / Energy_{simu}", 150, -BGO_SideXY, BGO_SideXY, 100, -0.1, 1);

        h_diff_bgo_simu_particle_direction_X = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_X", "X view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_particle_direction_Y = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_Y", "Y view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_extr_top_position_X = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_X", "TOP X ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
        h_diff_bgo_simu_extr_top_position_Y = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_Y", "TOP Y ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);

        h_diff_bgo_simu_particle_direction_X_after_bgofiducial = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_X_after_bgofiducial", "X view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_particle_direction_Y_after_bgofiducial = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_Y_after_bgofiducial", "Y view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_extr_top_position_X_after_bgofiducial = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_X_after_bgofiducial", "TOP X ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
        h_diff_bgo_simu_extr_top_position_Y_after_bgofiducial = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_Y_after_bgofiducial", "TOP Y ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
    
        h_diff_bgo_simu_particle_direction_X_bgofiducial_lastcut = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_X_bgofiducial_lastcut", "X view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_particle_direction_Y_bgofiducial_lastcut = std::make_shared<TH1D>("h_diff_bgo_simu_particle_direction_Y_bgofiducial_lastcut", "Y view - slope_{BGO} - slope_{simu}; slope_{BGO} - slope_{simu} (deg); entries", 500, -40, 40);
        h_diff_bgo_simu_extr_top_position_X_bgofiducial_lastcut = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_X_bgofiducial_lastcut", "TOP X ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
        h_diff_bgo_simu_extr_top_position_Y_bgofiducial_lastcut = std::make_shared<TH1D>("h_diff_bgo_simu_extr_top_position_Y_bgofiducial_lastcut", "TOP Y ext - position_{BGO} - position_{simu}; position_{BGO} - position_{simu} [mm]; entries", 200, -200, 200);
    }
}

void histos::SetWeight(std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries, const double evt_corr_energy_gev) {
    if (simu_primaries != nullptr) {
        auto pdgid = simu_primaries->pvpart_pdg;
        switch (pdgid) {
            case 11:
                weight = pow(evt_corr_energy_gev, -2);
                break;
            case 2212:
                weight = pow(evt_corr_energy_gev, -1.7);
                break;
            default:
                std::cerr << "\n\nWrong particle ID [" << pdgid << "] --> ... Neither electron [11] or proton [2212] ... EXIT\n\n";
                exit(200);
                break;
        }
    }
}

const double histos::GetWeight() {
    return weight;
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
    h_energy_fraction_2D->Write();
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

        h_max_bar_position_simu_reco_energy_diff_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_in->Write();

        h_max_bar_position_simu_reco_energy_diff_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgoshower_out->Write();

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

    h_energy_fraction_after_bgofiducial->Write();
    h_energy_fraction_2D_after_bgofiducial->Write();
    h_energy_fraction_no_trigger_after_bgofiducial->Write();
    h_energy_fraction_large_angles_56_after_bgofiducial->Write();
    h_energy_fraction_large_angles_60_after_bgofiducial->Write();
    h_energy_fraction_large_angles_61_after_bgofiducial->Write();
    h_energy_fraction_large_angles_62_after_bgofiducial->Write();
    h_energy_fraction_large_angles_63_after_bgofiducial->Write();
    h_energy_fraction_large_angles_64_after_bgofiducial->Write();
    h_energy_fraction_large_angles_65_after_bgofiducial->Write();
    h_energy_fraction_large_angles_70_after_bgofiducial->Write();
    h_energy_fraction_large_angles_80_after_bgofiducial->Write();
    h_energy_fraction_large_angles_85_after_bgofiducial->Write();

    h_bar_energy_after_bgofiducial->Write();
    h_bar_energy_2D_after_bgofiducial->Write();
    h_bars_last_layer_10MeV_after_bgofiducial->Write();
    h_bars_last_layer_10MeV_2D_after_bgofiducial->Write();

    h_maxrms_after_bgofiducial->Write();
    h_maxrms_2D_after_bgofiducial->Write();
    h_maxrms_no_trigger_after_bgofiducial->Write();
    h_maxrms_2D_no_trigger_after_bgofiducial->Write();

    h_bgoshower_top_X_after_bgofiducial->Write();
    h_bgoshower_bottom_X_after_bgofiducial->Write();
    h_bgoshower_top_Y_after_bgofiducial->Write();
    h_bgoshower_bottom_Y_after_bgofiducial->Write();

    h_energy_fraction_sh_axis_contained_after_bgofiducial->Write();
    h_energy_fraction_sh_axis_not_contained_after_bgofiducial->Write();

    h_BGOrec_sumRms_flast_after_bgofiducial->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_20_100->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_100_250->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_250_500->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_500_1000->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_1000_3000->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_3000_5000->Write();
    h_BGOrec_sumRms_flast_after_bgofiducial_5000->Write();

    if (h_simu) {
        h_max_bar_position_simu_reco_energy_diff_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial->Write();

        h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_in->Write();

        h_max_bar_position_simu_reco_energy_diff_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_after_bgofiducial_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_after_bgofiducial_bgoshower_out->Write();

        h_bgoshower_top_X_simu_reco_energy_diff_after_bgofiducial->Write();
        h_bgoshower_bottom_X_simu_reco_energy_diff_after_bgofiducial->Write();
        h_bgoshower_top_Y_simu_reco_energy_diff_after_bgofiducial->Write();
        h_bgoshower_bottom_Y_simu_reco_energy_diff_after_bgofiducial->Write();

        h_diff_bgo_simu_particle_direction_X_after_bgofiducial->Write();
        h_diff_bgo_simu_particle_direction_Y_after_bgofiducial->Write();
        h_diff_bgo_simu_extr_top_position_X_after_bgofiducial->Write();
        h_diff_bgo_simu_extr_top_position_Y_after_bgofiducial->Write();
    }

    outfile->mkdir("BGO_fiducial_lastcut");
    outfile->cd("BGO_fiducial_lastcut");

    h_energy_fraction_bgofiducial_lastcut->Write();
    h_energy_fraction_2D_bgofiducial_lastcut->Write();
    h_energy_fraction_no_trigger_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_56_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_60_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_61_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_62_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_63_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_64_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_65_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_70_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_80_bgofiducial_lastcut->Write();
    h_energy_fraction_large_angles_85_bgofiducial_lastcut->Write();

    h_bar_energy_bgofiducial_lastcut->Write();
    h_bar_energy_2D_bgofiducial_lastcut->Write();
    h_bars_last_layer_10MeV_bgofiducial_lastcut->Write();
    h_bars_last_layer_10MeV_2D_bgofiducial_lastcut->Write();

    h_maxrms_bgofiducial_lastcut->Write();
    h_maxrms_2D_bgofiducial_lastcut->Write();
    h_maxrms_no_trigger_bgofiducial_lastcut->Write();
    h_maxrms_2D_no_trigger_bgofiducial_lastcut->Write();

    h_bgoshower_top_X_bgofiducial_lastcut->Write();
    h_bgoshower_bottom_X_bgofiducial_lastcut->Write();
    h_bgoshower_top_Y_bgofiducial_lastcut->Write();
    h_bgoshower_bottom_Y_bgofiducial_lastcut->Write();

    h_energy_fraction_sh_axis_contained_bgofiducial_lastcut->Write();
    h_energy_fraction_sh_axis_not_contained_bgofiducial_lastcut->Write();

    h_BGOrec_sumRms_flast_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_20_100_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_100_250_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_250_500_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_500_1000_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_1000_3000_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_3000_5000_bgofiducial_lastcut->Write();
    h_BGOrec_sumRms_flast_5000_bgofiducial_lastcut->Write();

    if (h_simu) {
        h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut->Write();

        h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_in->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_in->Write();

        h_max_bar_position_simu_reco_energy_diff_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_0_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_1_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_2_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_3_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_4_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_5_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_6_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_7_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_8_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_9_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_10_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_11_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_12_bgofiducial_lastcut_bgoshower_out->Write();
        h_max_bar_position_simu_reco_energy_diff_ly_13_bgofiducial_lastcut_bgoshower_out->Write();

        h_bgoshower_top_X_simu_reco_energy_diff_bgofiducial_lastcut->Write();
        h_bgoshower_bottom_X_simu_reco_energy_diff_bgofiducial_lastcut->Write();
        h_bgoshower_top_Y_simu_reco_energy_diff_bgofiducial_lastcut->Write();
        h_bgoshower_bottom_Y_simu_reco_energy_diff_bgofiducial_lastcut->Write();

        h_diff_bgo_simu_particle_direction_X_bgofiducial_lastcut->Write();
        h_diff_bgo_simu_particle_direction_Y_bgofiducial_lastcut->Write();
        h_diff_bgo_simu_extr_top_position_X_bgofiducial_lastcut->Write();
        h_diff_bgo_simu_extr_top_position_Y_bgofiducial_lastcut->Write();
    }

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

    outfile->mkdir("lateral_showering_lastcut");
    outfile->cd("lateral_showering_lastcut");

    h_bar_energy_lateral_showering_lastcut->Write();
    h_bar_energy_2D_lateral_showering_lastcut->Write();
    h_bars_last_layer_10MeV_lateral_showering_lastcut->Write();
    h_bars_last_layer_10MeV_2D_lateral_showering_lastcut->Write();

    h_maxrms_lateral_showering_lastcut->Write();
    h_maxrms_2D_lateral_showering_lastcut->Write();

    h_BGOrec_sumRms_flast_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_20_100_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_100_250_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_250_500_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_500_1000_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_1000_3000_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_3000_5000_lateral_showering_lastcut->Write();
    h_BGOrec_sumRms_flast_5000_lateral_showering_lastcut->Write();  

    outfile->mkdir("TrackSelection");
    outfile->cd("TrackSelection");

    h_STK_tracks->Write();
    h_STK_tracks_vs_energy->Write();

    h_STK_X_clusters->Write();
    h_STK_Y_clusters->Write();
    h_STK_X_holes->Write();
    h_STK_Y_holes->Write();

    h_STK_X_clusters_vs_energy->Write();
    h_STK_Y_clusters_vs_energy->Write();
    h_STK_X_holes_vs_energy->Write();
    h_STK_Y_holes_vs_energy->Write();

    h_STK_BGO_TOP_spatial_difference->Write();
    h_STK_BGO_TOP_spatial_X_difference->Write();
    h_STK_BGO_TOP_spatial_Y_difference->Write();
    h_STK_BGO_track_angular_difference->Write();

    h_STK_BGO_TOP_spatial_difference_3_clusters->Write();
    h_STK_BGO_TOP_spatial_X_difference_3_clusters->Write();
    h_STK_BGO_TOP_spatial_Y_difference_3_clusters->Write();
    h_STK_BGO_track_angular_difference_3_clusters->Write();

    h_STK_BGO_TOP_spatial_difference_4_clusters->Write();
    h_STK_BGO_TOP_spatial_X_difference_4_clusters->Write();
    h_STK_BGO_TOP_spatial_Y_difference_4_clusters->Write();
    h_STK_BGO_track_angular_difference_4_clusters->Write();

    h_STK_BGO_TOP_spatial_difference_5_clusters->Write();
    h_STK_BGO_TOP_spatial_X_difference_5_clusters->Write();
    h_STK_BGO_TOP_spatial_Y_difference_5_clusters->Write();
    h_STK_BGO_track_angular_difference_5_clusters->Write();

    h_STK_BGO_TOP_spatial_difference_best_track->Write();
    h_STK_BGO_TOP_spatial_X_difference_best_track->Write();
    h_STK_BGO_TOP_spatial_Y_difference_best_track->Write();
    h_STK_BGO_track_angular_difference_best_track->Write();

    h_BGOrec_sumRms_flast_after_track_selection->Write();
    h_BGOrec_sumRms_flast_after_track_selection_20_100->Write();
    h_BGOrec_sumRms_flast_after_track_selection_100_250->Write();
    h_BGOrec_sumRms_flast_after_track_selection_250_500->Write();
    h_BGOrec_sumRms_flast_after_track_selection_500_1000->Write();
    h_BGOrec_sumRms_flast_after_track_selection_1000_3000->Write();
    h_BGOrec_sumRms_flast_after_track_selection_3000_5000->Write();
    h_BGOrec_sumRms_flast_after_track_selection_5000->Write();

    outfile->mkdir("TrackSelection/BestTrack");
    outfile->cd("TrackSelection/BestTrack");

    h_STK_X_clusters_best_track->Write();
    h_STK_Y_clusters_best_track->Write();
    h_STK_X_holes_best_track->Write();
    h_STK_Y_holes_best_track->Write();

    h_STK_best_track_clusters->Write();

    h_STK_X_clusters_vs_energy_best_track->Write();
    h_STK_Y_clusters_vs_energy_best_track->Write();
    h_STK_X_holes_vs_energy_best_track->Write();
    h_STK_Y_holes_vs_energy_best_track->Write();

    h_STK_BGO_TOP_spatial_difference_3_clusters_best_track->Write();
    h_STK_BGO_TOP_spatial_X_difference_3_clusters_best_track->Write();
    h_STK_BGO_TOP_spatial_Y_difference_3_clusters_best_track->Write();
    h_STK_BGO_track_angular_difference_3_clusters_best_track->Write();

    h_STK_BGO_TOP_spatial_difference_4_clusters_best_track->Write();
    h_STK_BGO_TOP_spatial_X_difference_4_clusters_best_track->Write();
    h_STK_BGO_TOP_spatial_Y_difference_4_clusters_best_track->Write();
    h_STK_BGO_track_angular_difference_4_clusters_best_track->Write();

    h_STK_BGO_TOP_spatial_difference_5_clusters_best_track->Write();
    h_STK_BGO_TOP_spatial_X_difference_5_clusters_best_track->Write();
    h_STK_BGO_TOP_spatial_Y_difference_5_clusters_best_track->Write();
    h_STK_BGO_track_angular_difference_5_clusters_best_track->Write();

    outfile->mkdir("TrackSelection/Charge");
    outfile->cd("TrackSelection/Charge");

    h_STK_charge_X->Write();
    h_STK_charge_Y->Write();
    h_STK_charge->Write();
    h_STK_charge_2D->Write();

    h_STK_charge_X_3_clusters->Write();
    h_STK_charge_Y_3_clusters->Write();
    h_STK_charge_3_clusters->Write();
    h_STK_charge_2D_3_clusters->Write();

    h_STK_charge_X_4_clusters->Write();
    h_STK_charge_Y_4_clusters->Write();
    h_STK_charge_4_clusters->Write();
    h_STK_charge_2D_4_clusters->Write();

    h_STK_charge_X_5_clusters->Write();
    h_STK_charge_Y_5_clusters->Write();
    h_STK_charge_5_clusters->Write();
    h_STK_charge_2D_5_clusters->Write();

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

    h_PSD_STK_X_match_energy_int_psd_fiducial->Write();
    h_PSD_STK_Y_match_energy_int_psd_fiducial->Write();
    h_PSD_STK_X_match_100_250_psd_fiducial->Write();
    h_PSD_STK_Y_match_100_250_psd_fiducial->Write();
    h_PSD_STK_X_match_250_500_psd_fiducial->Write();
    h_PSD_STK_Y_match_250_500_psd_fiducial->Write();
    h_PSD_STK_X_match_500_1000_psd_fiducial->Write();
    h_PSD_STK_Y_match_500_1000_psd_fiducial->Write();
    h_PSD_STK_X_match_1000_5000_psd_fiducial->Write();
    h_PSD_STK_Y_match_1000_5000_psd_fiducial->Write();
    h_PSD_STK_X_match_5000_psd_fiducial->Write();
    h_PSD_STK_Y_match_5000_psd_fiducial->Write();

    h_PSD_X_clusters->Write();
    h_PSD_Y_clusters->Write();

    h_PSD_charge_X->Write();
    h_PSD_charge_Y->Write();
    h_PSD_charge->Write();
    h_PSD_charge_2D->Write();
    h_PSD_sum_of_XY_charges->Write();

    outfile->mkdir("PSD_charges");
    outfile->cd("PSD_charges");
    
    h_PSD_charge_X_nocut->Write();
    h_PSD_charge_Y_nocut->Write();
    h_PSD_charge_nocut->Write();
    h_PSD_charge_2D_nocut->Write();
    h_PSD_sum_of_XY_charges_nocut->Write();

    outfile->mkdir("PSD_charges_after_STK_charge_cut");
    outfile->cd("PSD_charges_after_STK_charge_cut");
    
    h_PSD_charge_X_STK_charge_cut->Write();
    h_PSD_charge_Y_STK_charge_cut->Write();
    h_PSD_charge_STK_charge_cut->Write();
    h_PSD_charge_2D_STK_charge_cut->Write();
    h_PSD_sum_of_XY_charges_STK_charge_cut->Write();

    outfile->mkdir("STK_charges");
    outfile->cd("STK_charges");

    h_STK_charge_X_nocut->Write();
    h_STK_charge_Y_nocut->Write();
    h_STK_charge_nocut->Write();
    h_STK_charge_2D_nocut->Write();

    outfile->mkdir("STK_charges_after_PSD_charge_cut");
    outfile->cd("STK_charges_after_PSD_charge_cut");

    h_STK_charge_X_PSD_charge_cut->Write();
    h_STK_charge_Y_PSD_charge_cut->Write();
    h_STK_charge_PSD_charge_cut->Write();
    h_STK_charge_2D_PSD_charge_cut->Write();

    outfile->Close();

    if (verbose) std::cout << "\n\nOutput file has been written [" << final_output_path << "]\n\n";  
}