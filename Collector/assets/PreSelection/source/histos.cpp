#include "histos.h"
#include "binning.h"

#include <iostream>

#include "Dmp/DmpStruct.h"

#include "TFile.h"

histos::histos(std::shared_ptr<energy_config> econfig) {

    energy_binning = econfig->GetEnergyBinning();
    energy_nbins = (int)energy_binning.size() - 1;

    bars_energy_bins = createLogBinning(1e-2, 1e+5, 1000);
    number_of_bars_last_layer = createLinearBinning(0, 22, 23);
    max_rms_bins = createLinearBinning(0, 3000, 1000);

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
    h_bars_last_layer_10MeV = std::make_shared<TH1D>("h_bars_last_layer_10MeV", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 23, 0, 22);
    h_bars_last_layer_10MeV_2D = std::make_shared<TH2D>("h_bars_last_layer_10MeV_2D", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);
    h_maxrms = std::make_shared<TH1D>("h_maxrms", "Max RMS Layers", 1000, 0, 3000);
    h_maxrms_2D = std::make_shared<TH2D>("h_maxrms_2D", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);
    h_maxrms_no_trigger = std::make_shared<TH1D>("h_maxrms_no_trigger", "Max RMS Layers - No Trigger", 1000, 0, 3000);
    h_maxrms_2D_no_trigger = std::make_shared<TH2D>("h_maxrms_2D_no_trigger", "Max RMS - [> 1% energy content] - No Trigger; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

    h_bgoshower_top_X = std::make_shared<TH1D>("h_bgoshower_top_X", "BGO Shower X TOP Projection; Top X Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_X = std::make_shared<TH1D>("h_bgoshower_bottom_X", "BGO Shower X BOTTOM Projection; BOTTOM X Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_top_Y = std::make_shared<TH1D>("h_bgoshower_top_Y", "BGO Shower X TOP Projection; Top Y Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
    h_bgoshower_bottom_Y = std::make_shared<TH1D>("h_bgoshower_bottom_Y", "BGO Shower Y BOTTOM Projection; BOTTOM Y Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
    
    h_energy_fraction_sh_axis_contained = std::make_shared<TH1D>("h_energy_fraction_sh_axis_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);
    h_energy_fraction_sh_axis_not_contained = std::make_shared<TH1D>("h_energy_fraction_sh_axis_not_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);

    PSD_STK_X_match_energy_int = std::make_shared<TH1D>("PSD_STK_X_match_energy_int", "PSD - STK match X view- All Energies; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_Y_match_energy_int = std::make_shared<TH1D>("PSD_STK_Y_match_energy_int", "PSD - STK match Y view - All Energies; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_X_match_100_250 = std::make_shared<TH1D>("PSD_STK_X_match_100_250", "PSD - STK match X view - 100-250 GeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_Y_match_100_250 = std::make_shared<TH1D>("PSD_STK_Y_match_100_250", "PSD - STK match Y view - 100-250 GeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_X_match_250_500 = std::make_shared<TH1D>("PSD_STK_X_match_250_500", "PSD - STK match X view - 250-500 GeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_Y_match_250_500 = std::make_shared<TH1D>("PSD_STK_Y_match_250_500", "PSD - STK match Y view - 250-500 GeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_X_match_500_1000 = std::make_shared<TH1D>("PSD_STK_X_match_500_1000", "PSD - STK match X view - 500GeV-1TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_Y_match_500_1000 = std::make_shared<TH1D>("PSD_STK_Y_match_500_1000", "PSD - STK match Y view - 500GeV-1TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_X_match_1000_5000 = std::make_shared<TH1D>("PSD_STK_X_match_1000_5000", "PSD - STK match X view - 1TeV-5TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_Y_match_1000_5000 = std::make_shared<TH1D>("PSD_STK_Y_match_1000_5000", "PSD - STK match Y view - 1TeV-5TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_X_match_5000 = std::make_shared<TH1D>("PSD_STK_X_match_5000", "PSD - STK match X view - >5TeV; #Delta X (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
    PSD_STK_Y_match_5000 = std::make_shared<TH1D>("PSD_STK_Y_match_5000", "PSD - STK match Y view - >5TeV; #Delta Y (Track-closest PSD hit) [mm]; Entries", 200, -BGO_SideXY, BGO_SideXY);
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

    outfile->mkdir("PSD_STK");
    outfile->cd("PSD_STK");

    PSD_STK_X_match_energy_int->Write();
    PSD_STK_Y_match_energy_int->Write();
    PSD_STK_X_match_100_250->Write();
    PSD_STK_Y_match_100_250->Write();
    PSD_STK_X_match_250_500->Write();
    PSD_STK_Y_match_250_500->Write();
    PSD_STK_X_match_500_1000->Write();
    PSD_STK_Y_match_500_1000->Write();
    PSD_STK_X_match_1000_5000->Write();
    PSD_STK_Y_match_1000_5000->Write();
    PSD_STK_X_match_5000->Write();
    PSD_STK_Y_match_5000->Write();

    outfile->Close();

    if (verbose) std::cout << "\n\nOutput file has been written [" << final_output_path << "]\n\n";  
}