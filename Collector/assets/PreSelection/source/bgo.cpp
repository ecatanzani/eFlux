#include "bgo.h"
#include "binning.h"
#include "energy_config.h"

#include "Dmp/DmpBgoContainer.h"
#include "Dmp/DmpStruct.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtSimuPrimaries.h"

#include "TH1D.h"
#include "TH2D.h"
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
    in_pars input_pars,
    std::shared_ptr<TChain> evtch) {

    if (input_pars.mc) bgofiducial_distributions_mc(input_pars, evtch);
    else bgofiducial_distributions_data(input_pars, evtch);    
}

void bgofiducial_distributions_mc(
    in_pars input_pars,
    std::shared_ptr<TChain> evtch) {

        // Load the energy config file
        std::unique_ptr<energy_config> econfig = std::make_unique<energy_config>(input_pars.config_wd);
        auto energy_binning = econfig->GetEnergyBinning();
        auto energy_nbins = (int)energy_binning.size() - 1;
        auto min_evt_energy = econfig->GetMinEvtEnergy();
        auto max_evt_energy = econfig->GetMaxEvtEnergy();
        
        // Register Header container
        std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
        evtch->SetBranchAddress("EventHeader", &evt_header);

        // Register BGO container
        std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
        evtch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

        // Register BGO REC container
        std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
        evtch->SetBranchAddress("DmpEvtBgoRec", &bgorec);

        // Register SimuPrimaries container
        std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries = std::make_shared<DmpEvtSimuPrimaries>();
        evtch->SetBranchAddress("DmpEvtSimuPrimaries", &simu_primaries);

        double layer_min_energy {0};           //Minimum energy per layer
        double gev {0.001};
        double evt_corr_energy_gev {0};
        double evt_energy_gev {0};
        double simu_energy_gev {0};
        auto nevents {evtch->GetEntries()};

        auto bars_energy_bins = createLogBinning(1e-2, 1e+5, 1000);
        auto number_of_bars_last_layer = createLinearBinning(0, 22, 23);
        auto max_rms_bins = createLinearBinning(0, 3000, 1000);

        std::unique_ptr<TH1D> h_energy_fraction = std::make_unique<TH1D>("h_energy_fraction", "Energy fraction; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_no_trigger = std::make_unique<TH1D>("h_energy_fraction_no_trigger", "Energy fraction - No triggered events; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_56 = std::make_unique<TH1D>("h_energy_fraction_large_angles_56", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_60 = std::make_unique<TH1D>("h_energy_fraction_large_angles_60", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_61 = std::make_unique<TH1D>("h_energy_fraction_large_angles_61", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_62 = std::make_unique<TH1D>("h_energy_fraction_large_angles_62", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_63 = std::make_unique<TH1D>("h_energy_fraction_large_angles_63", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_64 = std::make_unique<TH1D>("h_energy_fraction_large_angles_64", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_65 = std::make_unique<TH1D>("h_energy_fraction_large_angles_65", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_70 = std::make_unique<TH1D>("h_energy_fraction_large_angles_70", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_80 = std::make_unique<TH1D>("h_energy_fraction_large_angles_80", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_85 = std::make_unique<TH1D>("h_energy_fraction_large_angles_85", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);
        
        std::unique_ptr<TH2D> h_bar_energy = std::make_unique<TH2D>("h_bar_energy", "Mean energy bar; Reconstructed Energy [GeV]; Mean Bar Energy [GeV]", energy_nbins -1, &energy_binning[0], (int)bars_energy_bins.size() -1, &bars_energy_bins[0]);
        std::unique_ptr<TH1D> h_bars_last_layer_10MeV = std::make_unique<TH1D>("h_bars_last_layer_10MeV", "Number of bars on last layer - 10 MeV minimum; Number of BGO bars", 23, 0, 22);
        std::unique_ptr<TH2D> h_bars_last_layer_10MeV_2D = std::make_unique<TH2D>("h_bars_last_layer_10MeV_2D", "Number of bars on last layer - 10 MeV minimum; Reconstructed Energy [GeV]; Number of BGO bars", energy_nbins -1, &energy_binning[0], (int)number_of_bars_last_layer.size() -1, &number_of_bars_last_layer[0]);
        std::unique_ptr<TH1D> h_maxrms = std::make_unique<TH1D>("h_maxrms", "Max RMS Layers", 1000, 0, 3000);
        std::unique_ptr<TH2D> h_maxrms_2D = std::make_unique<TH2D>("h_maxrms_2D", "Max RMS - [> 1% energy content]; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);
        std::unique_ptr<TH1D> h_maxrms_no_trigger = std::make_unique<TH1D>("h_maxrms_no_trigger", "Max RMS Layers - No Trigger", 1000, 0, 3000);
        std::unique_ptr<TH2D> h_maxrms_2D_no_trigger = std::make_unique<TH2D>("h_maxrms_2D_no_trigger", "Max RMS - [> 1% energy content] - No Trigger; Reconstructed Energy [GeV]; RMS_{max}", energy_nbins -1, &energy_binning[0], (int)max_rms_bins.size()-1, &max_rms_bins[0]);

        std::unique_ptr<TH1D> h_bgoshower_top_X = std::make_unique<TH1D>("h_bgoshower_top_X", "BGO Shower X TOP Projection; Top X Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
        std::unique_ptr<TH1D> h_bgoshower_bottom_X = std::make_unique<TH1D>("h_bgoshower_bottom_X", "BGO Shower X BOTTOM Projection; BOTTOM X Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
        std::unique_ptr<TH1D> h_bgoshower_top_Y = std::make_unique<TH1D>("h_bgoshower_top_Y", "BGO Shower X TOP Projection; Top Y Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);
        std::unique_ptr<TH1D> h_bgoshower_bottom_Y = std::make_unique<TH1D>("h_bgoshower_bottom_Y", "BGO Shower Y BOTTOM Projection; BOTTOM Y Projection [mm]", 1000, -BGO_SideXY, BGO_SideXY);

        std::unique_ptr<TH1D> h_energy_fraction_sh_axis_contained = std::make_unique<TH1D>("h_energy_fraction_sh_axis_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);
        std::unique_ptr<TH1D> h_energy_fraction_sh_axis_not_contained = std::make_unique<TH1D>("h_energy_fraction_sh_axis_not_contained", "Energy fraction; Energy fraction; counts", 100, 0, 1);


        unsigned int step {10000};
        if (input_pars.verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;

        for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx) {

            evtch->GetEvent(evIdx);
            if (input_pars.verbose && !((evIdx+1)%step))
                std::cout << "\nNumber of processed events [" << evIdx+1 << "]";
            
            evt_energy_gev = bgorec->GetTotalEnergy()*gev;
            evt_corr_energy_gev = bgorec->GetElectronEcor()*gev;
            simu_energy_gev = simu_primaries->pvpart_ekin*gev;

            if (evt_corr_energy_gev>=min_evt_energy && evt_corr_energy_gev<=max_evt_energy) {

                std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
                bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), layer_min_energy);

                if (check_trigger(evt_header)) {
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        h_energy_fraction->Fill(elm_energy_fraction);
                    
                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(56*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_56->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(60*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_60->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(61*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_61->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(62*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_62->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(63*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_63->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(64*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_64->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(65*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_65->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(70*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_70->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(80*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_80->Fill(elm_energy_fraction);

                    if (bgorec->GetTrajectoryDirection2D().CosTheta()<cos(85*TMath::DegToRad()))
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_large_angles_85->Fill(elm_energy_fraction);

                    h_bar_energy->Fill(evt_corr_energy_gev, get_mean_bar_energy(bgoVault->GetLayerBarEnergies()));
                    h_bars_last_layer_10MeV->Fill(count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10));
                    h_bars_last_layer_10MeV_2D->Fill(evt_corr_energy_gev, count_bars_on_layer((bgoVault->GetLayerBarEnergies())[DAMPE_bgo_nLayers-1], 10));
                    h_maxrms->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
                    h_maxrms_2D->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));

                    std::vector<double> bgorec_proj_X, bgorec_proj_Y;
                    std::tie(bgorec_proj_X, bgorec_proj_Y) = get_shorew_axis_reco_projections(bgoVault->GetBGOslope(), bgoVault->GetBGOintercept());
                    h_bgoshower_top_X->Fill(bgorec_proj_X[0]);
                    h_bgoshower_bottom_X->Fill(bgorec_proj_X[1]);
                    h_bgoshower_top_Y->Fill(bgorec_proj_Y[0]);
                    h_bgoshower_bottom_Y->Fill(bgorec_proj_Y[1]);


                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                        if ((abs(bgorec_proj_X[0])<=BGO_SideXY && abs(bgorec_proj_X[1]<=BGO_SideXY)) && (abs(bgorec_proj_Y[0])<=BGO_SideXY && abs(bgorec_proj_Y[1])<=BGO_SideXY))
                            h_energy_fraction_sh_axis_contained->Fill(elm_energy_fraction);
                        else
                            h_energy_fraction_sh_axis_not_contained->Fill(elm_energy_fraction);
                }
                else {
                    if (evt_corr_energy_gev) {
                        for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                            h_energy_fraction_no_trigger->Fill(elm_energy_fraction);

                        h_maxrms_no_trigger->Fill(get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
                        h_maxrms_2D_no_trigger->Fill(evt_corr_energy_gev, get_max_rms(bgoVault->GetRmsLayer(), bgoVault->GetLayerEnergies(), bgorec->GetTotalEnergy()));
                    }
                }
            }
        }

        auto final_output_path = input_pars.output_wd+std::string("/bgofiducial.root");
        TFile* outfile = TFile::Open(final_output_path.c_str(), "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << final_output_path << "]\n\n";
            exit(100);
        }

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

        outfile->Close();

        if (input_pars.verbose) std::cout << "\n\nOutput file has been written [" << final_output_path << "]\n\n";   
    }

    void bgofiducial_distributions_data(
        in_pars input_pars,
        std::shared_ptr<TChain> evtch) {

        }