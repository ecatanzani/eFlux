#include "bgofiducial.h"

#include "Dmp/DmpBgoContainer.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtBgoRec.h"
#include "DmpEvtSimuPrimaries.h"

#include "TH1D.h"
#include "TMath.h"

#include <vector>
#include <fstream>

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
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	return trigger_HET;
}

void bgofiducial_distributions(
    const std::string output_path, 
    const std::string logs_dir,
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> evt_config,
    const bool mc, 
    const bool verbose) {

    if (mc) bgofiducial_distributions_mc(output_path, logs_dir, evtch, evt_config, verbose);
    else bgofiducial_distributions_data(output_path, logs_dir, evtch, evt_config, verbose);    
}

void bgofiducial_distributions_mc(
    const std::string output_path, 
    const std::string logs_dir,
    std::shared_ptr<TChain> evtch,
    std::shared_ptr<config> evt_config,
    const bool verbose)
    {
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
        double simu_energy_gev {0};
        auto nevents {evtch->GetEntries()};

        std::vector<unsigned int> ev_number_eratio_35;

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
        std::unique_ptr<TH1D> h_energy_fraction_large_angles_85 = std::make_unique<TH1D>("h_energy_fraction_large_angles_80", "Energy fraction Triggered - Large angles ; Energy fraction; counts", 100, 0, 1);

        unsigned int step {10000};
        if (verbose) std::cout << "\n\nNumber of events: " << nevents << std::endl;
        for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx) {
            evtch->GetEvent(evIdx);
            bool write_evt {true};
            if (verbose && !((evIdx+1)%step))
                std::cout << "\nNumber of processed events [" << evIdx+1 << "]";
            
            evt_corr_energy_gev = bgorec->GetElectronEcor()*gev;
            simu_energy_gev = simu_primaries->pvpart_ekin*gev;

            if (simu_energy_gev>=evt_config->GetMinEnergyRange() && simu_energy_gev<=evt_config->GetMaxEnergyRange()) {

                std::unique_ptr<DmpBgoContainer> bgoVault = std::make_unique<DmpBgoContainer>();
                bgoVault->scanBGOHits(bgohits, bgorec, bgorec->GetTotalEnergy(), layer_min_energy);

                if (check_trigger(evt_header)) {
                    for(auto&& elm_energy_fraction : bgoVault->GetFracLayer()) {
                        h_energy_fraction->Fill(elm_energy_fraction);
                        if (elm_energy_fraction>0.35 && write_evt && !logs_dir.empty()) {
                            ev_number_eratio_35.push_back(evIdx);
                            write_evt = false;
                        }
                    }

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

                }
                else {
                        if (evt_corr_energy_gev) {
                            for(auto&& elm_energy_fraction : bgoVault->GetFracLayer())
                                h_energy_fraction_no_trigger->Fill(elm_energy_fraction);
                        }
                }


            }
        }

        TFile* outfile = TFile::Open(output_path.c_str(), "RECREATE");
        if (!outfile->IsOpen()) {
            std::cerr << "\n\nError writing output file [" << output_path << "]\n\n";
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

        outfile->Close();

        std::ofstream log_energy_fraction_35;
        
        if (!logs_dir.empty()) {
            log_energy_fraction_35.open(logs_dir+std::string("/energy_fraction_35.log"));
            for (auto&& elm : ev_number_eratio_35)
                log_energy_fraction_35 << elm << std::endl;
            log_energy_fraction_35.close();
        }

        if (verbose) std::cout << "\n\nOutput file has been written [" << output_path << "]\n\n";   
    }

    void bgofiducial_distributions_data(
        const std::string output_path, 
        const std::string logs_dir,
        std::shared_ptr<TChain> evtch,
        std::shared_ptr<config> evt_config,
        const bool verbose) {

        }