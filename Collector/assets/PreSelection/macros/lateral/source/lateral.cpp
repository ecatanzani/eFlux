#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>

#include "TKey.h"
#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>

#include "energy_config.h"

void buildLateralDistributions(
    const char* input_list,
    const char* output_file,
    const bool verbose,
    const bool mc_flag,
    const unsigned int threads,
    const std::string energy_config_file,
    const double input_spectral_index) {

        auto parse_input_list = [](const char* input_list, const bool verbose) -> std::shared_ptr<TChain> {

            auto get_tree_name = [](const char* file) -> std::string {
                std::unique_ptr<TFile> input_file = std::make_unique<TFile>(file, "READ");
                if (!input_file->IsOpen()) {
                    std::cerr << "\n\nError reading input file [" << file << "]\n\n";
                    exit(100);
                }
                std::string tree_name{""};
                for (TObject* keyAsObject : *input_file->GetListOfKeys()) {
                    auto key = dynamic_cast<TKey*>(keyAsObject);
                    if (!strcmp(key->GetClassName(), "TTree"))
                        tree_name = static_cast<std::string>(key->GetName());
                }
                return tree_name.c_str();
            };

            auto parser = [get_tree_name](const char* input_list, const bool verbose) -> std::shared_ptr<TChain> {
                std::shared_ptr<TChain> evtch {nullptr};
                std::ifstream input_file(input_list);
                if (!input_file.is_open()) {
                    std::cerr << "\n\nError (100) reading input file list...[" << input_list << "]" << std::endl;
                    exit(100);
                }
                std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
                input_file.close();

                std::istringstream input_stream(input_string);

                std::string tmp_str;
                bool first_elm {true};
                while (input_stream >> tmp_str) {
                    if (first_elm) {
                        evtch = std::make_shared<TChain>(get_tree_name(tmp_str.c_str()).c_str(), "Collector Tree");
                        first_elm = false;
                    }
                    evtch->Add(tmp_str.c_str());
                    if (verbose)
                        std::cout << "\nAdding " << tmp_str << " to the chain ...";
                }
                return evtch;
            };

            return parser(input_list, verbose);
        };
        
        auto evtch = parse_input_list(input_list, verbose);

        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame fr(*evtch);

        auto createLogBinning = [](const double eMin, const double eMax, const std::size_t n_bins) {
            std::vector<double> binning(n_bins + 1, 0);
            double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
            for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
                binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

            return binning;
        };

        auto sumRms_binning = createLogBinning(10, 2e+3, 1e+2);
        auto flast_binning = createLogBinning(1e-5, 2e-1, 1e+3);

        // Set event weight

        std::shared_ptr<energy_config> econfig = std::make_shared<energy_config>(energy_config_file);
        auto energy_binning = econfig->GetEnergyBinning();
       
        auto get_weight = [&input_spectral_index, &energy_binning, &mc_flag] (const double energy_gev) -> double {
            return !mc_flag ? 1 : std::pow(energy_gev, input_spectral_index+1)*std::pow(energy_binning[0], fabs(input_spectral_index+1));
        };

        // Flast vs sumRMS - only trigger

        auto sumrms_flast_lateral_distributions_only_trigger = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("evtfilter_evt_triggered==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_only_trigger", "sumrms_flast_lateral_distributions_only_trigger; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // Flast vs sumRMS - after lateral removing cuts

        auto sumrms_flast_lateral_distributions_20_100 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_20_100", "sumrms_flast_lateral_distributions_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_100_250 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_100_250", "sumrms_flast_lateral_distributions_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_250_500 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_250_500", "sumrms_flast_lateral_distributions_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_lateral_distributions_500_1000 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_500_1000", "sumrms_flast_lateral_distributions_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_1000_3000 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_1000_3000", "sumrms_flast_lateral_distributions_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_3000 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_3000", "sumrms_flast_lateral_distributions_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_lateral_distributions = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions", "sumrms_flast_lateral_distributions; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // Flast vs sumRMS - after lateral removing cuts and sumRMS low energy cut

        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100", "sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250", "sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500", "sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000", "sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000", "sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_3000 = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_3000", "sumrms_flast_lateral_distributions_after_lowenergysumrms_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_lateral_distributions_after_lowenergysumrms = fr.Define("energy_corr_gev", "energy_corr*0.001")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms", "sumrms_flast_lateral_distributions_after_lowenergysumrms; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        std::unique_ptr<TFile> output = std::make_unique<TFile>(output_file, "RECREATE");
        if (output->IsZombie()) {
            std::cerr << "Error: could not create output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        output->mkdir("sumrms_flast_only_trigger");
        output->cd("sumrms_flast_only_trigger");

        sumrms_flast_lateral_distributions_only_trigger->Write();

        output->mkdir("sumrms_flast_after_remove_lateral");
        output->cd("sumrms_flast_after_remove_lateral");

        sumrms_flast_lateral_distributions_20_100->Write();
        sumrms_flast_lateral_distributions_100_250->Write();
        sumrms_flast_lateral_distributions_250_500->Write();
        sumrms_flast_lateral_distributions_500_1000->Write();
        sumrms_flast_lateral_distributions_1000_3000->Write();
        sumrms_flast_lateral_distributions_3000->Write();    
        sumrms_flast_lateral_distributions->Write();

        output->mkdir("sumrms_flast_after_remove_lateral_lowenergysumrms");
        output->cd("sumrms_flast_after_remove_lateral_lowenergysumrms");

        sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_3000->Write();    
        sumrms_flast_lateral_distributions_after_lowenergysumrms->Write();                                                                                                 

        output->Close();
    }

