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

#include "preselection.h"
#include "energy_config.h"

void preselection(
    const char* input_list,
    const char* output_file,
    const bool verbose,
    const bool mc_flag,
    const unsigned int threads,
    const char* energy_config_file,
    const double input_spectral_index) {

        // Parse energy config file
        std::shared_ptr<energy_config> econfig = std::make_shared<energy_config>(energy_config_file);
        auto energy_binning = econfig->GetEnergyBinning();

        // Parse input files
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

        // Create RDF
        ROOT::EnableImplicitMT(threads);
        ROOT::RDataFrame fr(*evtch);

        // Set event weight
        auto get_weight = [&input_spectral_index, &energy_binning, &mc_flag] (const double energy_gev) -> double {
            return !mc_flag ? 1 : std::pow(energy_gev, input_spectral_index+1)*std::pow(energy_binning[0], fabs(input_spectral_index+1));
        };

        // Energy filter
        auto energy_filter = [&energy_binning] (const double energy_gev) -> bool {
            return energy_gev >= energy_binning.front() && energy_gev <= energy_binning.back();
        };

        auto fr_energy_filtered = fr.Define("energy_corr_gev", "energy_corr*0.001").Filter(energy_filter, {"energy_corr_gev"});

        // Define binnings
        auto createLogBinning = [](const double eMin, const double eMax, const std::size_t n_bins) -> std::vector<double> {
            std::vector<double> binning(n_bins + 1, 0);
            double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
            for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
                binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

            return binning;
        };

        auto createLinearBinning = [] (const double eMin, const double eMax, const std::size_t n_bins) -> std::vector<double> {
            double h = (eMax - eMin) / static_cast<double>(n_bins);
            std::vector<double> xs(n_bins + 1);
            std::vector<double>::iterator x;
            double val;
            for (x = xs.begin(), val = eMin; x != xs.end(); ++x, val += h)
                *x = val;

            return xs;
        };

        auto sumRms_binning = createLogBinning(10, 2e+3, 2e+2);
        auto flast_binning = createLogBinning(1e-5, 2e-1, 2e+3);
        auto cosine_binning = createLinearBinning(0.3, 1, 100);

        // Flast vs sumRMS - only trigger
        auto sumrms_flast_lateral_distributions_only_trigger = fr_energy_filtered.Filter("evtfilter_evt_triggered==1")
                                                    .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                    .Histo2D({"sumrms_flast_lateral_distributions_only_trigger", "sumrms_flast_lateral_distributions_only_trigger; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // Flast vs sumRMS - after bgo fiducial cuts
        auto sumrms_flast_bgofiducial_distributions_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions_20_100", "sumrms_flast_bgofiducial_distributions_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_bgofiducial_distributions_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions_100_250", "sumrms_flast_bgofiducial_distributions_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_bgofiducial_distributions_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions_250_500", "sumrms_flast_bgofiducial_distributions_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_bgofiducial_distributions_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions_500_1000", "sumrms_flast_bgofiducial_distributions_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_bgofiducial_distributions_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions_1000_3000", "sumrms_flast_bgofiducial_distributions_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_bgofiducial_distributions_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions_3000", "sumrms_flast_bgofiducial_distributions_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_bgofiducial_distributions = fr_energy_filtered.Filter("evtfilter_BGO_fiducial_HET==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_bgofiducial_distributions", "sumrms_flast_bgofiducial_distributions; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // Flast vs sumRMS - after lateral removing cuts
        auto sumrms_flast_lateral_distributions_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_20_100", "sumrms_flast_lateral_distributions_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_100_250", "sumrms_flast_lateral_distributions_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_250_500", "sumrms_flast_lateral_distributions_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_lateral_distributions_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_500_1000", "sumrms_flast_lateral_distributions_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_1000_3000", "sumrms_flast_lateral_distributions_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_3000", "sumrms_flast_lateral_distributions_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_lateral_distributions = fr_energy_filtered.Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions", "sumrms_flast_lateral_distributions; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

         auto sumrms_cosine_lateral_distributions_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions_20_100", "sumrms_cosine_lateral_distributions_20_100; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");
        auto sumrms_cosine_lateral_distributions_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions_100_250", "sumrms_cosine_lateral_distributions_100_250; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");
        auto sumrms_cosine_lateral_distributions_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions_250_500", "sumrms_cosine_lateral_distributions_250_500; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");     
        auto sumrms_cosine_lateral_distributions_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions_500_1000", "sumrms_cosine_lateral_distributions_500_1000; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");
        auto sumrms_cosine_lateral_distributions_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions_1000_3000", "sumrms_cosine_lateral_distributions_1000_3000; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");
        auto sumrms_cosine_lateral_distributions_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions_3000", "sumrms_cosine_lateral_distributions_3000; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");

        auto sumrms_cosine_lateral_distributions = fr_energy_filtered.Filter("evtfilter_maxRms_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Define("bgorec_cosine", "BGOrec_trajectoryDirection2D.CosTheta()")
                                                .Histo2D({"sumrms_cosine_lateral_distributions", "sumrms_cosine_lateral_distributions; cos(#theta); sumRMS [mm]", static_cast<int>(cosine_binning.size()-1), cosine_binning.data(), static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data()}, "bgorec_cosine", "sumRms", "evt_w");

        // Flast vs sumRMS - after lateral removing cuts and sumRMS low energy cut
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100", "sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250", "sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500", "sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000", "sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000", "sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_lateral_distributions_after_lowenergysumrms_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms_3000", "sumrms_flast_lateral_distributions_after_lowenergysumrms_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_lateral_distributions_after_lowenergysumrms = fr_energy_filtered.Filter("evtfilter_sumrms_low_energy_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_lateral_distributions_after_lowenergysumrms", "sumrms_flast_lateral_distributions_after_lowenergysumrms; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // Flast vs sumRMS - after track selection cut
        auto sumrms_flast_trackselection_distributions_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions_20_100", "sumrms_flast_trackselection_distributions_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_trackselection_distributions_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions_100_250", "sumrms_flast_trackselection_distributions_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_trackselection_distributions_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions_250_500", "sumrms_flast_trackselection_distributions_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_trackselection_distributions_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions_500_1000", "sumrms_flast_trackselection_distributions_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_trackselection_distributions_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions_1000_3000", "sumrms_flast_trackselection_distributions_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_trackselection_distributions_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions_3000", "sumrms_flast_trackselection_distributions_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_trackselection_distributions = fr_energy_filtered.Filter("evtfilter_track_selection_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_trackselection_distributions", "sumrms_flast_trackselection_distributions; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // Flast vs sumRMS - after stk 1 RM cleaning cut
        auto sumrms_flast_stk1rm_distributions_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions_20_100", "sumrms_flast_stk1rm_distributions_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_stk1rm_distributions_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions_100_250", "sumrms_flast_stk1rm_distributions_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_stk1rm_distributions_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions_250_500", "sumrms_flast_stk1rm_distributions_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_stk1rm_distributions_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions_500_1000", "sumrms_flast_stk1rm_distributions_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_stk1rm_distributions_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions_1000_3000", "sumrms_flast_stk1rm_distributions_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_stk1rm_distributions_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions_3000", "sumrms_flast_stk1rm_distributions_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_stk1rm_distributions = fr_energy_filtered.Filter("evtfilter_stk_1rm_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_stk1rm_distributions", "sumrms_flast_stk1rm_distributions; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");


        // Flast vs sumRMS - after all cuts
        auto sumrms_flast_all_cuts_distributions_20_100 = fr_energy_filtered.Filter("energy_corr_gev >= 20 && energy_corr_gev < 100")
                                                .Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions_20_100", "sumrms_flast_all_cuts_distributions_20_100; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_all_cuts_distributions_100_250 = fr_energy_filtered.Filter("energy_corr_gev >= 100 && energy_corr_gev < 250")
                                                .Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions_100_250", "sumrms_flast_all_cuts_distributions_100_250; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_all_cuts_distributions_250_500 = fr_energy_filtered.Filter("energy_corr_gev >= 250 && energy_corr_gev < 500")
                                                .Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions_250_500", "sumrms_flast_all_cuts_distributions_250_500; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");     
        auto sumrms_flast_all_cuts_distributions_500_1000 = fr_energy_filtered.Filter("energy_corr_gev >= 500 && energy_corr_gev < 1000")
                                                .Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions_500_1000", "sumrms_flast_all_cuts_distributions_500_1000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_all_cuts_distributions_1000_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 1000 && energy_corr_gev < 3000")
                                                .Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions_1000_3000", "sumrms_flast_all_cuts_distributions_1000_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");
        auto sumrms_flast_all_cuts_distributions_3000 = fr_energy_filtered.Filter("energy_corr_gev >= 3000")
                                                .Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions_3000", "sumrms_flast_all_cuts_distributions_3000; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        auto sumrms_flast_all_cuts_distributions = fr_energy_filtered.Filter("evtfilter_all_cut==1")
                                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                                .Histo2D({"sumrms_flast_all_cuts_distributions", "sumrms_flast_all_cuts_distributions; sumRMS [mm]; F_{last}", static_cast<int>(sumRms_binning.size()-1), sumRms_binning.data(), static_cast<int>(flast_binning.size()-1), flast_binning.data()}, "sumRms", "fracLast", "evt_w");

        // PSD charges
        auto psd_charge_1D = fr_energy_filtered.Filter("evtfilter_psd_charge_measurement==1")
                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                .Filter("PSD_chargeX!=-999 && PSD_chargeY!=-999")
                                .Define("psdcharge", "0.5*(PSD_chargeX+PSD_chargeY)")
                                .Histo1D({"psd_charge_1D", "PSD Charge; PSD Charge; entries", 300, 0, 40}, "psdcharge", "evt_w");
        
        auto psd_charge_2D = fr_energy_filtered.Filter("evtfilter_psd_charge_measurement==1")
                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                .Filter("PSD_chargeX!=-999 && PSD_chargeY!=-999")
                                .Histo2D({"psd_charge_2D", "PSD Charge; PSD Charge X; PSD Charge Y", 400, 0, 30, 400, 0, 30}, "PSD_chargeX", "PSD_chargeY", "evt_w");

        // PSD and STK charges
        auto psd_stk_X_charge = fr_energy_filtered.Filter("evtfilter_psd_charge_measurement==1 && evtfilter_stk_charge_measurement==1")
                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                .Filter("PSD_chargeX!=-999 && STK_chargeX!=-999")
                                .Histo2D({"psd_stk_X_charge", "PSD vs STK Charge; STK X Charge; PSD X Charge",  500, 0, 100, 500, 0, 40}, "STK_chargeX", "PSD_chargeX", "evt_w");

        auto psd_stk_Y_charge = fr_energy_filtered.Filter("evtfilter_psd_charge_measurement==1 && evtfilter_stk_charge_measurement==1")
                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                .Filter("PSD_chargeY!=-999 && STK_chargeY!=-999")
                                .Histo2D({"psd_stk_Y_charge", "PSD vs STK Charge; STK Y Charge; PSD Y Charge",  500, 0, 100, 500, 0, 40}, "STK_chargeY", "PSD_chargeY", "evt_w");

        // STK charges
        auto stk_charge_1D = fr_energy_filtered.Filter("evtfilter_stk_charge_measurement==1")
                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                .Define("stkcharge", "0.5*(STK_chargeX+STK_chargeY)")
                                .Histo1D({"stk_charge_1D", "STK Charge; STK Charge; entries", 500, 0, 60}, "stkcharge", "evt_w");
        
        auto stk_charge_2D = fr_energy_filtered.Filter("evtfilter_stk_charge_measurement==1")
                                .Define("evt_w", get_weight, {"energy_corr_gev"})
                                .Filter("STK_chargeX!=-999 && STK_chargeY!=-999")
                                .Histo2D({"stk_charge_2D", "STK Charge; STK Charge X; STK Charge Y", 500, 0, 60, 500, 0, 60}, "STK_chargeX", "STK_chargeY", "evt_w");


        std::unique_ptr<TFile> output = std::make_unique<TFile>(output_file, "RECREATE");
        if (output->IsZombie()) {
            std::cerr << "Error: could not create output ROOT file [" << output_file << "]\n\n";
            exit(100);
        }

        output->mkdir("sumrms_flast_only_trigger");
        output->cd("sumrms_flast_only_trigger");

        sumrms_flast_lateral_distributions_only_trigger->Write();

        output->mkdir("sumrms_flast_after_bgofiducial");
        output->cd("sumrms_flast_after_bgofiducial");

        sumrms_flast_bgofiducial_distributions_20_100->Write();
        sumrms_flast_bgofiducial_distributions_100_250->Write();
        sumrms_flast_bgofiducial_distributions_250_500->Write();
        sumrms_flast_bgofiducial_distributions_500_1000->Write();
        sumrms_flast_bgofiducial_distributions_1000_3000->Write();
        sumrms_flast_bgofiducial_distributions_3000->Write();    
        sumrms_flast_bgofiducial_distributions->Write();

        output->mkdir("sumrms_flast_after_remove_lateral");
        output->cd("sumrms_flast_after_remove_lateral");

        sumrms_flast_lateral_distributions_20_100->Write();
        sumrms_flast_lateral_distributions_100_250->Write();
        sumrms_flast_lateral_distributions_250_500->Write();
        sumrms_flast_lateral_distributions_500_1000->Write();
        sumrms_flast_lateral_distributions_1000_3000->Write();
        sumrms_flast_lateral_distributions_3000->Write();    
        sumrms_flast_lateral_distributions->Write();

        sumrms_cosine_lateral_distributions_20_100->Write();
        sumrms_cosine_lateral_distributions_100_250->Write();
        sumrms_cosine_lateral_distributions_250_500->Write();
        sumrms_cosine_lateral_distributions_500_1000->Write();
        sumrms_cosine_lateral_distributions_1000_3000->Write();
        sumrms_cosine_lateral_distributions_3000->Write();    
        sumrms_cosine_lateral_distributions->Write();

        output->mkdir("sumrms_flast_after_remove_lateral_lowenergysumrms");
        output->cd("sumrms_flast_after_remove_lateral_lowenergysumrms");

        sumrms_flast_lateral_distributions_after_lowenergysumrms_20_100->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_100_250->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_250_500->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_500_1000->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_1000_3000->Write();
        sumrms_flast_lateral_distributions_after_lowenergysumrms_3000->Write();    
        sumrms_flast_lateral_distributions_after_lowenergysumrms->Write();

        output->mkdir("sumrms_flast_after_trackselection");
        output->cd("sumrms_flast_after_trackselection");

        sumrms_flast_trackselection_distributions_20_100->Write();
        sumrms_flast_trackselection_distributions_100_250->Write();
        sumrms_flast_trackselection_distributions_250_500->Write();
        sumrms_flast_trackselection_distributions_500_1000->Write();
        sumrms_flast_trackselection_distributions_1000_3000->Write();
        sumrms_flast_trackselection_distributions_3000->Write();    
        sumrms_flast_trackselection_distributions->Write();

        output->mkdir("sumrms_flast_after_stk1rm");
        output->cd("sumrms_flast_after_stk1rm");

        sumrms_flast_stk1rm_distributions_20_100->Write();
        sumrms_flast_stk1rm_distributions_100_250->Write();
        sumrms_flast_stk1rm_distributions_250_500->Write();
        sumrms_flast_stk1rm_distributions_500_1000->Write();
        sumrms_flast_stk1rm_distributions_1000_3000->Write();
        sumrms_flast_stk1rm_distributions_3000->Write();    
        sumrms_flast_stk1rm_distributions->Write();

        output->mkdir("sumrms_flast_after_all_cuts");
        output->cd("sumrms_flast_after_all_cuts");

        sumrms_flast_all_cuts_distributions_20_100->Write();
        sumrms_flast_all_cuts_distributions_100_250->Write();
        sumrms_flast_all_cuts_distributions_250_500->Write();
        sumrms_flast_all_cuts_distributions_500_1000->Write();
        sumrms_flast_all_cuts_distributions_1000_3000->Write();
        sumrms_flast_all_cuts_distributions_3000->Write();    
        sumrms_flast_all_cuts_distributions->Write();                                                                                              

        output->mkdir("psd_charge");
        output->cd("psd_charge");

        psd_charge_1D->Write();
        psd_charge_2D->Write();

        output->mkdir("psd_stk_charge");
        output->cd("psd_stk_charge");

        psd_stk_X_charge->Write();
        psd_stk_Y_charge->Write();

        output->mkdir("stk_charge");
        output->cd("stk_charge");

        stk_charge_1D->Write();
        stk_charge_2D->Write();

        output->Close();
    }

