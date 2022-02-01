#include "config.h"
#include "bdtutils.h"
#include "efficiency.h"
#include "list_parser.h"

#include <memory>

#include "TF1.h"
#include "TFile.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

inline std::tuple<std::shared_ptr<TF1>, std::shared_ptr<TF1>> extractCorrectionFunctions(
    const char* file_name, 
    const char* tf1_shift_name="shift_fit_function",
    const char* tf1_sigma_name="sigma_ratio_fit_function") {

        TFile *infile = TFile::Open(file_name, "READ");
        if (infile->IsZombie()) {
            std::cout << "Error opening input TF1 file: [" << file_name << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TF1> tf1_shift = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_shift_name)));
        std::shared_ptr<TF1> tf1_sigma = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_sigma_name)));

        return std::tuple<std::shared_ptr<TF1>, std::shared_ptr<TF1>>(tf1_shift, tf1_sigma);
    }

const std::tuple<std::vector<std::tuple<double, double>>, std::vector<std::tuple<double, double>>> compute_efficiency(
    std::shared_ptr<config> bdt_config,
    const std::string learning_method,
    const std::vector<std::tuple<double, unsigned int>> bdt_cuts,
    const std::string mc_file_list,
    const char* mc_correction_function,
    const unsigned int energy_bin,
    std::vector<float>& energy_binning,
    const bool verbose,
    const unsigned int threads) {

        if (verbose)
            std::cout << "\n\n*** Efficiency Facility***\n\n";

        // Parse MC file list
        std::shared_ptr<parser> mc_list_parser = std::make_unique<parser>(mc_file_list, verbose);
        
        long long int total_events {mc_list_parser->GetEvtTree()->GetEntries()};
        
        if (verbose)
            std::cout << "\n\nTotal number of events: " << total_events << "\n\n";

        TMVA::Tools::Instance();
        auto methods_map = GetTMVAMethods(learning_method);

        std::shared_ptr<TMVA::Reader> tmva_LE_reader = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_ME_reader = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_HE_reader = std::make_shared<TMVA::Reader>();

        // Declare BDT and DATA variables
        bdt_vars tmva_vars;
        data_vars vars;

        // Attach variables to reader
        addVariableToReader(tmva_LE_reader, tmva_vars);
        addVariableToReader(tmva_ME_reader, tmva_vars);
        addVariableToReader(tmva_HE_reader, tmva_vars);

        tmva_LE_reader->BookMVA(learning_method, (bdt_config->GetLEWeights()).c_str());
        tmva_ME_reader->BookMVA(learning_method, (bdt_config->GetMEWeights()).c_str());
        tmva_HE_reader->BookMVA(learning_method, (bdt_config->GetHEWeights()).c_str());

        // Lionk variables to the MC TChain
        double simu_energy {0};
        linkTreeVariables(mc_list_parser->GetEvtTree(), vars);
        mc_list_parser->GetEvtTree()->SetBranchAddress("simu_energy", &simu_energy);

        // Weight MC events
        const double signal_spectral_index = -3;
        auto get_weight = [signal_spectral_index, &energy_binning] (const double energy_gev) -> double {
            return std::pow(energy_gev, signal_spectral_index+1)*std::pow(energy_binning[0], fabs(signal_spectral_index+1));
        };

        // Build the efficiency tuple
        std::vector<std::tuple<double, double>> efficiency_values_method1 (bdt_cuts.size());
        std::vector<std::tuple<double, double>> efficiency_values_method2 (bdt_cuts.size());

        for (size_t idx=0; idx<efficiency_values_method1.size(); ++idx) {
            efficiency_values_method1[idx] = std::make_tuple(std::get<0>(bdt_cuts[idx]), 0);
            efficiency_values_method2[idx] = std::make_tuple(std::get<0>(bdt_cuts[idx]), 0);
        }

        // Extract the efficiency correction functions
        std::shared_ptr<TF1> shift_function, sigma_function;
        std::tie(shift_function, sigma_function) = extractCorrectionFunctions(mc_correction_function);

        // Loop on the events
        double gev {0.001};
        int kstep {10};
        
        auto evstatus_printer = [](unsigned int evidx, int &kstep, long long int total_events) {
            auto percentage = ((evidx + 1) / (double)total_events) * 100;
            if (floor(percentage) != 0 && ((int)floor(percentage) % kstep) == 0) {
                std::cout << "\n" << (int)percentage << " %\t | \tProcessed " << evidx + 1 << " events / " << total_events;
                kstep += 10;
            }
        };

        auto get_bdt_cut_method1 = [shift_function, sigma_function] (const double cut, const double energy_gev) -> double {
            return cut - shift_function->Eval(energy_gev)/sigma_function->Eval(energy_gev);
        };
        auto get_bdt_cut_method2 = [shift_function, sigma_function] (const double cut, const double energy_gev) -> double {
            return cut - shift_function->Eval(energy_gev)*sigma_function->Eval(energy_gev);
        };

        double n_total {0};
        double tmva_classifier {-100};

        if (verbose)
            std::cout << "\n\n *** Computing efficiency for energy bin: " << energy_bin << " *** \n\n";

        for (unsigned int evidx=0; evidx<mc_list_parser->GetEvtTree()->GetEntries(); ++evidx) {
            
            // Parse event
            mc_list_parser->GetEvtTree()->GetEntry(evidx);
            sync_vars(vars, tmva_vars);
            tmva_classifier = -999;

            if (verbose)
                evstatus_printer(evidx, kstep, total_events);

            // Filter the event to be contained in the selected energy bin
            if (vars.energy_bin == (int)energy_bin) {
                
                if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100) 
                    tmva_classifier = tmva_LE_reader->EvaluateMVA(learning_method.c_str());
                
                else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<1000)
                    tmva_classifier = tmva_ME_reader->EvaluateMVA(learning_method.c_str());
                
                else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<=10000)
                    tmva_classifier = tmva_HE_reader->EvaluateMVA(learning_method.c_str());

                for (auto &cut : efficiency_values_method1)
                    if (tmva_classifier < get_bdt_cut_method1(std::get<0>(cut), vars.evt_corr_energy*gev))
                        std::get<1>(cut) += get_weight(vars.evt_corr_energy*gev);
                
                for (auto &cut : efficiency_values_method2)
                    if (tmva_classifier < get_bdt_cut_method2(std::get<0>(cut), vars.evt_corr_energy*gev))
                        std::get<1>(cut) += get_weight(vars.evt_corr_energy*gev);
            
                n_total += get_weight(vars.evt_corr_energy*gev);
            }
        }

        for (size_t idx=0; idx<efficiency_values_method1.size(); ++idx) {
            std::get<1>(efficiency_values_method1[idx]) /= n_total;
            std::get<1>(efficiency_values_method2[idx]) /= n_total;
        }

        return std::make_tuple(efficiency_values_method1, efficiency_values_method2);
    }