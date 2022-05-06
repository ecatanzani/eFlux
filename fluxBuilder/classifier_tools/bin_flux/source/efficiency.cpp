#include "config.h"
#include "bdtutils.h"
#include "efficiency.h"
#include "list_parser.h"

#include <memory>

#include "TF1.h"
#include "TFile.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TGraphErrors.h"

inline std::tuple<std::tuple<std::shared_ptr<TF1>, std::shared_ptr<TF1>>, std::tuple<std::shared_ptr<TGraphErrors>, std::shared_ptr<TGraphErrors>>> extractCorrectionFunctions(
    const char* file_name, 
    const char* tf1_shift_name="shift_fit_function",
    const char* tf1_sigma_name="sigma_ratio_fit_function",
    const char* gr_shift_name="gr_bin_shift",
    const char* gr_sigma_ratio_name="gr_bin_ratio") {

        TFile *infile = TFile::Open(file_name, "READ");
        if (infile->IsZombie()) {
            std::cout << "Error opening input TF1 file: [" << file_name << "]\n\n";
            exit(100);
        }

        std::shared_ptr<TF1> tf1_shift = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_shift_name)));
        std::shared_ptr<TF1> tf1_sigma = std::shared_ptr<TF1>(static_cast<TF1*>(infile->Get(tf1_sigma_name)));
        std::shared_ptr<TGraphErrors> gr_shift = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(gr_shift_name)));
        std::shared_ptr<TGraphErrors> gr_sigma = std::shared_ptr<TGraphErrors>(static_cast<TGraphErrors*>(infile->Get(gr_sigma_ratio_name)));

        return std::make_tuple(std::make_tuple(tf1_shift, tf1_sigma), std::make_tuple(gr_shift, gr_sigma));
    }

const std::vector<std::tuple<double, double, double>> compute_efficiency(
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

        /*
        std::shared_ptr<TMVA::Reader> tmva_LE_reader = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_ME_reader = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_HE_reader = std::make_shared<TMVA::Reader>();
        */

        std::shared_ptr<TMVA::Reader> tmva_reader_10_100 = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_reader_100_250 = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_reader_250_500 = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_reader_500_1000 = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_reader_1000_3000 = std::make_shared<TMVA::Reader>();
        std::shared_ptr<TMVA::Reader> tmva_reader_3000 = std::make_shared<TMVA::Reader>();

        // Declare BDT and DATA variables
        bdt_vars tmva_vars;
        data_vars vars;

        // Attach variables to reader
        /*
        addVariableToReader(tmva_LE_reader, tmva_vars);
        addVariableToReader(tmva_ME_reader, tmva_vars);
        addVariableToReader(tmva_HE_reader, tmva_vars);
        */

        addVariableToReader(tmva_reader_10_100, tmva_vars);
        addVariableToReader(tmva_reader_100_250, tmva_vars);
        addVariableToReader(tmva_reader_250_500, tmva_vars);
        addVariableToReader(tmva_reader_500_1000, tmva_vars);
        addVariableToReader(tmva_reader_1000_3000, tmva_vars);
        addVariableToReader(tmva_reader_3000, tmva_vars);

        /*
        tmva_LE_reader->BookMVA(learning_method, (bdt_config->GetLEWeights()).c_str());
        tmva_ME_reader->BookMVA(learning_method, (bdt_config->GetMEWeights()).c_str());
        tmva_HE_reader->BookMVA(learning_method, (bdt_config->GetHEWeights()).c_str());
        */

        tmva_reader_10_100->BookMVA(learning_method.c_str(), (bdt_config->GetBDTWeights_10_100()).c_str());
        tmva_reader_100_250->BookMVA(learning_method.c_str(), (bdt_config->GetBDTWeights_100_250()).c_str());
        tmva_reader_250_500->BookMVA(learning_method.c_str(), (bdt_config->GetBDTWeights_250_500()).c_str());
        tmva_reader_500_1000->BookMVA(learning_method.c_str(), (bdt_config->GetBDTWeights_500_1000()).c_str());
        tmva_reader_1000_3000->BookMVA(learning_method.c_str(), (bdt_config->GetBDTWeights_1000_3000()).c_str());
        tmva_reader_3000->BookMVA(learning_method.c_str(), (bdt_config->GetBDTWeights_3000()).c_str());

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
        std::vector<std::tuple<double, double, double>> efficiency_values (bdt_cuts.size());
        
        for (size_t idx=0; idx<efficiency_values.size(); ++idx)
            efficiency_values[idx] = std::make_tuple(std::get<0>(bdt_cuts[idx]), 0, 0);

        // Extract the efficiency correction functions
        auto corrections = extractCorrectionFunctions(mc_correction_function);

        auto shift_function = std::get<0>(std::get<0>(corrections));
        auto sigma_function = std::get<1>(std::get<0>(corrections));
        auto gr_shift = std::get<0>(std::get<1>(corrections));
        auto gr_sigma = std::get<1>(std::get<1>(corrections));

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

        auto get_bdt_cut = [shift_function, sigma_function, gr_shift, gr_sigma] (const double cut, const double energy_gev, const int energy_bin) -> double {
            double cut_corr {0};
            if (energy_bin<=33)
                cut_corr = (cut - gr_shift->GetPointY(energy_bin-1))/gr_sigma->GetPointY(energy_bin-1); 
            else
                cut_corr = (cut - shift_function->Eval(energy_gev))/sigma_function->Eval(energy_gev);
            return cut_corr;
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
                
                /*
                if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100) 
                    tmva_classifier = tmva_LE_reader->EvaluateMVA(learning_method.c_str());
                
                else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<1000)
                    tmva_classifier = tmva_ME_reader->EvaluateMVA(learning_method.c_str());
                
                else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<=10000)
                    tmva_classifier = tmva_HE_reader->EvaluateMVA(learning_method.c_str());
                */

                if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100)
                    tmva_classifier = tmva_reader_10_100->EvaluateMVA(learning_method.c_str());
                else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<250)
                    tmva_classifier = tmva_reader_100_250->EvaluateMVA(learning_method.c_str());
                else if (vars.evt_corr_energy*gev>=250 && vars.evt_corr_energy*gev<500)
                    tmva_classifier = tmva_reader_250_500->EvaluateMVA(learning_method.c_str());
                else if (vars.evt_corr_energy*gev>=500 && vars.evt_corr_energy*gev<1000)
                    tmva_classifier = tmva_reader_500_1000->EvaluateMVA(learning_method.c_str());
                else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<3000)
                    tmva_classifier = tmva_reader_1000_3000->EvaluateMVA(learning_method.c_str());
                else if (vars.evt_corr_energy*gev>=3000 && vars.evt_corr_energy*gev<=10000)
                    tmva_classifier = tmva_reader_3000->EvaluateMVA(learning_method.c_str());

                for (auto &cut : efficiency_values)
                    if (tmva_classifier < get_bdt_cut(std::get<0>(cut), vars.evt_corr_energy*gev, vars.energy_bin))
                        std::get<1>(cut) += get_weight(vars.evt_corr_energy*gev);
            
                n_total += get_weight(vars.evt_corr_energy*gev);
            }
        }

        for (auto& eff_tuple : efficiency_values) {
            auto success_try = n_total - std::get<1>(eff_tuple);
            auto fail_try = std::get<1>(eff_tuple);

            // Compute the error with falt bayesian prior
            std::get<2>(eff_tuple) = sqrt((success_try+1)*(fail_try+1)/pow(success_try+fail_try+2, 3));
            // Compute the efficiency
            std::get<1>(eff_tuple) = 1 - std::get<1>(eff_tuple)/n_total;
        }

        return efficiency_values;
    }