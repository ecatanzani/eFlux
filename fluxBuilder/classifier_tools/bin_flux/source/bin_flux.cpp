#include "energy_config.h"
#include "list_parser.h"
#include "efficiency.h"
#include "background.h"
#include "bin_flux.h"
#include "bdtutils.h"
#include "config.h"
#include "cuts.h"

#include <memory>
#include <vector>
#include <tuple>

#include "TKey.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

inline std::tuple<double, double> getAcceptanceInBin(const char* acceptance_file, unsigned int energy_bin, const bool verbose) {
    TFile* infile {TFile::Open(acceptance_file, "READ")};
    if (infile->IsZombie()) {
        std::cerr << "\n\nError opening acceptance input file [" << acceptance_file << "]" << std::endl;
        exit(100);
    }

    TIter nextkey(infile->GetListOfKeys());
    TKey *key {nullptr};
    std::shared_ptr<TH1D> myhisto;
    while ((key=static_cast<TKey*>(nextkey())))  {
        TObject *obj {key->ReadObj()};
        if (obj->IsA()->InheritsFrom(TH1D::Class())) {
            myhisto = std::shared_ptr<TH1D>(static_cast<TH1D*>(obj));
            break;
        }
    }

    if (verbose)
        std::cout << "\nFound acceptance TH1D in input file [" << myhisto->GetName() << "] --> [" << acceptance_file << "]\n\n";

    return std::make_tuple(static_cast<double>(myhisto->GetBinContent(energy_bin)), static_cast<double>(myhisto->GetBinError(energy_bin)));
}

inline double wtsydp(const float minene, const float maxene, const float index) {
    float dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}

void bin_flux(in_args input_args) {

    std::shared_ptr<parser> list_parser = std::make_shared<parser>(input_args.input_list, input_args.verbose);
    std::shared_ptr<config> bdt_config = std::make_shared<config>(input_args.bdt_config_file);
    std::shared_ptr<energy_config> _energy_config = std::make_shared<energy_config>(input_args.energy_config_file);
    
    auto energy_binning = _energy_config->GetEnergyBinning();

    long long int total_events {list_parser->GetEvtTree()->GetEntries()};

    if (input_args.verbose) {
        bdt_config->PrintWeights();
        std::cout << "\n\nTotal number of events: " << total_events << "\n\n";
    }

    TMVA::Tools::Instance();
    auto methods_map = GetTMVAMethods(input_args.learning_method);

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
    tmva_LE_reader->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetLEWeights()).c_str());
    tmva_ME_reader->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetMEWeights()).c_str());
    tmva_HE_reader->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetHEWeights()).c_str());
    */

    tmva_reader_10_100->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetBDTWeights_10_100()).c_str());
    tmva_reader_100_250->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetBDTWeights_100_250()).c_str());
    tmva_reader_250_500->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetBDTWeights_250_500()).c_str());
    tmva_reader_500_1000->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetBDTWeights_500_1000()).c_str());
    tmva_reader_1000_3000->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetBDTWeights_1000_3000()).c_str());
    tmva_reader_3000->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetBDTWeights_3000()).c_str());

    linkTreeVariables(list_parser->GetEvtTree(), vars);

    // Extract the acceptance in the corresponding energy bin
    auto tuple_acceptance {getAcceptanceInBin(input_args.acceptance_file, input_args.energy_bin, input_args.verbose)};
    auto acceptance {std::get<0>(tuple_acceptance)};
    auto acceptance_err {std::get<1>(tuple_acceptance)};

    // Get the energy bin width
    auto energy_bin_width {_energy_config->GetEnergyBinWidth(input_args.energy_bin)};

    // Get the energy of the bin according to WTSYDP
    const float spectral_index = -3;
    auto energy_wtsydp {wtsydp(energy_binning[input_args.energy_bin-1], energy_binning[input_args.energy_bin], spectral_index)};

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

    std::unique_ptr<cuts> cl_cuts = std::make_unique<cuts>();
    auto bdt_cuts = cl_cuts->GetBDT();
    auto xtrl_cuts = cl_cuts->GetXTRL();
    auto bdt_xtrl_cuts = cl_cuts->GetBDT_XTRL();

    double tmva_classifier {-100};

    if (input_args.verbose)
        std::cout << "\n *** Computing fluxes for energy bin: " << input_args.energy_bin << " *** \n\n";

    for (unsigned int evidx=0; evidx<list_parser->GetEvtTree()->GetEntries(); ++evidx) {

        // Parse event
        list_parser->GetEvtTree()->GetEntry(evidx);
        sync_vars(vars, tmva_vars);
        tmva_classifier = -999;

        if (input_args.verbose)
            evstatus_printer(evidx, kstep, total_events);

        // Filter the event to be contained in the selected energy bin
        if (vars.energy_bin==(int)input_args.energy_bin) {

            /*
            if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100)
                tmva_classifier = tmva_LE_reader->EvaluateMVA(input_args.learning_method.c_str());
            
            else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<1000)
                tmva_classifier = tmva_ME_reader->EvaluateMVA(input_args.learning_method.c_str());
            
            else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<=10000)
                tmva_classifier = tmva_HE_reader->EvaluateMVA(input_args.learning_method.c_str());
            */

            if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100)
                tmva_classifier = tmva_reader_10_100->EvaluateMVA(input_args.learning_method.c_str());
            else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<250)
                tmva_classifier = tmva_reader_100_250->EvaluateMVA(input_args.learning_method.c_str());
            else if (vars.evt_corr_energy*gev>=250 && vars.evt_corr_energy*gev<500)
                tmva_classifier = tmva_reader_250_500->EvaluateMVA(input_args.learning_method.c_str());
            else if (vars.evt_corr_energy*gev>=500 && vars.evt_corr_energy*gev<1000)
                tmva_classifier = tmva_reader_500_1000->EvaluateMVA(input_args.learning_method.c_str());
            else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<3000)
                tmva_classifier = tmva_reader_1000_3000->EvaluateMVA(input_args.learning_method.c_str());
            else if (vars.evt_corr_energy*gev>=3000 && vars.evt_corr_energy*gev<=10000)
                tmva_classifier = tmva_reader_3000->EvaluateMVA(input_args.learning_method.c_str());

            // Loop over the bdt cuts tuples
            for (auto &cut : bdt_cuts)
                if (tmva_classifier > std::get<0>(cut))
                    std::get<1>(cut) += 1;
        
            // Loop over the xtrl cuts tuples 
            for (auto &cut : xtrl_cuts)
                if (vars.xtrl < std::get<0>(cut))
                    std::get<1>(cut) += 1;

            // Loop over the bdt_xtrl cuts tuples
            for (auto &cut : bdt_xtrl_cuts)
                if (tmva_classifier > std::get<0>(cut) && vars.xtrl < std::get<1>(cut))
                    std::get<2>(cut) += 1;
        }
    }
    
    // Calculate the signal efficiency for the given bdt cuts
    auto efficiencies = compute_efficiency(
        bdt_config,
        input_args.learning_method,
        bdt_cuts, 
        input_args.e_simu_input_list, 
        input_args.eff_corr_function, 
        input_args.energy_bin,
        energy_binning,
        input_args.verbose, 
        input_args.threads);

    // Estimate the background contamnination a given energy bin
    auto background_contamination = estimate_background(
        input_args.p_background_fit, 
        bdt_cuts, 
        input_args.energy_bin, 
        input_args.verbose, 
        input_args.threads);

    // Build the final graphs

    // BDT TGraph
    std::unique_ptr<TGraph> gr_flux_bdt                                 = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_bdt                              = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_bdt_eff_corr                        = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_bdt_eff_corr                     = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_bdt_eff_corr_b_subtracted           = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_bdt_eff_corr_b_subtracted        = std::make_unique<TGraph>();
    
    // BDT TGraphErrors
    std::unique_ptr<TGraphErrors> gr_flux_bdt_we                        = std::make_unique<TGraphErrors>();
    std::unique_ptr<TGraphErrors> gr_flux_E3_bdt_we                     = std::make_unique<TGraphErrors>();
    std::unique_ptr<TGraphErrors> gr_flux_bdt_eff_corr_we               = std::make_unique<TGraphErrors>();
    std::unique_ptr<TGraphErrors> gr_flux_E3_bdt_eff_corr_we            = std::make_unique<TGraphErrors>();

    // BDT efficiency
    std::unique_ptr<TGraph> gr_efficiency                               = std::make_unique<TGraph>();
    std::unique_ptr<TGraphErrors> gr_efficiency_we                      = std::make_unique<TGraphErrors>();

    // BDT background estimation
    std::unique_ptr<TGraph> gr_background_contamination                 = std::make_unique<TGraph>();

    // XTRL TGraph
    std::unique_ptr<TGraph> gr_flux_xtrl                                = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_xtrl                             = std::make_unique<TGraph>();

    // XTRL TGraphErrors
    std::unique_ptr<TGraphErrors> gr_flux_xtrl_we                       = std::make_unique<TGraphErrors>();
    std::unique_ptr<TGraphErrors> gr_flux_E3_xtrl_we                    = std::make_unique<TGraphErrors>();

    // BDT XTRL TGraph
    std::unique_ptr<TGraph2D> gr_flux_bdt_xtrl                          = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_E3_bdt_xtrl                       = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_bdt_xtrl_eff_corr                 = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_E3_bdt_xtrl_eff_corr              = std::make_unique<TGraph2D>();
    
    
    // Fill graphs
    for (size_t idx = 0; idx < bdt_cuts.size(); ++idx) {

        // BDT

        // Build the flux
        double flux             = std::get<1>(bdt_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        double flux_err         = sqrt(
                                        pow(sqrt(std::get<1>(bdt_cuts[idx]))/input_args.exposure*acceptance*energy_bin_width, 2) + 
                                        pow(std::get<1>(bdt_cuts[idx])*acceptance_err/(input_args.exposure*energy_bin_width*pow(acceptance, 2)), 2)
                                );
        double flux_E3          = flux*pow(energy_wtsydp, 3);
        double flux_E3_err      = flux_err*pow(energy_wtsydp, 3);

        // Build the flux subtracting the background
        double flux_b_sub       = (std::get<1>(bdt_cuts[idx]) - std::get<1>(background_contamination[idx]))/(input_args.exposure*acceptance*energy_bin_width);
        
        gr_flux_bdt             ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux);
        gr_flux_bdt_we          ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux);
        gr_flux_bdt_we          ->SetPointError((int)idx, 0, flux_err);
        gr_flux_E3_bdt          ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3);
        gr_flux_E3_bdt_we       ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3);
        gr_flux_E3_bdt_we       ->SetPointError((int)idx, 0, flux_E3_err);
        
        double eff              = std::get<1>(efficiencies[idx]);
        double eff_err          = std::get<2>(efficiencies[idx]);

        if (eff) {
            double flux_ec              = flux/eff;
            double flux_E3_ec           = flux_ec*pow(energy_wtsydp, 3);

            double flux_ec_b_sub        = flux_b_sub/eff;
            double flux_E3_ec_b_sub     = flux_ec_b_sub*pow(energy_wtsydp, 3); 

            flux_err                    = sqrt(
                                                pow(sqrt(std::get<1>(bdt_cuts[idx]))/input_args.exposure*acceptance*energy_bin_width*eff, 2) + 
                                                pow(std::get<1>(bdt_cuts[idx])*acceptance_err/(input_args.exposure*energy_bin_width*eff*pow(acceptance, 2)), 2) +
                                                pow(std::get<1>(bdt_cuts[idx])*eff_err/(input_args.exposure*energy_bin_width*acceptance*pow(eff, 2)), 2)
                                        );
            flux_E3_err                 = flux_err*pow(energy_wtsydp, 3);
            
            gr_flux_bdt_eff_corr                        ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_ec);
            
            gr_flux_bdt_eff_corr_we                     ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_ec);
            gr_flux_bdt_eff_corr_we                     ->SetPointError((int)idx, 0, flux_err);
            
            gr_flux_E3_bdt_eff_corr                     ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3_ec);

            gr_flux_E3_bdt_eff_corr_we                  ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3_ec);
            gr_flux_E3_bdt_eff_corr_we                  ->SetPointError((int)idx, 0, flux_E3_err);

            gr_flux_bdt_eff_corr_b_subtracted           ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_ec_b_sub);
            gr_flux_E3_bdt_eff_corr_b_subtracted        ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3_ec_b_sub);
        }

        gr_efficiency                   ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), eff);
        gr_efficiency_we                ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), eff);
        gr_efficiency_we                ->SetPointError((int)idx, 0, eff_err);

        gr_background_contamination     ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), std::get<1>(background_contamination[idx]));

        // XTRL
        flux                            = std::get<1>(xtrl_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        flux_E3                         = flux*pow(energy_wtsydp, 3);
        flux_err                        = sqrt(
                                                pow(sqrt(std::get<1>(xtrl_cuts[idx]))/input_args.exposure*acceptance*energy_bin_width, 2) + 
                                                pow(std::get<1>(xtrl_cuts[idx])*acceptance_err/(input_args.exposure*energy_bin_width*pow(acceptance, 2)), 2)
                                        );
        flux_E3_err                     = flux_err*pow(energy_wtsydp, 3);

        gr_flux_xtrl                    ->SetPoint((int)idx, std::get<0>(xtrl_cuts[idx]), flux);
        gr_flux_xtrl_we                 ->SetPoint((int)idx, std::get<0>(xtrl_cuts[idx]), flux);
        gr_flux_xtrl_we                 ->SetPointError((int)idx, 0, flux_err);
        gr_flux_E3_xtrl                 ->SetPoint((int)idx, std::get<0>(xtrl_cuts[idx]), flux_E3); 
        gr_flux_E3_xtrl_we              ->SetPoint((int)idx, std::get<0>(xtrl_cuts[idx]), flux_E3);
        gr_flux_E3_xtrl_we              ->SetPointError((int)idx, 0, flux_E3_err);   
    }
    
    for (size_t idx = 0; idx < bdt_xtrl_cuts.size(); ++idx) {

        double flux             = std::get<2>(bdt_xtrl_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        double flux_E3          = flux*pow(energy_wtsydp, 3);

        gr_flux_bdt_xtrl                        ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux);
        gr_flux_E3_bdt_xtrl                     ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_E3);

        // Get the efficiencies corresponding to the BDT cut
        double eff {0};
        for (auto&& eff_tuple : efficiencies)
            if (std::get<0>(bdt_xtrl_cuts[idx]) == std::get<0>(eff_tuple)) {
                eff = std::get<1>(eff_tuple);
                break;
            }

        if (eff) {
            double flux_ec        = flux/eff;
            double flux_E3_ec     = flux_ec*pow(energy_wtsydp, 3);

            gr_flux_bdt_xtrl_eff_corr       ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_ec);
            gr_flux_E3_bdt_xtrl_eff_corr    ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_E3_ec);
        }
    }

    // Set graphs properties
    gr_flux_bdt->SetName((std::string("gr_flux_bdt_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_we->SetName((std::string("gr_flux_bdt_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_we->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_we->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_we->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_eff_corr->SetName((std::string("gr_flux_bdt_eff_corr_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_eff_corr->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_eff_corr_b_subtracted->SetName((std::string("gr_flux_bdt_eff_corr_b_subtracted_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_b_subtracted->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_b_subtracted->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_eff_corr_b_subtracted->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_eff_corr_we->SetName((std::string("gr_flux_bdt_eff_corr_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_we->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_we->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_eff_corr_we->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");
    
    gr_flux_xtrl->SetName((std::string("gr_flux_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_xtrl->GetXaxis()->SetTitle("XTRL");
    gr_flux_xtrl->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_xtrl_we->SetName((std::string("gr_flux_xtrl_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_xtrl_we->SetTitle((std::string("All-Electron Flux vs. XTRL cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_xtrl_we->GetXaxis()->SetTitle("XTRL");
    gr_flux_xtrl_we->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");
    
    gr_flux_E3_bdt->SetName((std::string("gr_flux_E3_bdt_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_we->SetName((std::string("gr_flux_E3_bdt_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_we->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_we->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_we->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_eff_corr->SetName((std::string("gr_flux_E3_bdt_eff_corr_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_eff_corr->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_eff_corr_we->SetName((std::string("gr_flux_E3_bdt_eff_corr_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_we->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_we->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_eff_corr_we->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_eff_corr_b_subtracted->SetName((std::string("gr_flux_E3_bdt_eff_corr_b_subtracted_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_b_subtracted->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_b_subtracted->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_eff_corr_b_subtracted->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_xtrl->SetName((std::string("gr_flux_E3_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_xtrl->GetXaxis()->SetTitle("XTRL");
    gr_flux_E3_xtrl->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_xtrl_we->SetName((std::string("gr_flux_E3_xtrl_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_xtrl_we->SetTitle((std::string("All-Electron Flux vs. XTRL cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_xtrl_we->GetXaxis()->SetTitle("XTRL");
    gr_flux_E3_xtrl_we->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_xtrl->SetName((std::string("gr_flux_bdt_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_xtrl->GetYaxis()->SetTitle("XTRL");
    gr_flux_bdt_xtrl->GetZaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_xtrl_eff_corr->SetName((std::string("gr_flux_bdt_xtrl_eff_corr_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl_eff_corr->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl_eff_corr->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_xtrl_eff_corr->GetYaxis()->SetTitle("XTRL");
    gr_flux_bdt_xtrl_eff_corr->GetZaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_xtrl->SetName((std::string("gr_flux_E3_bdt_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_xtrl->GetYaxis()->SetTitle("XTRL");
    gr_flux_E3_bdt_xtrl->GetZaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_xtrl_eff_corr->SetName((std::string("gr_flux_E3_bdt_xtrl_eff_corr_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl_eff_corr->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl_eff_corr->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_xtrl_eff_corr->GetYaxis()->SetTitle("XTRL");
    gr_flux_E3_bdt_xtrl_eff_corr->GetZaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_efficiency->SetName((std::string("gr_efficiency_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency->SetTitle((std::string("BDT efficiency for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency->GetXaxis()->SetTitle("BDT");
    gr_efficiency->GetYaxis()->SetTitle("Efficiency");

    gr_efficiency_we->SetName((std::string("gr_efficiency_we_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency_we->SetTitle((std::string("BDT efficiency for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency_we->GetXaxis()->SetTitle("BDT");
    gr_efficiency_we->GetYaxis()->SetTitle("Efficiency");

    gr_background_contamination->SetName((std::string("gr_background_contamination_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_background_contamination->SetTitle((std::string("BDT background contamination for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_background_contamination->GetXaxis()->SetTitle("BDT");
    gr_background_contamination->GetYaxis()->SetTitle("Background estimation");

    // Write to output file
    TFile *outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (outfile->IsZombie()) {
        std::cerr << "Error writing output file " << input_args.output_path << std::endl;
        exit(1);
    }

    std::string out_dir_name = "energybin_" + std::to_string(input_args.energy_bin);
    outfile->mkdir(out_dir_name.c_str());
    outfile->cd(out_dir_name.c_str());

    gr_flux_bdt                             ->Write();
    gr_flux_E3_bdt                          ->Write();
    gr_flux_bdt_eff_corr                    ->Write();
    gr_flux_E3_bdt_eff_corr                 ->Write();
    gr_flux_bdt_eff_corr_b_subtracted       ->Write();
    gr_flux_E3_bdt_eff_corr_b_subtracted    ->Write();

    gr_flux_bdt_we                          ->Write();
    gr_flux_E3_bdt_we                       ->Write();
    gr_flux_bdt_eff_corr_we                 ->Write();
    gr_flux_E3_bdt_eff_corr_we              ->Write();
    
    gr_efficiency                           ->Write();
    gr_efficiency_we                        ->Write();

    gr_background_contamination             ->Write();

    gr_flux_xtrl                            ->Write();
    gr_flux_E3_xtrl                         ->Write();

    gr_flux_xtrl_we                         ->Write();  
    gr_flux_E3_xtrl_we                      ->Write();

    gr_flux_bdt_xtrl                        ->Write();
    gr_flux_bdt_xtrl_eff_corr               ->Write();
    gr_flux_E3_bdt_xtrl                     ->Write();
    gr_flux_E3_bdt_xtrl_eff_corr            ->Write();
    
    outfile->Close();
}