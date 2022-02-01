#include "energy_config.h"
#include "list_parser.h"
#include "efficiency.h"
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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

inline double getAcceptanceInBin(const char* acceptance_file, unsigned int energy_bin, const bool verbose) {
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

    return static_cast<double>(myhisto->GetBinContent(energy_bin));
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

    tmva_LE_reader->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetLEWeights()).c_str());
    tmva_ME_reader->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetMEWeights()).c_str());
    tmva_HE_reader->BookMVA(input_args.learning_method.c_str(), (bdt_config->GetHEWeights()).c_str());
    
    linkTreeVariables(list_parser->GetEvtTree(), vars);

    // Extract the acceptance in the corresponding energy bin
    auto acceptance {getAcceptanceInBin(input_args.acceptance_file, input_args.energy_bin, input_args.verbose)};

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

            if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100)
                tmva_classifier = tmva_LE_reader->EvaluateMVA(input_args.learning_method.c_str());
            
            else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<1000)
                tmva_classifier = tmva_ME_reader->EvaluateMVA(input_args.learning_method.c_str());
            
            else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<=10000)
                tmva_classifier = tmva_HE_reader->EvaluateMVA(input_args.learning_method.c_str());

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
        input_args.simu_input_list, 
        input_args.eff_corr_function, 
        input_args.energy_bin,
        energy_binning,
        input_args.verbose, 
        input_args.threads);

    // Build the final graphs
    std::unique_ptr<TGraph> gr_flux_bdt                                 = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_bdt_eff_corr_method1                = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_bdt_eff_corr_method2                = std::make_unique<TGraph>();

    std::unique_ptr<TGraph> gr_flux_xtrl                                = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_bdt                              = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_bdt_eff_corr_method1             = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_bdt_eff_corr_method2             = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_flux_E3_xtrl                             = std::make_unique<TGraph>();

    std::unique_ptr<TGraph2D> gr_flux_bdt_xtrl                          = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_bdt_xtrl_eff_corr_method1         = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_bdt_xtrl_eff_corr_method2         = std::make_unique<TGraph2D>();

    std::unique_ptr<TGraph2D> gr_flux_E3_bdt_xtrl                       = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_E3_bdt_xtrl_eff_corr_method1      = std::make_unique<TGraph2D>();
    std::unique_ptr<TGraph2D> gr_flux_E3_bdt_xtrl_eff_corr_method2      = std::make_unique<TGraph2D>();

    std::unique_ptr<TGraph> gr_efficiency_method1                       = std::make_unique<TGraph>();
    std::unique_ptr<TGraph> gr_efficiency_method2                       = std::make_unique<TGraph>();
    
    // Fill graphs
    for (size_t idx = 0; idx < bdt_cuts.size(); ++idx) {

        // BDT
        double flux             = std::get<1>(bdt_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        double flux_E3          = flux*pow(energy_wtsydp, 3);
        
        gr_flux_bdt                         ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux);
        gr_flux_E3_bdt                      ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3);
        
        double eff_method1      = 1-std::get<1>((std::get<0>(efficiencies))[idx]);
        double eff_method2      = 1-std::get<1>((std::get<1>(efficiencies))[idx]); 

        if (eff_method1) {
            double flux_ec_1        = flux/eff_method1;
            double flux_E3_ec_1     = flux_ec_1*pow(energy_wtsydp, 3);
            gr_flux_bdt_eff_corr_method1        ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_ec_1);
            gr_flux_E3_bdt_eff_corr_method1     ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3_ec_1);
        }

        if (eff_method2) {
            double flux_ec_2        = flux/eff_method2;
            double flux_E3_ec_2     = flux_ec_2*pow(energy_wtsydp, 3);
            gr_flux_bdt_eff_corr_method2        ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_ec_2);
            gr_flux_E3_bdt_eff_corr_method2     ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), flux_E3_ec_2);
        }

        gr_efficiency_method1                   ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), eff_method1);
        gr_efficiency_method2                   ->SetPoint((int)idx, std::get<0>(bdt_cuts[idx]), eff_method2);

        // XTRL
        flux                    = std::get<1>(xtrl_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        flux_E3                 = flux*pow(energy_wtsydp, 3);

        gr_flux_xtrl                        ->SetPoint((int)idx, std::get<0>(xtrl_cuts[idx]), flux);
        gr_flux_E3_xtrl                     ->SetPoint((int)idx, std::get<0>(xtrl_cuts[idx]), flux_E3);    
    }
    
    for (size_t idx = 0; idx < bdt_xtrl_cuts.size(); ++idx) {

        double flux             = std::get<2>(bdt_xtrl_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        double flux_E3          = flux*pow(energy_wtsydp, 3);

        gr_flux_bdt_xtrl                        ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux);
        gr_flux_E3_bdt_xtrl                     ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_E3);

        // Get the efficiencies corresponding to the BDT cut
        double eff_method1 {0}, eff_method2 {0};
        for (auto&& eff : std::get<0>(efficiencies))
            if (std::get<0>(bdt_xtrl_cuts[idx]) == std::get<0>(eff)) {
                eff_method1 = 1-std::get<1>(eff);
                break;
            }

        for (auto&& eff : std::get<1>(efficiencies))
            if (std::get<0>(bdt_xtrl_cuts[idx]) == std::get<0>(eff)) {
                eff_method2 = 1-std::get<1>(eff);
                break;
            }

        if (eff_method1) {
            double flux_ec_1        = flux/eff_method1;
            double flux_E3_ec_1     = flux_ec_1*pow(energy_wtsydp, 3);

            gr_flux_bdt_xtrl_eff_corr_method1       ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_ec_1);
            gr_flux_E3_bdt_xtrl_eff_corr_method1    ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_E3_ec_1);
        }

        if (eff_method2) {
            double flux_ec_2        = flux/eff_method2;
            double flux_E3_ec_2     = flux_ec_2*pow(energy_wtsydp, 3);
            gr_flux_bdt_xtrl_eff_corr_method2       ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_ec_2);
            gr_flux_E3_bdt_xtrl_eff_corr_method2    ->SetPoint((int)idx, std::get<0>(bdt_xtrl_cuts[idx]), std::get<1>(bdt_xtrl_cuts[idx]), flux_E3_ec_2);
        }
    }

    // Set graphs properties
    gr_flux_bdt->SetName((std::string("gr_flux_bdt_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_eff_corr_method1->SetName((std::string("gr_flux_bdt_eff_corr_method1_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_method1->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_method1->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_eff_corr_method1->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_eff_corr_method2->SetName((std::string("gr_flux_bdt_eff_corr_method2_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_method2->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_eff_corr_method2->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_eff_corr_method2->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");
    
    gr_flux_xtrl->SetName((std::string("gr_flux_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_xtrl->GetXaxis()->SetTitle("XTRL");
    gr_flux_xtrl->GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");
    
    gr_flux_E3_bdt->SetName((std::string("gr_flux_E3_bdt_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_eff_corr_method1->SetName((std::string("gr_flux_E3_bdt_eff_corr_method1_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_method1->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_method1->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_eff_corr_method1->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_eff_corr_method2->SetName((std::string("gr_flux_E3_bdt_eff_corr_method2_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_method2->SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_eff_corr_method2->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_eff_corr_method2->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");
    
    gr_flux_E3_xtrl->SetName((std::string("gr_flux_E3_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_xtrl->GetXaxis()->SetTitle("XTRL");
    gr_flux_E3_xtrl->GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_xtrl->SetName((std::string("gr_flux_bdt_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_xtrl->GetYaxis()->SetTitle("XTRL");
    gr_flux_bdt_xtrl->GetZaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_xtrl_eff_corr_method1->SetName((std::string("gr_flux_bdt_xtrl_eff_corr_method1_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl_eff_corr_method1->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl_eff_corr_method1->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_xtrl_eff_corr_method1->GetYaxis()->SetTitle("XTRL");
    gr_flux_bdt_xtrl_eff_corr_method1->GetZaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_bdt_xtrl_eff_corr_method2->SetName((std::string("gr_flux_bdt_xtrl_eff_corr_method2_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl_eff_corr_method2->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_bdt_xtrl_eff_corr_method2->GetXaxis()->SetTitle("BDT");
    gr_flux_bdt_xtrl_eff_corr_method2->GetYaxis()->SetTitle("XTRL");
    gr_flux_bdt_xtrl_eff_corr_method2->GetZaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_xtrl->SetName((std::string("gr_flux_E3_bdt_xtrl_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_xtrl->GetYaxis()->SetTitle("XTRL");
    gr_flux_E3_bdt_xtrl->GetZaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_xtrl_eff_corr_method1->SetName((std::string("gr_flux_E3_bdt_xtrl_eff_corr_method1_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl_eff_corr_method1->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl_eff_corr_method1->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_xtrl_eff_corr_method1->GetYaxis()->SetTitle("XTRL");
    gr_flux_E3_bdt_xtrl_eff_corr_method1->GetZaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_flux_E3_bdt_xtrl_eff_corr_method2->SetName((std::string("gr_flux_E3_bdt_xtrl_eff_corr_method2_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl_eff_corr_method2->SetTitle((std::string("All-Electron Flux vs. XTRL and BDT cuts for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3_bdt_xtrl_eff_corr_method2->GetXaxis()->SetTitle("BDT");
    gr_flux_E3_bdt_xtrl_eff_corr_method2->GetYaxis()->SetTitle("XTRL");
    gr_flux_E3_bdt_xtrl_eff_corr_method2->GetZaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    gr_efficiency_method1->SetName((std::string("gr_efficiency_method1_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency_method1->SetTitle((std::string("BDT efficiency for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency_method1->GetXaxis()->SetTitle("BDT");
    gr_efficiency_method1->GetYaxis()->SetTitle("Efficiency");

    gr_efficiency_method2->SetName((std::string("gr_efficiency_method2_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency_method2->SetTitle((std::string("BDT efficiency for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_efficiency_method2->GetXaxis()->SetTitle("BDT");
    gr_efficiency_method2->GetYaxis()->SetTitle("Efficiency");

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
    gr_flux_bdt_eff_corr_method1            ->Write();
    gr_flux_bdt_eff_corr_method2            ->Write();

    gr_flux_xtrl                            ->Write();
    gr_flux_E3_bdt                          ->Write();
    gr_flux_E3_bdt_eff_corr_method1         ->Write();
    gr_flux_E3_bdt_eff_corr_method2         ->Write();
    gr_flux_E3_xtrl                         ->Write();

    gr_flux_bdt_xtrl                        ->Write();
    gr_flux_bdt_xtrl_eff_corr_method1       ->Write();
    gr_flux_bdt_xtrl_eff_corr_method2       ->Write();
    gr_flux_E3_bdt_xtrl                     ->Write();
    gr_flux_E3_bdt_xtrl_eff_corr_method1    ->Write();
    gr_flux_E3_bdt_xtrl_eff_corr_method2    ->Write();

    gr_efficiency_method1                   ->Write();
    gr_efficiency_method2                   ->Write();
    
    outfile->Close();
}