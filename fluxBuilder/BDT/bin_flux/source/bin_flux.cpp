#include "energy_config.h"
#include "list_parser.h"
#include "bdtutils.h"
#include "bin_flux.h"
#include "config.h"

#include <memory>
#include <vector>
#include <tuple>

#include "TKey.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

inline void addVariableToReader(std::shared_ptr<TMVA::Reader> tmva_reader, bdt_vars &tmva_vars) {
    tmva_reader->AddVariable("rmslayer_norm_1", &tmva_vars.rmslayer_norm_1);
    tmva_reader->AddVariable("rmslayer_norm_2", &tmva_vars.rmslayer_norm_2);
    tmva_reader->AddVariable("rmslayer_norm_3", &tmva_vars.rmslayer_norm_3);
    tmva_reader->AddVariable("rmslayer_norm_4", &tmva_vars.rmslayer_norm_4);
    tmva_reader->AddVariable("rmslayer_norm_5", &tmva_vars.rmslayer_norm_5);
    tmva_reader->AddVariable("rmslayer_norm_6", &tmva_vars.rmslayer_norm_6);
    tmva_reader->AddVariable("rmslayer_norm_7", &tmva_vars.rmslayer_norm_7);
    tmva_reader->AddVariable("rmslayer_norm_8", &tmva_vars.rmslayer_norm_8);
    tmva_reader->AddVariable("rmslayer_norm_9", &tmva_vars.rmslayer_norm_9);
    tmva_reader->AddVariable("rmslayer_norm_10", &tmva_vars.rmslayer_norm_10);
    tmva_reader->AddVariable("rmslayer_norm_11", &tmva_vars.rmslayer_norm_11);
    tmva_reader->AddVariable("rmslayer_norm_12", &tmva_vars.rmslayer_norm_12);
    tmva_reader->AddVariable("rmslayer_norm_13", &tmva_vars.rmslayer_norm_13);
    tmva_reader->AddVariable("rmslayer_norm_14", &tmva_vars.rmslayer_norm_14);

    tmva_reader->AddVariable("fraclayer_norm_1", &tmva_vars.fraclayer_norm_1);
    tmva_reader->AddVariable("fraclayer_norm_2", &tmva_vars.fraclayer_norm_2);
    tmva_reader->AddVariable("fraclayer_norm_3", &tmva_vars.fraclayer_norm_3);
    tmva_reader->AddVariable("fraclayer_norm_4", &tmva_vars.fraclayer_norm_4);
    tmva_reader->AddVariable("fraclayer_norm_5", &tmva_vars.fraclayer_norm_5);
    tmva_reader->AddVariable("fraclayer_norm_6", &tmva_vars.fraclayer_norm_6);
    tmva_reader->AddVariable("fraclayer_norm_7", &tmva_vars.fraclayer_norm_7);
    tmva_reader->AddVariable("fraclayer_norm_8", &tmva_vars.fraclayer_norm_8);
    tmva_reader->AddVariable("fraclayer_norm_9", &tmva_vars.fraclayer_norm_9);
    tmva_reader->AddVariable("fraclayer_norm_10", &tmva_vars.fraclayer_norm_10);
    tmva_reader->AddVariable("fraclayer_norm_11", &tmva_vars.fraclayer_norm_11);
    tmva_reader->AddVariable("fraclayer_norm_12", &tmva_vars.fraclayer_norm_12);
    tmva_reader->AddVariable("fraclayer_norm_13", &tmva_vars.fraclayer_norm_13);
    tmva_reader->AddVariable("fraclayer_norm_14", &tmva_vars.fraclayer_norm_14);

    tmva_reader->AddVariable("sumrms_norm", &tmva_vars.sumrms_norm);
    tmva_reader->AddVariable("fraclastlayer_norm", &tmva_vars.fraclastlayer_norm);
    tmva_reader->AddVariable("xtrl_norm", &tmva_vars.xtrl_norm);
    tmva_reader->AddSpectator("xtrl", &tmva_vars.xtrl);
}

inline void linkTreeVariables(std::shared_ptr<TChain> evtch, data_vars &vars) {
    evtch->SetBranchAddress("rmslayer_norm_1", &vars.rmslayer_norm_1);
    evtch->SetBranchAddress("rmslayer_norm_2", &vars.rmslayer_norm_2);
    evtch->SetBranchAddress("rmslayer_norm_3", &vars.rmslayer_norm_3);
    evtch->SetBranchAddress("rmslayer_norm_4", &vars.rmslayer_norm_4);
    evtch->SetBranchAddress("rmslayer_norm_5", &vars.rmslayer_norm_5);
    evtch->SetBranchAddress("rmslayer_norm_6", &vars.rmslayer_norm_6);
    evtch->SetBranchAddress("rmslayer_norm_7", &vars.rmslayer_norm_7);
    evtch->SetBranchAddress("rmslayer_norm_8", &vars.rmslayer_norm_8);
    evtch->SetBranchAddress("rmslayer_norm_9", &vars.rmslayer_norm_9);
    evtch->SetBranchAddress("rmslayer_norm_10", &vars.rmslayer_norm_10);
    evtch->SetBranchAddress("rmslayer_norm_11", &vars.rmslayer_norm_11);
    evtch->SetBranchAddress("rmslayer_norm_12", &vars.rmslayer_norm_12);
    evtch->SetBranchAddress("rmslayer_norm_13", &vars.rmslayer_norm_13);
    evtch->SetBranchAddress("rmslayer_norm_14", &vars.rmslayer_norm_14);

    evtch->SetBranchAddress("fraclayer_norm_1", &vars.fraclayer_norm_1);
    evtch->SetBranchAddress("fraclayer_norm_2", &vars.fraclayer_norm_2);
    evtch->SetBranchAddress("fraclayer_norm_3", &vars.fraclayer_norm_3);
    evtch->SetBranchAddress("fraclayer_norm_4", &vars.fraclayer_norm_4);
    evtch->SetBranchAddress("fraclayer_norm_5", &vars.fraclayer_norm_5);
    evtch->SetBranchAddress("fraclayer_norm_6", &vars.fraclayer_norm_6);
    evtch->SetBranchAddress("fraclayer_norm_7", &vars.fraclayer_norm_7);
    evtch->SetBranchAddress("fraclayer_norm_8", &vars.fraclayer_norm_8);
    evtch->SetBranchAddress("fraclayer_norm_9", &vars.fraclayer_norm_9);
    evtch->SetBranchAddress("fraclayer_norm_10", &vars.fraclayer_norm_10);
    evtch->SetBranchAddress("fraclayer_norm_11", &vars.fraclayer_norm_11);
    evtch->SetBranchAddress("fraclayer_norm_12", &vars.fraclayer_norm_12);
    evtch->SetBranchAddress("fraclayer_norm_13", &vars.fraclayer_norm_13);
    evtch->SetBranchAddress("fraclayer_norm_14", &vars.fraclayer_norm_14);

    evtch->SetBranchAddress("sumrms_norm", &vars.sumrms_norm);
    evtch->SetBranchAddress("fraclastlayer_norm", &vars.fraclastlayer_norm);
    evtch->SetBranchAddress("xtrl_norm", &vars.xtrl_norm);
    evtch->SetBranchAddress("xtrl", &vars.xtrl);

    evtch->SetBranchAddress("energy_corr", &vars.evt_corr_energy);

    evtch->SetBranchAddress("energy_bin", &vars.energy_bin);

}

inline void sync_vars(const data_vars &vars, bdt_vars &tmva_vars) {
    tmva_vars.rmslayer_norm_1 = static_cast<float>(vars.rmslayer_norm_1);
    tmva_vars.rmslayer_norm_2 = static_cast<float>(vars.rmslayer_norm_2);
    tmva_vars.rmslayer_norm_3 = static_cast<float>(vars.rmslayer_norm_3);
    tmva_vars.rmslayer_norm_4 = static_cast<float>(vars.rmslayer_norm_4);
    tmva_vars.rmslayer_norm_5 = static_cast<float>(vars.rmslayer_norm_5);
    tmva_vars.rmslayer_norm_6 = static_cast<float>(vars.rmslayer_norm_6);
    tmva_vars.rmslayer_norm_7 = static_cast<float>(vars.rmslayer_norm_7);
    tmva_vars.rmslayer_norm_8 = static_cast<float>(vars.rmslayer_norm_8);
    tmva_vars.rmslayer_norm_9 = static_cast<float>(vars.rmslayer_norm_9);
    tmva_vars.rmslayer_norm_10 = static_cast<float>(vars.rmslayer_norm_10);
    tmva_vars.rmslayer_norm_11 = static_cast<float>(vars.rmslayer_norm_11);
    tmva_vars.rmslayer_norm_12 = static_cast<float>(vars.rmslayer_norm_12);
    tmva_vars.rmslayer_norm_13 = static_cast<float>(vars.rmslayer_norm_13);
    tmva_vars.rmslayer_norm_14 = static_cast<float>(vars.rmslayer_norm_14);

    tmva_vars.fraclayer_norm_1 = static_cast<float>(vars.fraclayer_norm_1);
    tmva_vars.fraclayer_norm_2 = static_cast<float>(vars.fraclayer_norm_2);
    tmva_vars.fraclayer_norm_3 = static_cast<float>(vars.fraclayer_norm_3);
    tmva_vars.fraclayer_norm_4 = static_cast<float>(vars.fraclayer_norm_4);
    tmva_vars.fraclayer_norm_5 = static_cast<float>(vars.fraclayer_norm_5);
    tmva_vars.fraclayer_norm_6 = static_cast<float>(vars.fraclayer_norm_6);
    tmva_vars.fraclayer_norm_7 = static_cast<float>(vars.fraclayer_norm_7);
    tmva_vars.fraclayer_norm_8 = static_cast<float>(vars.fraclayer_norm_8);
    tmva_vars.fraclayer_norm_9 = static_cast<float>(vars.fraclayer_norm_9);
    tmva_vars.fraclayer_norm_10 = static_cast<float>(vars.fraclayer_norm_10);
    tmva_vars.fraclayer_norm_11 = static_cast<float>(vars.fraclayer_norm_11);
    tmva_vars.fraclayer_norm_12 = static_cast<float>(vars.fraclayer_norm_12);
    tmva_vars.fraclayer_norm_13 = static_cast<float>(vars.fraclayer_norm_13);
    tmva_vars.fraclayer_norm_14 = static_cast<float>(vars.fraclayer_norm_14);

    tmva_vars.sumrms_norm = static_cast<float>(vars.sumrms_norm);
    tmva_vars.fraclastlayer_norm = static_cast<float>(vars.fraclastlayer_norm);
    tmva_vars.xtrl_norm = static_cast<float>(vars.xtrl_norm);
    tmva_vars.xtrl = static_cast<float>(vars.xtrl);
}

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
        std::cout << "\nFound acceptance TH1D in input file [" << myhisto->GetName() << "] --> [" << acceptance_file << "]";

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

    // Declare BDT variables
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

    
    size_t n_cuts {1000};
    const double bdt_cut_le {-0.6};
    const double bdt_cut_he {0.6};

    auto get_bdt_cuts = [bdt_cut_le, bdt_cut_he](const size_t size) -> std::vector<std::tuple<double, unsigned int>> {
        std::vector<std::tuple<double, unsigned int>> bdt_cuts (size);
        for (unsigned int idx = 0; idx < size; ++idx) {
            double bdt_cut = bdt_cut_le + (bdt_cut_he - bdt_cut_le) * idx / (size - 1);
            bdt_cuts[idx] = std::make_tuple(bdt_cut, 0);
        }
        return bdt_cuts;
    };

    auto my_cuts = get_bdt_cuts(n_cuts);

    double tmva_classifier {-100};

    auto in_energy_bin = [&energy_binning] (const double energy_gev, const unsigned int energy_bin) -> bool {
        bool in_bin {false};
        if (energy_gev>=energy_binning[energy_bin-1] && energy_gev<energy_binning[energy_bin]) {
            in_bin = true;
        }
        return in_bin;
    };

    for (unsigned int evidx=0; evidx<list_parser->GetEvtTree()->GetEntries(); ++evidx) {

        // Parse event
        list_parser->GetEvtTree()->GetEntry(evidx);
        sync_vars(vars, tmva_vars);
        tmva_classifier = -999;

        if (input_args.verbose)
            evstatus_printer(evidx, kstep, total_events);

        // Filter the event to be contained in the selected energy bin
        if (in_energy_bin(vars.evt_corr_energy*gev, input_args.energy_bin)) {

            if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100) {
                tmva_classifier = tmva_LE_reader->EvaluateMVA(input_args.learning_method.c_str());
            }
            else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<1000) {
                tmva_classifier = tmva_ME_reader->EvaluateMVA(input_args.learning_method.c_str());
            }
            else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<=10000) {
                tmva_classifier = tmva_HE_reader->EvaluateMVA(input_args.learning_method.c_str());
            }

            // Loop over the cuts tuples
            for (auto &cut : my_cuts) {
                if (tmva_classifier > std::get<0>(cut)) {
                    std::get<1>(cut) += 1;
                }
            }
        }
    }

    // Build the fluxes and cuts vectors
    std::vector<double> fluxes (n_cuts), fluxes_E3 (n_cuts), bdt_cuts (n_cuts);
    for (size_t idx = 0; idx < n_cuts; ++idx) {
        fluxes[idx] = std::get<1>(my_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width);
        fluxes_E3[idx] = (std::get<1>(my_cuts[idx])/(input_args.exposure*acceptance*energy_bin_width))*pow(energy_wtsydp, 3);
        bdt_cuts[idx] = std::get<0>(my_cuts[idx]);
    }

    TGraph gr_flux(n_cuts, &bdt_cuts[0], &fluxes[0]);
    gr_flux.SetName((std::string("gr_flux_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux.SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux.GetXaxis()->SetTitle("BDT");
    gr_flux.GetYaxis()->SetTitle("Flux [s^{-1}E^{-1}st^{-1}]");

    TGraph gr_flux_E3(n_cuts, &bdt_cuts[0], &fluxes_E3[0]);
    gr_flux_E3.SetName((std::string("gr_flux_E3_energybin_") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3.SetTitle((std::string("All-Electron Flux vs. BDT cut for energy bin ") + std::to_string(input_args.energy_bin)).c_str());
    gr_flux_E3.GetXaxis()->SetTitle("BDT");
    gr_flux_E3.GetYaxis()->SetTitle("Flux E^{3}*[s^{-1}E^{-1}st^{-1}]");

    TFile *outfile = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (outfile->IsZombie()) {
        std::cerr << "Error writing output file " << input_args.output_path << std::endl;
        exit(1);
    }

    gr_flux.Write();
    gr_flux_E3.Write();

    outfile->Close();

}