#include "list_parser.h"
#include "bdtutils.h"
#include "bdtinfo.h"
#include "config.h"

#include <memory>
#include <vector>

#include "TFile.h"
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

    evtch->SetBranchAddress("fraclayer_norm_1", &vars.rmslayer_norm_1);
    evtch->SetBranchAddress("fraclayer_norm_2", &vars.rmslayer_norm_2);
    evtch->SetBranchAddress("fraclayer_norm_3", &vars.rmslayer_norm_3);
    evtch->SetBranchAddress("fraclayer_norm_4", &vars.rmslayer_norm_4);
    evtch->SetBranchAddress("fraclayer_norm_5", &vars.rmslayer_norm_5);
    evtch->SetBranchAddress("fraclayer_norm_6", &vars.rmslayer_norm_6);
    evtch->SetBranchAddress("fraclayer_norm_7", &vars.rmslayer_norm_7);
    evtch->SetBranchAddress("fraclayer_norm_8", &vars.rmslayer_norm_8);
    evtch->SetBranchAddress("fraclayer_norm_9", &vars.rmslayer_norm_9);
    evtch->SetBranchAddress("fraclayer_norm_10", &vars.rmslayer_norm_10);
    evtch->SetBranchAddress("fraclayer_norm_11", &vars.rmslayer_norm_11);
    evtch->SetBranchAddress("fraclayer_norm_12", &vars.rmslayer_norm_12);
    evtch->SetBranchAddress("fraclayer_norm_13", &vars.rmslayer_norm_13);
    evtch->SetBranchAddress("fraclayer_norm_14", &vars.rmslayer_norm_14);

    evtch->SetBranchAddress("sumrms_norm", &vars.sumrms_norm);
    evtch->SetBranchAddress("fraclastlayer_norm", &vars.fraclastlayer_norm);
    evtch->SetBranchAddress("xtrl_norm", &vars.xtrl_norm);
    evtch->SetBranchAddress("xtrl", &vars.xtrl);

    evtch->SetBranchAddress("energy_corr", &vars.evt_corr_energy);

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

void ExtractBDTInfo(in_args input_args)
{
    std::shared_ptr<parser> list_parser = std::make_shared<parser>(input_args.input_list, input_args.verbose);
    std::shared_ptr<config> sw_config = std::make_shared<config>(input_args.config_dir);

    if (input_args.verbose) {
        sw_config->PrintActiveFilters();
        sw_config->PrintWeights();
        std::cout << "\n\nTotal number of events: " << list_parser->GetEvtTree()->GetEntries();
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

    tmva_LE_reader->BookMVA(input_args.learning_method.c_str(), (sw_config->GetLEWeights()).c_str());
    tmva_ME_reader->BookMVA(input_args.learning_method.c_str(), (sw_config->GetMEWeights()).c_str());
    tmva_HE_reader->BookMVA(input_args.learning_method.c_str(), (sw_config->GetHEWeights()).c_str());
    
    linkTreeVariables(list_parser->GetEvtTree(), vars);

    // Clone TChain
    std::shared_ptr<TTree> electron_tree = std::shared_ptr<TTree>(static_cast<TTree*>(list_parser->GetEvtTree()->CloneTree(0)));

    // Loop on the events
    double tmva_classifier;
    double gev {0.001};
    bool is_electron;
    for (unsigned int evidx=0; evidx<list_parser->GetEvtTree()->GetEntries(); ++evidx) {
        list_parser->GetEvtTree()->GetEntry(evidx);
        sync_vars(vars, tmva_vars);
        is_electron = false;

        if (vars.evt_corr_energy*gev>=10 && vars.evt_corr_energy*gev<100) {
            tmva_classifier = tmva_LE_reader->EvaluateMVA(input_args.learning_method.c_str());
            if (tmva_classifier>sw_config->GetLEClassifierCut())
                is_electron = true;
        }
        else if (vars.evt_corr_energy*gev>=100 && vars.evt_corr_energy*gev<1000) {
            tmva_classifier = tmva_ME_reader->EvaluateMVA(input_args.learning_method.c_str());
            if (tmva_classifier>sw_config->GetMEClassifierCut())
                is_electron = true;
        }
        else if (vars.evt_corr_energy*gev>=1000 && vars.evt_corr_energy*gev<=10000) {
            tmva_classifier = tmva_HE_reader->EvaluateMVA(input_args.learning_method.c_str());
            if (tmva_classifier>sw_config->GetHEClassifierCut())
                is_electron = true;
        }
        
        if (is_electron)
            electron_tree->Fill();
    }

    TFile* output_file = TFile::Open(input_args.output_path.c_str(), "RECREATE");
    if (output_file->IsZombie()) {
        std::cout << "Error creating output file: [" << input_args.output_path << "]\n\n";
        exit(100);
    }

    electron_tree->Write();
    output_file->Close();

    if (input_args.verbose)
        std::cout << "\nOutput TFile has been written [" << input_args.output_path << "]\n";
}