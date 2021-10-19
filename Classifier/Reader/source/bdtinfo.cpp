#include "list_parser.h"
#include "bdtutils.h"
#include "bdtinfo.h"
#include "config.h"

#include <memory>
#include <vector>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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

    float rmslayer_norm_1 {0};
    float rmslayer_norm_2 {0};
    float rmslayer_norm_3 {0};
    float rmslayer_norm_4 {0};
    float rmslayer_norm_5 {0};
    float rmslayer_norm_6 {0};
    float rmslayer_norm_7 {0};
    float rmslayer_norm_8 {0};
    float rmslayer_norm_9 {0};
    float rmslayer_norm_10 {0};
    float rmslayer_norm_11 {0};
    float rmslayer_norm_12 {0};
    float rmslayer_norm_13 {0};
    float rmslayer_norm_14 {0};
    float fraclayer_norm_1 {0};
    float fraclayer_norm_2 {0};
    float fraclayer_norm_3 {0};
    float fraclayer_norm_4 {0};
    float fraclayer_norm_5 {0};
    float fraclayer_norm_6 {0};
    float fraclayer_norm_7 {0};
    float fraclayer_norm_8 {0};
    float fraclayer_norm_9 {0};
    float fraclayer_norm_10 {0};
    float fraclayer_norm_11 {0};
    float fraclayer_norm_12 {0};
    float fraclayer_norm_13 {0};
    float fraclayer_norm_14 {0};
    float sumrms_norm {0};
    float fraclastlayer_norm {0};
    float xtrl_norm {0};
    float xtrl {0};

    tmva_LE_reader->AddVariable("rmslayer_norm_1", &rmslayer_norm_1);
    tmva_LE_reader->AddVariable("rmslayer_norm_2", &rmslayer_norm_2);
    tmva_LE_reader->AddVariable("rmslayer_norm_3", &rmslayer_norm_3);
    tmva_LE_reader->AddVariable("rmslayer_norm_4", &rmslayer_norm_4);
    tmva_LE_reader->AddVariable("rmslayer_norm_5", &rmslayer_norm_5);
    tmva_LE_reader->AddVariable("rmslayer_norm_6", &rmslayer_norm_6);
    tmva_LE_reader->AddVariable("rmslayer_norm_7", &rmslayer_norm_7);
    tmva_LE_reader->AddVariable("rmslayer_norm_8", &rmslayer_norm_8);
    tmva_LE_reader->AddVariable("rmslayer_norm_9", &rmslayer_norm_9);
    tmva_LE_reader->AddVariable("rmslayer_norm_10", &rmslayer_norm_10);
    tmva_LE_reader->AddVariable("rmslayer_norm_11", &rmslayer_norm_11);
    tmva_LE_reader->AddVariable("rmslayer_norm_12", &rmslayer_norm_12);
    tmva_LE_reader->AddVariable("rmslayer_norm_13", &rmslayer_norm_13);
    tmva_LE_reader->AddVariable("rmslayer_norm_14", &rmslayer_norm_14);

    tmva_LE_reader->AddVariable("fraclayer_norm_1", &fraclayer_norm_1);
    tmva_LE_reader->AddVariable("fraclayer_norm_2", &fraclayer_norm_2);
    tmva_LE_reader->AddVariable("fraclayer_norm_3", &fraclayer_norm_3);
    tmva_LE_reader->AddVariable("fraclayer_norm_4", &fraclayer_norm_4);
    tmva_LE_reader->AddVariable("fraclayer_norm_5", &fraclayer_norm_5);
    tmva_LE_reader->AddVariable("fraclayer_norm_6", &fraclayer_norm_6);
    tmva_LE_reader->AddVariable("fraclayer_norm_7", &fraclayer_norm_7);
    tmva_LE_reader->AddVariable("fraclayer_norm_8", &fraclayer_norm_8);
    tmva_LE_reader->AddVariable("fraclayer_norm_9", &fraclayer_norm_9);
    tmva_LE_reader->AddVariable("fraclayer_norm_10", &fraclayer_norm_10);
    tmva_LE_reader->AddVariable("fraclayer_norm_11", &fraclayer_norm_11);
    tmva_LE_reader->AddVariable("fraclayer_norm_12", &fraclayer_norm_12);
    tmva_LE_reader->AddVariable("fraclayer_norm_13", &fraclayer_norm_13);
    tmva_LE_reader->AddVariable("fraclayer_norm_14", &fraclayer_norm_14);

    tmva_LE_reader->AddVariable("sumrms_norm", &sumrms_norm);
    tmva_LE_reader->AddVariable("fraclastlayer_norm", &fraclastlayer_norm);
    tmva_LE_reader->AddVariable("xtrl_norm", &xtrl_norm);
    tmva_LE_reader->AddVariable("xtrl", &xtrl);

    tmva_ME_reader->AddVariable("rmslayer_norm_1", &rmslayer_norm_1);
    tmva_ME_reader->AddVariable("rmslayer_norm_2", &rmslayer_norm_2);
    tmva_ME_reader->AddVariable("rmslayer_norm_3", &rmslayer_norm_3);
    tmva_ME_reader->AddVariable("rmslayer_norm_4", &rmslayer_norm_4);
    tmva_ME_reader->AddVariable("rmslayer_norm_5", &rmslayer_norm_5);
    tmva_ME_reader->AddVariable("rmslayer_norm_6", &rmslayer_norm_6);
    tmva_ME_reader->AddVariable("rmslayer_norm_7", &rmslayer_norm_7);
    tmva_ME_reader->AddVariable("rmslayer_norm_8", &rmslayer_norm_8);
    tmva_ME_reader->AddVariable("rmslayer_norm_9", &rmslayer_norm_9);
    tmva_ME_reader->AddVariable("rmslayer_norm_10", &rmslayer_norm_10);
    tmva_ME_reader->AddVariable("rmslayer_norm_11", &rmslayer_norm_11);
    tmva_ME_reader->AddVariable("rmslayer_norm_12", &rmslayer_norm_12);
    tmva_ME_reader->AddVariable("rmslayer_norm_13", &rmslayer_norm_13);
    tmva_ME_reader->AddVariable("rmslayer_norm_14", &rmslayer_norm_14);

    tmva_ME_reader->AddVariable("fraclayer_norm_1", &fraclayer_norm_1);
    tmva_ME_reader->AddVariable("fraclayer_norm_2", &fraclayer_norm_2);
    tmva_ME_reader->AddVariable("fraclayer_norm_3", &fraclayer_norm_3);
    tmva_ME_reader->AddVariable("fraclayer_norm_4", &fraclayer_norm_4);
    tmva_ME_reader->AddVariable("fraclayer_norm_5", &fraclayer_norm_5);
    tmva_ME_reader->AddVariable("fraclayer_norm_6", &fraclayer_norm_6);
    tmva_ME_reader->AddVariable("fraclayer_norm_7", &fraclayer_norm_7);
    tmva_ME_reader->AddVariable("fraclayer_norm_8", &fraclayer_norm_8);
    tmva_ME_reader->AddVariable("fraclayer_norm_9", &fraclayer_norm_9);
    tmva_ME_reader->AddVariable("fraclayer_norm_10", &fraclayer_norm_10);
    tmva_ME_reader->AddVariable("fraclayer_norm_11", &fraclayer_norm_11);
    tmva_ME_reader->AddVariable("fraclayer_norm_12", &fraclayer_norm_12);
    tmva_ME_reader->AddVariable("fraclayer_norm_13", &fraclayer_norm_13);
    tmva_ME_reader->AddVariable("fraclayer_norm_14", &fraclayer_norm_14);

    tmva_ME_reader->AddVariable("sumrms_norm", &sumrms_norm);
    tmva_ME_reader->AddVariable("fraclastlayer_norm", &fraclastlayer_norm);
    tmva_ME_reader->AddVariable("xtrl_norm", &xtrl_norm);
    tmva_ME_reader->AddVariable("xtrl", &xtrl);

    tmva_HE_reader->AddVariable("rmslayer_norm_1", &rmslayer_norm_1);
    tmva_HE_reader->AddVariable("rmslayer_norm_2", &rmslayer_norm_2);
    tmva_HE_reader->AddVariable("rmslayer_norm_3", &rmslayer_norm_3);
    tmva_HE_reader->AddVariable("rmslayer_norm_4", &rmslayer_norm_4);
    tmva_HE_reader->AddVariable("rmslayer_norm_5", &rmslayer_norm_5);
    tmva_HE_reader->AddVariable("rmslayer_norm_6", &rmslayer_norm_6);
    tmva_HE_reader->AddVariable("rmslayer_norm_7", &rmslayer_norm_7);
    tmva_HE_reader->AddVariable("rmslayer_norm_8", &rmslayer_norm_8);
    tmva_HE_reader->AddVariable("rmslayer_norm_9", &rmslayer_norm_9);
    tmva_HE_reader->AddVariable("rmslayer_norm_10", &rmslayer_norm_10);
    tmva_HE_reader->AddVariable("rmslayer_norm_11", &rmslayer_norm_11);
    tmva_HE_reader->AddVariable("rmslayer_norm_12", &rmslayer_norm_12);
    tmva_HE_reader->AddVariable("rmslayer_norm_13", &rmslayer_norm_13);
    tmva_HE_reader->AddVariable("rmslayer_norm_14", &rmslayer_norm_14);

    tmva_HE_reader->AddVariable("fraclayer_norm_1", &fraclayer_norm_1);
    tmva_HE_reader->AddVariable("fraclayer_norm_2", &fraclayer_norm_2);
    tmva_HE_reader->AddVariable("fraclayer_norm_3", &fraclayer_norm_3);
    tmva_HE_reader->AddVariable("fraclayer_norm_4", &fraclayer_norm_4);
    tmva_HE_reader->AddVariable("fraclayer_norm_5", &fraclayer_norm_5);
    tmva_HE_reader->AddVariable("fraclayer_norm_6", &fraclayer_norm_6);
    tmva_HE_reader->AddVariable("fraclayer_norm_7", &fraclayer_norm_7);
    tmva_HE_reader->AddVariable("fraclayer_norm_8", &fraclayer_norm_8);
    tmva_HE_reader->AddVariable("fraclayer_norm_9", &fraclayer_norm_9);
    tmva_HE_reader->AddVariable("fraclayer_norm_10", &fraclayer_norm_10);
    tmva_HE_reader->AddVariable("fraclayer_norm_11", &fraclayer_norm_11);
    tmva_HE_reader->AddVariable("fraclayer_norm_12", &fraclayer_norm_12);
    tmva_HE_reader->AddVariable("fraclayer_norm_13", &fraclayer_norm_13);
    tmva_HE_reader->AddVariable("fraclayer_norm_14", &fraclayer_norm_14);

    tmva_HE_reader->AddVariable("sumrms_norm", &sumrms_norm);
    tmva_HE_reader->AddVariable("fraclastlayer_norm", &fraclastlayer_norm);
    tmva_HE_reader->AddVariable("xtrl_norm", &xtrl_norm);
    tmva_HE_reader->AddVariable("xtrl", &xtrl);

    tmva_LE_reader->BookMVA(input_args.learning_method.c_str(), (sw_config->GetLEWeights()).c_str());
    tmva_ME_reader->BookMVA(input_args.learning_method.c_str(), (sw_config->GetMEWeights()).c_str());
    tmva_HE_reader->BookMVA(input_args.learning_method.c_str(), (sw_config->GetHEWeights()).c_str());

    


}