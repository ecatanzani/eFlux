#ifndef BDT_SUPPORT_H
#define BDT_SUPPORT_H

#include <string>
#include <memory>

#include "vars.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

struct bdt_weights {
    std::string le_weights;
    std::string me_weights;
    std::string he_weights;
};

struct classifier_vars {
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
};

extern void check_bdt_learnign_method(const std::string method);
extern std::string parse_config_file(std::string bdt_config_file);
extern const bdt_weights get_config_info(std::string parsed_config);
extern void link_reader_vars(
    std::shared_ptr<TMVA::Reader> reader,
    classifier_vars &tmva_vars);
extern void bookMVA(
    std::shared_ptr<TMVA::Reader> reader, 
    const std::string weights,
    const std::string method);

extern void sync_vars(
    const bdt_vars &transformed_vars, 
    classifier_vars &tmva_vars);

#endif